//------------------------------------------------------------------------
/** The Darcy star operator uses a different matrix.
 * In theory, the Darcy matrix has to be recomputed on every time step
 * since the porosity is changing through time and the mass matrix term
 * (\beta(\eta) U, V)_{\Omega_1} depends on \eta. A sensitivity analysis may 
 * suggest that it is not necessary to recompute the matrix. See (step-43)
 */
//------------------------------------------------------------------------
template <int dim>
void SD<dim>::darcy_star (MPI_BlockVector &darcy_target,
                         const Vector<double> &flux_source)
{

  TimerOutput::Scope timer_section(computing_timer, "Darcy star");
  pcout << ">>>Darcy star " << std::endl;

  darcy_rhs    = 0;
  
  FEFaceValues<dim> darcy_fe_face_values (darcy_fe, face_q_rule,
                                    update_values    | update_normal_vectors |
                                    update_quadrature_points  | update_JxW_values);

  //------------------------FLUX
  FEValues<dim-1,dim> flux_fe_values (flux_fe, face_q_rule, update_values |
                                                            update_quadrature_points |
                                                            update_JxW_values);

  const unsigned int   dofs_per_cell   = darcy_fe.dofs_per_cell;

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  std::vector<Point<dim> >             face_q_points(n_face_q_points);

  //std::vector<double> div_values      (n_q_points);
  std::vector<double> flux_values (n_face_q_points);

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  typename DoFHandler<dim>::active_cell_iterator
  cell = darcy_dof_handler.begin_active(),
  endc = darcy_dof_handler.end();

  //pcout << "*********** DARCY///////////////// " << std::endl;

  for( auto flux_cell  = flux_dof_handler.begin_active();
            flux_cell != flux_dof_handler.end(); ++flux_cell)
  {
    local_rhs = 0;
    flux_fe_values.reinit(flux_cell);
    darcy_flux_map.
    get_source_basis_functions_values(flux_fe_values,
                                      i_sqrt_flux_cell_measure, 
                                      flux_cell, flux_source, 
                                      local_dof_indices, local_rhs);

    //cell->get_dof_indices (local_dof_indices);

    for(unsigned int i = 0; i < dofs_per_cell; ++i)
      darcy_rhs(local_dof_indices[i]) += local_rhs(i);
  }

  //darcy_matrix = 0;
  //darcy_matrix.copy_from(darcy_template_matrix); 
  //darcy_zero_constraints.condense(darcy_matrix, darcy_rhs);
  //darcy_solver.initialize(darcy_matrix);
  //darcy_target = 0;
  //darcy_solver.vmult(darcy_target, darcy_rhs);
  //darcy_zero_constraints.distribute (darcy_target);

}//end_star
//------------------------------------------------------------------------
/** In this version, we enforce weakly the boundary conditions. This implies
 * that the only kind of boundary conditions that we consider are of the
 * homogeneous type.
 */
//template <int dim>
//void SD<dim>::darcy_operator (MPI_BlockVector &darcy_target,
                             //const Vector<double> &lambda_source)
template <int dim>
void SD<dim>::darcy_operator ()
{
  TimerOutput::Scope timer_section(computing_timer, "Darcy operator");
  pcout << ">>>Darcy operator " << std::endl;

  darcy_rhs = 0;
  
  const unsigned int   dofs_per_cell   = darcy_fe.dofs_per_cell;
  const unsigned int   faces_per_cell  = GeometryInfo<dim>::faces_per_cell;

  FEValues<dim>     darcy_fe_values (darcy_fe, q_rule,
                             update_values    | update_gradients |
                             update_quadrature_points  | update_JxW_values);

  FEFaceValues<dim> darcy_fe_face_values (darcy_fe, face_q_rule,
                                    update_values    | update_normal_vectors |
                                    update_quadrature_points  | update_JxW_values);

  FEValues<dim-1,dim>   lambda_fe_values (lambda_fe, face_q_rule, update_values);   

  //----------------------------------------------------------------
  std::vector<double> lambda_face_values (n_face_q_points);
  //std::vector<Point<dim> > face_q_points (n_face_q_points);
  //----------------------------------------------------------------

  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  darcy_exact_solution.set_time(time);

  std::vector<Tensor<1,dim> >   phi_u       (dofs_per_cell);
  std::vector<double>           phi_p       (dofs_per_cell);
  std::vector<double>           div_values  (n_q_points);
  //std::vector<double> div_face_values  (n_face_q_points);
  std::vector<double> p_face_values    (n_face_q_points);
  std::vector<Tensor<1,dim> >   rhs_values  (n_q_points);
  //std::vector<Tensor<1,dim> >     grad_div_values  (n_q_points);

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  typename DoFHandler<dim>::active_cell_iterator
  cell = darcy_dof_handler.begin_active(),
  endc = darcy_dof_handler.end();

  for (; cell!=endc; ++cell)
  if(cell->is_locally_owned())
  {
    darcy_fe_values.reinit (cell);
    local_rhs = 0;

    darcy_exact_solution.get_div_list (darcy_fe_values.get_quadrature_points(),
                                div_values);
    darcy_exact_solution.get_rhs_list (darcy_fe_values.get_quadrature_points(),
                                      rhs_values);

    for (unsigned int q=0; q<n_q_points; ++q)
    {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            phi_u[k]   = darcy_fe_values[velocities].value (k, q);
            phi_p[k]   = darcy_fe_values[pressure].value (k, q);
          }
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {

          local_rhs(i) += ( - phi_p[i] * div_values[q]  
                            + rhs_values[q] * phi_u[i]
                            //- grad_div * (grad_div_values[q] * phi_u[i]) 
                          )
                            *   darcy_fe_values.JxW(q);
        }
    }//end_q_points

    if (cell->at_boundary())
    for (unsigned int face_no=0; face_no < faces_per_cell; ++face_no)
    {
      if (darcy_bc_are_enforced_weakly)
      if (cell->face(face_no)->boundary_id() == wall_id)
      {
        darcy_fe_face_values.reinit (cell, face_no);

        //darcy_exact_solution.get_div_list (darcy_fe_face_values.get_quadrature_points(),
                                    //div_face_values);
        darcy_exact_solution.get_p_list (darcy_fe_face_values.get_quadrature_points(),
                                    p_face_values);
        for (unsigned int q=0; q<n_face_q_points; ++q)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          local_rhs(i) +=  (darcy_fe_face_values[velocities].value (i, q) *
                            darcy_fe_face_values.normal_vector(q) *
                            ( 
                             //div_face_values[q] 
                             - p_face_values[q] 
                            ) * 
                            darcy_fe_face_values.JxW(q));
        }
      }//end_darcy_cell_at_boundary
      if (cell->face(face_no)->boundary_id() == interface_id)
      {
        darcy_fe_face_values.reinit (cell, face_no);
        //face_q_points = darcy_fe_face_values.get_quadrature_points();
        //darcy_lambda_map.
        //get_function_values(cell, lambda_source, face_q_points, lambda_face_values);
        darcy_lambda_map.
        get_function_values(cell, lambda_vector, lambda_fe_values, lambda_face_values);

        for (unsigned int q=0; q<n_face_q_points; ++q)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
         //MINUS REMOVED
          local_rhs(i) +=  (darcy_fe_face_values[velocities].value (i, q) *
                            darcy_fe_face_values.normal_vector(q) *
                            lambda_face_values[q] *
                            darcy_fe_face_values.JxW(q));
        }//face_q_points

      }//end_on_the_interface
    }//end_boundary_contribution

    cell->get_dof_indices (local_dof_indices);


    darcy_constraints.distribute_local_to_global (local_rhs,
                                                  local_dof_indices,
                                                  darcy_rhs);
  }//end_cell

    darcy_rhs.compress (VectorOperation::add);
    darcy_fgmres_solver ();

    {
      double my_rhs_l2_norm = darcy_rhs.l2_norm();
      pcout << "darcy rhs norm: " << std::setprecision(10) << my_rhs_l2_norm << std::endl;
      double my_sol_l2_norm = darcy_solution.l2_norm();
      pcout << "darcy sol norm: " << std::setprecision(10) << my_sol_l2_norm << std::endl;
    }


}//end_operator 
//------------------------------------------------------------------------
/** See the comment for the Darcy star operator. Homogeneous boundary conditions
 * are required in this case.
 */
template <int dim>
void SD<dim>::darcy_prime (MPI_BlockVector &darcy_target,
                             const Vector<double> &lambda_source)
{
  TimerOutput::Scope timer_section(computing_timer, "Darcy prime");
  pcout << ">>>Darcy prime " << std::endl;

  darcy_rhs    = 0;
  
  //QGauss<dim>   q_rule(darcy_degree+3);
  //QGauss<dim-1> face_q_rule(darcy_degree+3);
  const unsigned int   dofs_per_cell   = darcy_fe.dofs_per_cell;
  //const unsigned int   n_q_points      = q_rule.size();
  //const unsigned int   n_face_q_points = face_q_rule.size();
  const unsigned int   faces_per_cell  = GeometryInfo<dim>::faces_per_cell;

  //FEValues<dim> darcy_fe_values (darcy_fe, q_rule,
                           //update_values    | update_gradients |
                           //update_quadrature_points  | update_JxW_values);

  FEFaceValues<dim> darcy_fe_face_values (darcy_fe, face_q_rule,
                                    update_values    | update_normal_vectors |
                                    update_quadrature_points  | update_JxW_values);

  //FEValues<dim-1,dim> lambda_fe_values (lambda_fe, face_q_rule, update_values);
                           //update_JxW_values);
  //Functions::FEFieldFunction<dim, DoFHandler<dim-1,dim>, Vector<double> > 
   //lambda_field (lambda_dof_handler, lambda_source);
  std::vector<double> lambda_face_values (n_face_q_points);
  std::vector<Point<dim> > face_q_points (n_face_q_points);


  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  // With all this in place, we can go on with the loop over all cells. The
  // body of this loop has been discussed in the introduction, and will not
  // be commented any further here:
  typename DoFHandler<dim>::active_cell_iterator
  cell = darcy_dof_handler.begin_active(),
  endc = darcy_dof_handler.end();

  //typename DoFHandler<dim-1,dim>::active_cell_iterator lambda_cell;

  for (; cell!=endc; ++cell)
    {
      local_rhs = 0;

      if (cell->at_boundary())
      for (unsigned int face_no=0; face_no < faces_per_cell; ++face_no)
      {

        if (cell->face(face_no)->boundary_id() == wall_id)
          {
          }
        else if (cell->face(face_no)->boundary_id() == interface_id)
          {
            darcy_fe_face_values.reinit (cell, face_no);


            face_q_points = darcy_fe_face_values.get_quadrature_points();
            darcy_lambda_map.
            get_function_values(cell, lambda_source, face_q_points, lambda_face_values);

            for (unsigned int q=0; q<n_face_q_points; ++q)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
               //MINUS REMOVED
                local_rhs(i) +=  (darcy_fe_face_values[velocities].value (i, q) *
                                  darcy_fe_face_values.normal_vector(q) *
                                  lambda_face_values[q] *
                                  darcy_fe_face_values.JxW(q));
              }

          }
      }//end_boundary_contribution

      cell->get_dof_indices (local_dof_indices);
      for(unsigned int i = 0; i < dofs_per_cell; ++i)
        darcy_rhs(local_dof_indices[i]) += local_rhs(i);

      //darcy_constraints.distribute_local_to_global (local_matrix, local_rhs,
                                              //local_dof_indices,
                                              //darcy_matrix, darcy_rhs);
    }//end_cell


    //darcy_matrix.copy_from(darcy_template_matrix); 
    //darcy_zero_constraints.condense(darcy_matrix, darcy_rhs);
    //darcy_solver.initialize(darcy_matrix);
    //darcy_target = 0;
    //darcy_solver.vmult(darcy_target, darcy_rhs);
    //darcy_zero_constraints.distribute (darcy_target);

}//end_prime
