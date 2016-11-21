/** The Stokes star operator corresponds to the (N')* operator
 * that arises in the article:
 * Approximation of the Stokes–Darcy System by Optimization
 * (see page: 10)
 * This operator uses homogeneous boundary conditions
 * Use stokes_zero_constraints and condense is one of the options
 * CHECK: An alternative is to use a direct incorporation into the matrix,
 * but this would require the use of two matrices. In parallel this
 * is an acceptable price. The original version uses condense, but
 * we opt in this case to use additional matrices, expecting that
 * in a time dependent probelm the cost of the CG-LSQ algorithm is not
 * high.
 */
template <int dim>

void SD<dim>::stokes_star (MPI_BlockVector &stokes_target,
                          const Vector<double> &flux_source)
{
  TimerOutput::Scope timer_section(computing_timer, "Stokes star");
  pcout << ">>>Stokes star" << std::endl;

  stokes_rhs=0;



  //------------------------FLUX
  FEValues<dim-1,dim> flux_fe_values (flux_fe, face_q_rule, update_values |
                                                            update_quadrature_points |
                                                            update_JxW_values);

  const unsigned int   dofs_per_cell   = stokes_fe.dofs_per_cell;
  //const unsigned int   faces_per_cell  = GeometryInfo<dim>::faces_per_cell;

  //const unsigned int   n_q_points      = q_rule.size();
  //const unsigned int   n_face_q_points      = face_q_rule.size();

  //FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  //const RightHandSide<dim>          right_hand_side;
  stokes_exact_solution.set_time(time);
  ////----make general
  //Tensor<1,dim>                     tangent_vector;
  //tangent_vector[0] = 1;
  //tangent_vector[1] = 0;
  //----
  std::vector<double>               flux_values (n_face_q_points);
  std::vector<Point<dim> >          face_q_points (n_face_q_points);
                                                

  // Next, we need two objects that work as extractors for the FEValues
  // object. Their use is explained in detail in the report on @ref
  // vector_valued :
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  //std::vector<SymmetricTensor<2,dim> > symgrad_phi_u (dofs_per_cell);
  //std::vector<Tensor<2,dim> >          grad_phi_u  (dofs_per_cell);
  std::vector<Tensor<1,dim> >          phi_u       (dofs_per_cell);
  //std::vector<double>                  div_phi_u   (dofs_per_cell);
  //std::vector<double>                  phi_p       (dofs_per_cell);
  //pcout << "Before exit" << std::endl;
  //exit(0);
  //pcout << "After exit" << std::endl;

  //pcout << "*********** STOKES///////////////// " << std::endl;

  for( auto flux_cell  = flux_dof_handler.begin_active();
            flux_cell != flux_dof_handler.end(); ++flux_cell)
  {
    local_rhs = 0;
    flux_fe_values.reinit(flux_cell);
    stokes_flux_map.
    get_source_basis_functions_values(flux_fe_values,
                                      i_sqrt_flux_cell_measure, 
                                      flux_cell, flux_source, 
                                      local_dof_indices, local_rhs);


    for(unsigned int i = 0; i < dofs_per_cell; ++i)
      stokes_rhs(local_dof_indices[i]) += local_rhs(i);
  }

  //Solve stokes matrix Ax = b, where
  //A = stokes_matrix
  //x = stokes_target ---> stokes_vector
  //b = stokes_rhs

  //stokes_matrix.copy_from(stokes_template_matrix); 
  //stokes_zero_constraints.condense(stokes_matrix, stokes_rhs);
  //stokes_solver.initialize(stokes_matrix);
  //stokes_target = 0;
  //stokes_solver.vmult(stokes_target, stokes_rhs);
  //stokes_zero_constraints.distribute (stokes_target);

}//end_stokes_star
//------------------------------------------------------
template <int dim>
void SD<dim>::stokes_operator (MPI_BlockVector &stokes_target,
                              const Vector<double> &lambda_source)
{
  TimerOutput::Scope timer_section(computing_timer, "Stokes operator");
  pcout << ">>>Stokes operator" << std::endl;

  stokes_rhs=0;

  const unsigned int   dofs_per_cell   = stokes_fe.dofs_per_cell;
  const unsigned int   faces_per_cell  = GeometryInfo<dim>::faces_per_cell;

  FEFaceValues<dim> stokes_fe_face_values (stokes_fe, face_q_rule,
                           update_normal_vectors | 
                           update_values    |
                           update_quadrature_points  |
                           update_JxW_values);

  FEValues<dim> stokes_fe_values (stokes_fe, q_rule,
                           update_values    |
                           update_quadrature_points  |
                           update_JxW_values);

  FEValues<dim-1,dim>   lambda_fe_values (lambda_fe, face_q_rule, update_values);   

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  std::vector<double>       lambda_face_values (n_face_q_points);
  //std::vector<Point<dim> >  face_q_points(n_face_q_points);
  Vector<double>            local_rhs (dofs_per_cell);

  stokes_exact_solution.set_time(time);
  std::vector<Tensor<1,dim> >       rhs_values (n_q_points);
  std::vector<double>               div_values (n_q_points);

  //----make general
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  std::vector<Tensor<1,dim> >          phi_u       (dofs_per_cell);
  std::vector<double>                  phi_p       (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
  cell = stokes_dof_handler.begin_active(),
  endc = stokes_dof_handler.end();

  typename DoFHandler<dim-1,dim>::active_cell_iterator lambda_cell;

  for (; cell!=endc; ++cell)
  if(cell->is_locally_owned())
    {
      stokes_fe_values.reinit (cell);
      local_rhs = 0;

      stokes_exact_solution.get_rhs_list (stokes_fe_values.get_quadrature_points(),
                                        rhs_values);
      stokes_exact_solution.get_div_list (stokes_fe_values.get_quadrature_points(),
                                        div_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            phi_u[k] = stokes_fe_values[velocities].value (k, q);
            phi_p[k] = stokes_fe_values[pressure].value (k, q);
          }

          for (unsigned int i=0; i<dofs_per_cell; ++i)
              local_rhs(i) += (
                                 phi_u[i] * rhs_values[q] 
                               - phi_p[i] * div_values[q] //Why + does not work?
                              )  
                                * stokes_fe_values.JxW(q);

        }//end_q_points

      //NEW
        if (cell->at_boundary())
        for (unsigned int face_no=0; face_no < faces_per_cell; ++face_no)
        if (cell->face(face_no)->boundary_id() == interface_id)
        {
          stokes_fe_face_values.reinit (cell, face_no);

          //face_q_points = stokes_fe_face_values.get_quadrature_points();
          stokes_lambda_map.
          get_function_values(cell, lambda_source, lambda_fe_values, lambda_face_values);

          for (unsigned int q=0; q<n_face_q_points; ++q)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
           //MINUS REMOVED
            local_rhs(i) += (stokes_fe_face_values[velocities].value (i, q) *
                              stokes_fe_face_values.normal_vector(q) *
                              lambda_face_values[q] *
                              stokes_fe_face_values.JxW(q));
          }

        }//end_interface_contribution


      cell->get_dof_indices (local_dof_indices);

      stokes_constraints.distribute_local_to_global (local_rhs, local_dof_indices, stokes_rhs);
    }//end_cell

    //stokes_matrix = 0;
    //stokes_matrix.copy_from(stokes_template_matrix); 
    //stokes_constraints.condense(stokes_matrix, stokes_rhs);
    //stokes_solver.initialize(stokes_matrix);
    //stokes_target = 0;
    //stokes_solver.vmult(stokes_target, stokes_rhs);
    //stokes_constraints.distribute (stokes_target);
    double my_l2_norm = stokes_rhs.l2_norm();
    pcout << "Stokes norm: " << std::setprecision(10) << my_l2_norm << std::endl;

}//end_stokes_op
//-------------------------------------------------------------------
template <int dim>
void SD<dim>::stokes_prime (MPI_BlockVector &stokes_target,
                           const Vector<double> &lambda_source)
{
  TimerOutput::Scope timer_section(computing_timer, "Stokes prime");
  pcout << ">>>Stokes prime" << std::endl;

  stokes_rhs=0;

  //QGauss<dim>   q_rule(stokes_degree+3);
  //QGauss<dim-1>   face_q_rule(stokes_degree+3);
  const unsigned int   dofs_per_cell   = stokes_fe.dofs_per_cell;
  const unsigned int   faces_per_cell  = GeometryInfo<dim>::faces_per_cell;

  //const unsigned int   n_q_points      = q_rule.size();
  //const unsigned int   n_face_q_points      = face_q_rule.size();

  FEFaceValues<dim> stokes_fe_face_values (stokes_fe, face_q_rule,
                           update_normal_vectors | 
                           update_values    |
                           update_quadrature_points  |
                           update_JxW_values);

  std::vector<double>       lambda_face_values (n_face_q_points);
  std::vector<Point<dim> >  face_q_points      (n_face_q_points);


  Vector<double>            local_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  //----make general

  // Next, we need two objects that work as extractors for the FEValues
  // object. Their use is explained in detail in the report on @ref
  // vector_valued :
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);


  typename DoFHandler<dim>::active_cell_iterator
  cell = stokes_dof_handler.begin_active(),
  endc = stokes_dof_handler.end();

  typename DoFHandler<dim-1,dim>::active_cell_iterator lambda_cell;

  for (; cell!=endc; ++cell)
    {
      //stokes_fe_values.reinit (cell);
      local_rhs = 0;

      //NEW
        if (cell->at_boundary())
        for (unsigned int face_no=0; face_no < faces_per_cell; ++face_no)
        if (cell->face(face_no)->boundary_id() == interface_id)
        {
          stokes_fe_face_values.reinit (cell, face_no);

          face_q_points = stokes_fe_face_values.get_quadrature_points();
          stokes_lambda_map.
          get_function_values(cell, lambda_source, face_q_points, lambda_face_values);

          for (unsigned int q=0; q<n_face_q_points; ++q)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
             //MINUS REMOVED
              local_rhs(i) += (stokes_fe_face_values[velocities].value (i, q) *
                                stokes_fe_face_values.normal_vector(q) *
                                lambda_face_values[q] *
                                stokes_fe_face_values.JxW(q));
            }

        }//end_interface_contribution


      cell->get_dof_indices (local_dof_indices);
      for(unsigned int i = 0; i < dofs_per_cell; ++i)
        stokes_rhs(local_dof_indices[i]) += local_rhs(i);

      //stokes_constraints.distribute_local_to_global (local_matrix, local_rhs,
                                              //local_dof_indices,
                                              //stokes_matrix, stokes_rhs);
    }//end_cell

    //stokes_matrix.copy_from(stokes_template_matrix); 
    //stokes_zero_constraints.condense(stokes_matrix, stokes_rhs);
    //stokes_solver.initialize(stokes_matrix);
    //stokes_target = 0;
    //stokes_solver.vmult(stokes_target, stokes_rhs);
    //stokes_zero_constraints.distribute (stokes_target);

  
}//end_stokes_prime
//---------------------------------------------
