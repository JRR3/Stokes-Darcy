//------------------------------------------------------
//template <int dim>
//void SD<dim>::stokes_operator (MPI_BlockVector &stokes_target,
                              //const Vector<double> &lambda_source)
template <int dim>
void SD<dim>::stokes_operator ()
{
  TimerOutput::Scope timer_section(computing_timer, "Stokes operator");
  pcout << ">>>Stokes operator" << std::endl;

  stokes_rhs=0;

  const unsigned int   dofs_per_cell   = stokes_fe.dofs_per_cell;
  const unsigned int   faces_per_cell  = GeometryInfo<dim>::faces_per_cell;

  if (rebuild_stokes_matrix)
    stokes_matrix=0;

  stokes_rhs=0;

  FEValues<dim>     stokes_fe_values (stokes_fe, q_rule,
                                      update_values    |
                                      update_quadrature_points  |
                                      update_JxW_values |
                                      ((rebuild_stokes_matrix == true)
                                       ?
                                       update_gradients
                                       :
                                       UpdateFlags(0)));

  FEFaceValues<dim> stokes_fe_face_values (stokes_fe, face_q_rule,
                           update_normal_vectors | 
                           update_values    |
                           update_quadrature_points  |
                           update_JxW_values);


  FEValues<dim-1,dim>   lambda_fe_values (lambda_fe, face_q_rule, update_values);   

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  std::vector<double>       lambda_face_values (n_face_q_points);
  FullMatrix<double>        local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>            local_rhs (dofs_per_cell);

  stokes_exact_solution.set_time(time);
  std::vector<Tensor<1,dim> >       rhs_values (n_q_points);
  std::vector<double>               div_values (n_q_points);

  std::vector<Tensor<1,dim> >          phi_u       (dofs_per_cell);
  std::vector<SymmetricTensor<2,dim> > grads_phi_u (dofs_per_cell);
  std::vector<double>                  div_phi_u   (dofs_per_cell);
  std::vector<double>                  phi_p       (dofs_per_cell);

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  typename DoFHandler<dim>::active_cell_iterator
  cell = stokes_dof_handler.begin_active(),
  endc = stokes_dof_handler.end();

  for (; cell!=endc; ++cell)
  if(cell->is_locally_owned())
    {
      stokes_fe_values.reinit (cell);
      local_matrix = 0;
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

          if (rebuild_stokes_matrix)
          {
            grads_phi_u[k] = stokes_fe_values[velocities].symmetric_gradient(k,q);
            div_phi_u[k]   = stokes_fe_values[velocities].divergence (k, q);
          }
        }

        if (rebuild_stokes_matrix)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
            local_matrix(i,j) += (2  * viscosity * 
                                (grads_phi_u[i] * grads_phi_u[j])
                                - div_phi_u[i] * phi_p[j]
                                - phi_p[i] * div_phi_u[j])
                               * stokes_fe_values.JxW(q);

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
          get_function_values(cell, lambda_vector, lambda_fe_values, lambda_face_values);

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

      if (rebuild_stokes_matrix)
        stokes_constraints.distribute_local_to_global (local_matrix,
                                                       local_rhs,
                                                       local_dof_indices,
                                                       stokes_matrix,
                                                       stokes_rhs);
      else
        stokes_constraints.distribute_local_to_global (local_rhs,
                                                       local_dof_indices,
                                                       stokes_rhs);

    }//end_cell

    if (rebuild_stokes_matrix)
      stokes_matrix.compress(VectorOperation::add);

    stokes_rhs.compress(VectorOperation::add);

    stokes_solver ();

    {
      double my_rhs_l2_norm = stokes_rhs.l2_norm();
      double my_sol_l2_norm = stokes_solution.l2_norm();
      pcout << "Stokes rhs norm: " << std::setprecision(10) << my_rhs_l2_norm << std::endl;
      pcout << "Stokes sol norm: " << std::setprecision(10) << my_sol_l2_norm << std::endl;
    }

}//end_stokes_op
