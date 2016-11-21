template <int dim>
void SD<dim>::assemble_darcy_system ()
{

  if(rebuild_darcy_matrix == false)
    return;

  TimerOutput::Scope timer_section(computing_timer, "Darcy assembly");

  darcy_matrix = 0;
  //darcy_rhs    = 0;

  FEValues<dim> darcy_fe_values (darcy_fe, q_rule,
                                 update_values    |
                                 update_gradients |
                                 update_quadrature_points  |
                                 update_JxW_values);

  FEValues<dim> porosity_fe_values (porosity_fe, q_rule,
                                      update_values);

  //FEFaceValues<dim> darcy_fe_face_values (darcy_fe, face_q_rule,
                                          //update_values    | update_normal_vectors |
                                          //update_quadrature_points  | update_JxW_values);

  const unsigned int   dofs_per_cell   = darcy_fe.dofs_per_cell;

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  //Vector<double>       local_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  std::vector<Tensor<2,dim> >       k_inverse_values (n_q_points);

  //std::vector<double>               old_porosity_values (n_q_points);

  std::vector<Tensor<1,dim> >       phi_u (dofs_per_cell);
  std::vector<double>               div_phi_u (dofs_per_cell);
  std::vector<double>               phi_p (dofs_per_cell);

  const FEValuesExtractors::Vector  velocities (0);
  const FEValuesExtractors::Scalar  pressure (dim);

  typename DoFHandler<dim>::active_cell_iterator
  cell = darcy_dof_handler.begin_active(),
  endc = darcy_dof_handler.end();

  typename DoFHandler<dim>::active_cell_iterator
  porosity_cell = porosity_dof_handler.begin_active();

  for (; cell!=endc; ++cell, ++porosity_cell)
  if(cell->is_locally_owned())
  {
    darcy_fe_values.reinit (cell);
    porosity_fe_values.reinit (porosity_cell);

    local_matrix = 0;
    //local_rhs = 0;

    //porosity_fe_values.get_function_values (old_porosity_solution, old_porosity_values);

    //pressure_right_hand_side.value_list (darcy_fe_values.get_quadrature_points(),
                                         //pressure_rhs_values);
    k_inverse.value_list (darcy_fe_values.get_quadrature_points(),
                          k_inverse_values);

    for (unsigned int q=0; q<n_q_points; ++q)
    {
      for (unsigned int k=0; k<dofs_per_cell; ++k)
      {
        phi_u[k]     = darcy_fe_values[velocities].value (k,q);
        div_phi_u[k] = darcy_fe_values[velocities].divergence (k,q);
        phi_p[k]     = darcy_fe_values[pressure].value (k,q);
      }
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
        for (unsigned int j=0; j<=i; ++j)
        {
          local_matrix(i,j) += (phi_u[i] * k_inverse_values[q] * phi_u[j]
                                + grad_div * div_phi_u[i] * div_phi_u[j]
                                - div_phi_u[i] * phi_p[j]
                                - phi_p[i] * div_phi_u[j])
                               * darcy_fe_values.JxW(q);
        }

      }
    }//for_each_q_point

    for (unsigned int i=0; i<dofs_per_cell; ++i)
      for (unsigned int j=i+1; j<dofs_per_cell; ++j)
        local_matrix(i,j) = local_matrix(j,i);

    cell->get_dof_indices (local_dof_indices);

    darcy_constraints.distribute_local_to_global (local_matrix,
                                                  local_dof_indices,
                                                  darcy_matrix);

  }//for_each_cell

  darcy_matrix.compress(VectorOperation::add);

  //double matrix_norm = 0;
  //for(unsigned int i = 0; i < dim; ++i)
  //for(unsigned int j = 0; j < dim; ++j)
    //matrix_norm += darcy_matrix.block(i,j).frobenius_norm();
//
  //pcout << "Darcy matrix norm is: " << matrix_norm << std::endl;
    
}

