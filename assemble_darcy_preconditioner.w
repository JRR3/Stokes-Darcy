template <int dim>
void SD<dim>::assemble_darcy_preconditioner ()
{
  //TimerOutput::Scope timer_section(computing_timer, "Darcy preconditioner");

  //pcout << "   Rebuilding darcy preconditioner..." << std::endl;

  darcy_preconditioner_matrix = 0;

  FEValues<dim>     darcy_fe_values (darcy_fe, q_rule,
                                     update_JxW_values |
                                     update_values |
                                     update_gradients |
                                     update_quadrature_points);

  const unsigned int   dofs_per_cell   = darcy_fe.dofs_per_cell;

  std::vector<Tensor<2,dim> >       k_inverse_values (n_q_points);


  FullMatrix<double>                    local_matrix (dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index>  local_dof_indices (dofs_per_cell);

  std::vector<Tensor<1,dim> > phi_u   (dofs_per_cell);
  std::vector<Tensor<1,dim> > grad_phi_p (dofs_per_cell);

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  typename DoFHandler<dim>::active_cell_iterator
  cell = darcy_dof_handler.begin_active(),
  endc = darcy_dof_handler.end();

  for (; cell!=endc; ++cell)
  if(cell->is_locally_owned())
  {
    darcy_fe_values.reinit (cell);

    local_matrix = 0;

    k_inverse.value_list (darcy_fe_values.get_quadrature_points(),
                          k_inverse_values);

    for (unsigned int q=0; q<n_q_points; ++q)
      {

        const Tensor<2,dim> permeability     = invert(k_inverse_values[q]);

        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            phi_u[k]       = darcy_fe_values[velocities].value (k,q);
            grad_phi_p[k]  = darcy_fe_values[pressure].gradient (k,q);
          }

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
              local_matrix(i,j) += (k_inverse_values[q] * 
                                    phi_u[i] * phi_u[j]
                                    +
                                    permeability * 
                                    grad_phi_p[i] * grad_phi_p[j])
                                   * darcy_fe_values.JxW(q);
            }
      }

    cell->get_dof_indices (local_dof_indices);
    darcy_preconditioner_constraints.distribute_local_to_global (local_matrix,
                                                                 local_dof_indices,
                                                                 darcy_preconditioner_matrix);
  }//end_cell
  darcy_preconditioner_matrix.compress(VectorOperation::add);

  double my_f_norm = darcy_preconditioner_matrix.block(0,0).frobenius_norm();
  my_f_norm += darcy_preconditioner_matrix.block(1,1).frobenius_norm();
  pcout << "The F norm of darcy prec is: " << my_f_norm << std::endl;
}
