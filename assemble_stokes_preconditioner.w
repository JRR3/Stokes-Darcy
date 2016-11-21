/** This function assembles the stokes preconditioner for the
 * matrix with the original boundary conditions and the 
 * matrix that uses the zero boundary conditions.
 */
template <int dim>
void SD<dim>::assemble_stokes_preconditioner ()
{
  //TimerOutput::Scope timer_section(computing_timer, "Stokes preconditioner");
  //This function is called in build_stokes_preconditioner();

  stokes_preconditioner_matrix = 0;

  FEValues<dim>     stokes_fe_values (stokes_fe, q_rule,
                                      update_JxW_values |
                                      update_values |
                                      update_gradients);

  const unsigned int   dofs_per_cell   = stokes_fe.dofs_per_cell;

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  std::vector<Tensor<2,dim> > grad_phi_u (dofs_per_cell);
  std::vector<double>         phi_p      (dofs_per_cell);

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

    for (unsigned int q=0; q<n_q_points; ++q)
    {
      for (unsigned int k=0; k<dofs_per_cell; ++k)
      {
        grad_phi_u[k] = stokes_fe_values[velocities].gradient(k,q);
        phi_p[k]      = stokes_fe_values[pressure].value (k, q);
      }

      for (unsigned int i=0; i<dofs_per_cell; ++i)
      for (unsigned int j=0; j<dofs_per_cell; ++j)
          local_matrix(i,j) += (scalar_product (grad_phi_u[i], grad_phi_u[j])
                                +
                                phi_p[i] * phi_p[j])
                               * stokes_fe_values.JxW(q);
    }

    cell->get_dof_indices (local_dof_indices);

    stokes_constraints.distribute_local_to_global (local_matrix,
                                                   local_dof_indices,
                                                   stokes_preconditioner_matrix);
    //INCLUDED
    stokes_zero_bc_constraints.distribute_local_to_global (local_matrix,
                                                   local_dof_indices,
                                                   stokes_zero_bc_preconditioner_matrix);
  }//end_cell

  
  stokes_preconditioner_matrix.compress(VectorOperation::add);
  stokes_zero_bc_preconditioner_matrix.compress(VectorOperation::add);
}//end_assemble_stokes_preconditioner

