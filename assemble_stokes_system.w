template <int dim>
void StokesProblem<dim>::assemble_stokes_system ()
{
  computing_timer.enter_section ("Stokes matrix assembly");
  pcout << "   Assembling..." << std::flush;

  if (rebuild_stokes_matrix == true)
    stokes_matrix=0;

  stokes_rhs=0;

  const QGauss<dim> quadrature_formula (parameters.stokes_degree+2);
  FEValues<dim>     stokes_fe_values (stokes_fe, quadrature_formula,
                                      update_values    |
                                      update_quadrature_points  |
                                      update_JxW_values |
                                      (rebuild_stokes_matrix == true
                                       ?
                                       update_gradients
                                       :
                                       UpdateFlags(0)));

  //FEValues<dim>     temperature_fe_values (temperature_fe, quadrature_formula,
                                           //update_values);

  const unsigned int   dofs_per_cell   = stokes_fe.dofs_per_cell;
  const unsigned int   n_q_points      = quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs    (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  // Next we need a vector that will contain the values of the temperature
  // solution at the previous time level at the quadrature points to
  // assemble the source term in the right hand side of the momentum
  // equation. Let's call this vector <code>old_solution_values</code>.
  //
  // The set of vectors we create next hold the evaluations of the basis
  // functions as well as their gradients and symmetrized gradients that
  // will be used for creating the matrices. Putting these into their own
  // arrays rather than asking the FEValues object for this information each
  // time it is needed is an optimization to accelerate the assembly
  // process, see step-22 for details.
  //
  // The last two declarations are used to extract the individual blocks
  // (velocity, pressure, temperature) from the total FE system.
  //std::vector<double>               old_temperature_values(n_q_points);

  std::vector<Tensor<1,dim> >          phi_u       (dofs_per_cell);
  std::vector<SymmetricTensor<2,dim> > grads_phi_u (dofs_per_cell);
  std::vector<double>                  div_phi_u   (dofs_per_cell);
  std::vector<double>                  phi_p       (dofs_per_cell);

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  // Now start the loop over all cells in the problem. We are working on two
  // different DoFHandlers for this assembly routine, so we must have two
  // different cell iterators for the two objects in use. This might seem a
  // bit peculiar, since both the Stokes system and the temperature system
  // use the same grid, but that's the only way to keep degrees of freedom
  // in sync. The first statements within the loop are again all very
  // familiar, doing the update of the finite element data as specified by
  // the update flags, zeroing out the local arrays and getting the values
  // of the old solution at the quadrature points. Then we are ready to loop
  // over the quadrature points on the cell.
  typename DoFHandler<dim>::active_cell_iterator
  cell = stokes_dof_handler.begin_active(),
  endc = stokes_dof_handler.end();

  //typename DoFHandler<dim>::active_cell_iterator
  //temperature_cell = temperature_dof_handler.begin_active();

  for (; cell!=endc; ++cell)
  if(cell->is_locally_owned())
    {
      stokes_fe_values.reinit (cell);
      //temperature_fe_values.reinit (temperature_cell);

      local_matrix = 0;
      local_rhs = 0;

      //temperature_fe_values.get_function_values (old_temperature_solution,
                                                 //old_temperature_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          //const double old_temperature = old_temperature_values[q];

          // Next we extract the values and gradients of basis functions
          // relevant to the terms in the inner products. As shown in
          // step-22 this helps accelerate assembly.
          //
          // Once this is done, we start the loop over the rows and columns
          // of the local matrix and feed the matrix with the relevant
          // products. The right hand side is filled with the forcing term
          // driven by temperature in direction of gravity (which is
          // vertical in our example).  Note that the right hand side term
          // is always generated, whereas the matrix contributions are only
          // updated when it is requested by the
          // <code>rebuild_matrices</code> flag.
          for (unsigned int k=0; k<dofs_per_cell; ++k)
            {
              phi_u[k] = stokes_fe_values[velocities].value (k,q);
              if (rebuild_stokes_matrix)
                {
                  grads_phi_u[k] = stokes_fe_values[velocities].symmetric_gradient(k,q);
                  div_phi_u[k]   = stokes_fe_values[velocities].divergence (k, q);
                  phi_p[k]       = stokes_fe_values[pressure].value (k, q);
                }
            }

          if (rebuild_stokes_matrix)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                local_matrix(i,j) += (2 * parameters.viscosity * 
                                      (grads_phi_u[i] * grads_phi_u[j])
                                      - div_phi_u[i] * phi_p[j]
                                      - phi_p[i] * div_phi_u[j])
                                     * stokes_fe_values.JxW(q);

          const Point<dim> gravity_vector = parameters.gravity * parameters.density *
                                 ( (dim == 2) ? (Point<dim> (0,1)) :
                                        (Point<dim> (0,0,1)) );
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            local_rhs(i) += phi_u[i] * gravity_vector * stokes_fe_values.JxW(q);
            //local_rhs(i) = 0;
        }

      // The last step in the loop over all cells is to enter the local
      // contributions into the global matrix and vector structures to the
      // positions specified in <code>local_dof_indices</code>.  Again, we
      // let the ConstraintMatrix class do the insertion of the cell matrix
      // elements to the global matrix, which already condenses the hanging
      // node constraints.
      cell->get_dof_indices (local_dof_indices);

      if (rebuild_stokes_matrix == true)
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

    if (rebuild_stokes_matrix == true)
      stokes_matrix.compress(VectorOperation::add);

    stokes_rhs.compress(VectorOperation::add);

    rebuild_stokes_matrix = false;

    pcout << std::endl;
    computing_timer.exit_section();

}//end_stokes_system

