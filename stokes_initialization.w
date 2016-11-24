//-------------------------------------------------------------
template<>
void SD<2>::create_stokes_grid()
{
    TimerOutput::Scope timer_section(computing_timer, "Stokes grid");

    GridGenerator::hyper_rectangle (stokes_domain, Point<2>(0,1), Point<2>(1,2));
    typename Triangulation<2>::active_cell_iterator 
      cell = stokes_domain.begin_active();
    cell->face(0)->set_boundary_id(wall_id);
    cell->face(1)->set_boundary_id(wall_id);
    cell->face(2)->set_boundary_id(interface_id);
    cell->face(3)->set_boundary_id(wall_id);
}
//-------------------------------------------------------------
//-------------------------------------------------------------
template<>
void SD<3>::create_stokes_grid()
{
    TimerOutput::Scope timer_section(computing_timer, "Stokes grid");

    GridGenerator::hyper_rectangle (stokes_domain, Point<3>(0,0,1), Point<3>(1,1,2));
    typename Triangulation<3>::active_cell_iterator cell = stokes_domain.begin_active();

    cell->face(0)->set_boundary_id(wall_id);
    cell->face(1)->set_boundary_id(wall_id);
    cell->face(2)->set_boundary_id(wall_id);
    cell->face(3)->set_boundary_id(wall_id);
    cell->face(4)->set_boundary_id(interface_id);
    cell->face(5)->set_boundary_id(wall_id);
}
//-------------------------------------------------------------
template<int dim>
void SD<dim>::create_stokes_zero_bc_constraints(IndexSet &stokes_relevant_set)
{
  stokes_zero_bc_constraints.clear ();

  stokes_zero_bc_constraints.reinit (stokes_relevant_set);

  FEValuesExtractors::Vector velocities(0);
  DoFTools::make_hanging_node_constraints (stokes_dof_handler,
                                           stokes_zero_bc_constraints);
  //Velocity
  VectorTools::interpolate_boundary_values (stokes_dof_handler,
                                            wall_id,
                                            ZeroFunction<dim> (dim+1),
                                            stokes_zero_bc_constraints,
                                            stokes_fe.component_mask(velocities));
  stokes_zero_bc_constraints.close ();
}
//-------------------------------------------------------
template<int dim>
void SD<dim>::create_stokes_constraints(IndexSet &stokes_relevant_set)
{
  stokes_constraints.clear ();

  stokes_constraints.reinit (stokes_relevant_set);

  FEValuesExtractors::Vector velocities(0);
  DoFTools::make_hanging_node_constraints (stokes_dof_handler,
                                           stokes_constraints);
  stokes_exact_solution.set_time(time);

  //Velocity
  VectorTools::interpolate_boundary_values (stokes_dof_handler,
                                            wall_id,
                                            stokes_exact_solution,
                                            stokes_constraints,
                                            stokes_fe.component_mask(velocities));
  stokes_constraints.close ();
}
//-------------------------------------------------------------
template <int dim>
void SD<dim>::setup_stokes_dofs ()
{
  TimerOutput::Scope timer_section(computing_timer, "Stokes setup dofs");

  std::vector<unsigned int> stokes_sub_blocks (dim+1,0);
  stokes_sub_blocks[dim] = 1;
  stokes_dof_handler.distribute_dofs (stokes_fe);
  DoFRenumbering::component_wise (stokes_dof_handler, stokes_sub_blocks);

  //stress_dof_handler.distribute_dofs (stress_fe);

  std::vector<types::global_dof_index> stokes_dofs_per_block (2);
  DoFTools::count_dofs_per_block (stokes_dof_handler, stokes_dofs_per_block,
                                  stokes_sub_blocks);

  const unsigned int n_u = stokes_dofs_per_block[0],
                     n_p = stokes_dofs_per_block[1];
                     //n_s = stress_dof_handler.n_dofs();

  std::vector<IndexSet> stokes_partitioning, stokes_relevant_partitioning;
  //IndexSet stress_partitioning (n_s), stress_relevant_partitioning (n_s);
  IndexSet stokes_relevant_set;
  {
    IndexSet stokes_index_set = stokes_dof_handler.locally_owned_dofs();
    stokes_partitioning.push_back(stokes_index_set.get_view(0,n_u));
    stokes_partitioning.push_back(stokes_index_set.get_view(n_u,n_u+n_p));

    DoFTools::extract_locally_relevant_dofs (stokes_dof_handler,
                                             stokes_relevant_set);
    stokes_relevant_partitioning.push_back(stokes_relevant_set.get_view(0,n_u));
    stokes_relevant_partitioning.push_back(stokes_relevant_set.get_view(n_u,n_u+n_p));

    //stress_partitioning = stress_dof_handler.locally_owned_dofs();
    //DoFTools::extract_locally_relevant_dofs (stress_dof_handler,
                                             //stress_relevant_partitioning);
  }

  //Create boundary condition matrices
  create_stokes_constraints        ( stokes_relevant_set );
  create_stokes_zero_bc_constraints( stokes_relevant_set );

  setup_stokes_matrix         (stokes_partitioning, stokes_relevant_partitioning);
  setup_stokes_preconditioner (stokes_partitioning, stokes_relevant_partitioning);
  //setup_temperature_matrices (temperature_partitioning,
                              //temperature_relevant_partitioning);
  setup_stokes_zero_bc_matrix         (stokes_partitioning, stokes_relevant_partitioning);
  setup_stokes_zero_bc_preconditioner (stokes_partitioning, stokes_relevant_partitioning);

  //stokes_rhs.reinit (stokes_partitioning, stokes_relevant_partitioning,
                     //MPI_COMM_WORLD, true);
  stokes_rhs.reinit       (stokes_partitioning, MPI_COMM_WORLD);
  stokes_solution.reinit  (stokes_partitioning, MPI_COMM_WORLD);
  stokes_vector.reinit    (stokes_partitioning, stokes_relevant_partitioning, MPI_COMM_WORLD);

  rebuild_stokes_matrix                      = true;
  rebuild_stokes_preconditioner              = true;
  rebuild_stokes_zero_bc_matrix              = true;
  rebuild_stokes_zero_bc_preconditioner      = true;
}

//-----------------------------------------------------------------------------
template <int dim>
void SD<dim>::
setup_stokes_matrix (const std::vector<IndexSet> &stokes_partitioning,
                     const std::vector<IndexSet> &stokes_relevant_partitioning)
{
  //TimerOutput::Scope timer_section(computing_timer, "Stokes setup matrix");
  //Called inside Stokes setup dofs

  stokes_matrix.clear ();

  TrilinosWrappers::BlockSparsityPattern sp(stokes_partitioning, stokes_partitioning,
                                            stokes_relevant_partitioning,
                                            MPI_COMM_WORLD);

  Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
  for (unsigned int c=0; c<dim+1; ++c)
    for (unsigned int d=0; d<dim+1; ++d)
      if (! ((c==dim) && (d==dim)))
        coupling[c][d] = DoFTools::always;
      else
        coupling[c][d] = DoFTools::none;

  DoFTools::make_sparsity_pattern (stokes_dof_handler,
                                   coupling, sp,
                                   stokes_constraints, false,
                                   Utilities::MPI::
                                   this_mpi_process(MPI_COMM_WORLD));
  sp.compress();

  stokes_matrix.reinit (sp);
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
template <int dim>
void SD<dim>::
setup_stokes_zero_bc_matrix (const std::vector<IndexSet> &stokes_partitioning,
                             const std::vector<IndexSet> &stokes_relevant_partitioning)
{
  //TimerOutput::Scope timer_section(computing_timer, "Stokes setup matrix");
  //Called inside Stokes setup dofs

  stokes_zero_bc_matrix.clear ();

  TrilinosWrappers::BlockSparsityPattern sp(stokes_partitioning, stokes_partitioning,
                                            stokes_relevant_partitioning,
                                            MPI_COMM_WORLD);

  Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
  for (unsigned int c=0; c<dim+1; ++c)
    for (unsigned int d=0; d<dim+1; ++d)
      if (! ((c==dim) && (d==dim)))
        coupling[c][d] = DoFTools::always;
      else
        coupling[c][d] = DoFTools::none;

  DoFTools::make_sparsity_pattern (stokes_dof_handler,
                                   coupling, sp,
                                   stokes_zero_bc_constraints, false,
                                   Utilities::MPI::
                                   this_mpi_process(MPI_COMM_WORLD));
  sp.compress();

  stokes_zero_bc_matrix.reinit (sp);
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------
template <int dim>
void SD<dim>::
setup_stokes_preconditioner (const std::vector<IndexSet> &stokes_partitioning,
                             const std::vector<IndexSet> &stokes_relevant_partitioning)
{
  stokes_Amg_preconditioner.reset ();
  stokes_Mp_preconditioner.reset ();

  stokes_preconditioner_matrix.clear ();

  TrilinosWrappers::BlockSparsityPattern sp(stokes_partitioning, stokes_partitioning,
                                            stokes_relevant_partitioning,
                                            MPI_COMM_WORLD);

  Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
  for (unsigned int c=0; c<dim+1; ++c)
  for (unsigned int d=0; d<dim+1; ++d)
      if (c == d)
        coupling[c][d] = DoFTools::always;
      else
        coupling[c][d] = DoFTools::none;

  DoFTools::make_sparsity_pattern (stokes_dof_handler,
                                   coupling, sp,
                                   stokes_constraints, false,
                                   Utilities::MPI::
                                   this_mpi_process(MPI_COMM_WORLD));
  sp.compress();

  stokes_preconditioner_matrix.reinit (sp);
}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
template <int dim>
void SD<dim>::
setup_stokes_zero_bc_preconditioner (const std::vector<IndexSet> &stokes_partitioning,
                                     const std::vector<IndexSet> &stokes_relevant_partitioning)
{
  stokes_zero_bc_Amg_preconditioner.reset ();
  stokes_zero_bc_Mp_preconditioner.reset ();

  stokes_zero_bc_preconditioner_matrix.clear ();

  TrilinosWrappers::BlockSparsityPattern sp(stokes_partitioning, stokes_partitioning,
                                            stokes_relevant_partitioning,
                                            MPI_COMM_WORLD);

  Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
  for (unsigned int c=0; c<dim+1; ++c)
  for (unsigned int d=0; d<dim+1; ++d)
      if (c == d)
        coupling[c][d] = DoFTools::always;
      else
        coupling[c][d] = DoFTools::none;

  DoFTools::make_sparsity_pattern (stokes_dof_handler,
                                   coupling, sp,
                                   stokes_zero_bc_constraints, false,
                                   Utilities::MPI::
                                   this_mpi_process(MPI_COMM_WORLD));
  sp.compress();

  stokes_zero_bc_preconditioner_matrix.reinit (sp);
}//end_stokes_zero_bc_preconditioner
