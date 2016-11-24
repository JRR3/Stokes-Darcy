//-------------------------------------------------------------
//-------------------------------------------------------------
template<>
void SD<2>::create_darcy_grid()
{
    TimerOutput::Scope timer_section(computing_timer, "Darcy grid");

    GridGenerator::hyper_rectangle (darcy_domain, Point<2>(0,0), Point<2>(1,1));
    typename Triangulation<2>::active_cell_iterator cell = darcy_domain.begin_active();
    if(!darcy_bc_are_enforced_weakly)
    {
      cell->face(0)->set_boundary_id(1);
      cell->face(1)->set_boundary_id(2);
      cell->face(2)->set_boundary_id(3);
      cell->face(3)->set_boundary_id(interface_id);
    }
    else
    {
      cell->face(0)->set_boundary_id(wall_id);
      cell->face(1)->set_boundary_id(wall_id);
      cell->face(2)->set_boundary_id(wall_id);
      cell->face(3)->set_boundary_id(interface_id);
    }
}
//-------------------------------------------------------------
//-------------------------------------------------------------
template<>
void SD<3>::create_darcy_grid()
{
    TimerOutput::Scope timer_section(computing_timer, "Darcy grid");

    GridGenerator::hyper_rectangle (darcy_domain, Point<3>(0,0,0), Point<3>(1,1,1));
    typename Triangulation<3>::active_cell_iterator cell = darcy_domain.begin_active();

    //Assume Darcy boundary conditions are enforced weakly
    cell->face(0)->set_boundary_id(wall_id);
    cell->face(1)->set_boundary_id(wall_id);
    cell->face(2)->set_boundary_id(wall_id);
    cell->face(3)->set_boundary_id(wall_id);
    cell->face(4)->set_boundary_id(wall_id);
    cell->face(5)->set_boundary_id(interface_id);
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
template <int dim>
void SD<dim>::
setup_darcy_matrix ( const std::vector<IndexSet> &darcy_partitioning,
                     const std::vector<IndexSet> &darcy_relevant_partitioning)
{
  //TimerOutput::Scope timer_section(computing_timer, "Stokes setup matrix");
  //Called inside Stokes setup dofs

  darcy_matrix.clear ();

  TrilinosWrappers::BlockSparsityPattern sp(darcy_partitioning, darcy_partitioning,
                                            darcy_relevant_partitioning,
                                            MPI_COMM_WORLD);

  Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
  for (unsigned int c=0; c<dim+1; ++c)
    for (unsigned int d=0; d<dim+1; ++d)
      if (! ((c==dim) && (d==dim)))
        coupling[c][d] = DoFTools::always;
      else
        coupling[c][d] = DoFTools::none;

  DoFTools::make_sparsity_pattern (darcy_dof_handler,
                                   coupling, sp,
                                   darcy_constraints, false,
                                   Utilities::MPI::
                                   this_mpi_process(MPI_COMM_WORLD));
  sp.compress();

  darcy_matrix.reinit (sp);
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
template <int dim>
void SD<dim>::setup_porosity_matrix ()
{
  porosity_locally_owned_dofs = porosity_dof_handler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (porosity_dof_handler, porosity_locally_relevant_dofs);

  //Porosity boundary conditions
  {
    porosity_constraints.clear ();
    porosity_constraints.reinit (porosity_locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (porosity_dof_handler, porosity_constraints);
    //VectorTools::interpolate_boundary_values (porosity_dof_handler,
                                              //0,
                                              //ZeroFunction<dim>(),
                                              //porosity_constraints);
    porosity_constraints.close ();
  }

  porosity_matrix.clear();
  DynamicSparsityPattern dsp (porosity_locally_relevant_dofs);

  DoFTools::make_sparsity_pattern (porosity_dof_handler, dsp,
                                   porosity_constraints, false);
  SparsityTools::distribute_sparsity_pattern (dsp, 
     porosity_dof_handler.n_locally_owned_dofs_per_processor(),
                                              MPI_COMM_WORLD,
                                              porosity_locally_relevant_dofs);

  porosity_matrix.reinit (porosity_locally_owned_dofs,
                          porosity_locally_owned_dofs,
                          dsp,
                          MPI_COMM_WORLD);

  porosity_rhs.reinit (porosity_locally_owned_dofs, MPI_COMM_WORLD);
  //Porosity rhs
  //
}
//-----------------------------------------------------------------------------
template <int dim>
void SD<dim>::setup_darcy_dofs ()
{
  TimerOutput::Scope timer_section(computing_timer, "Darcy setup dofs");

  std::vector<unsigned int> darcy_block_component (dim+1,0);
  darcy_block_component[dim] = 1;
  {
    darcy_dof_handler.distribute_dofs (darcy_fe);
    //DoFRenumbering::Cuthill_McKee (darcy_dof_handler);
    DoFRenumbering::component_wise (darcy_dof_handler, darcy_block_component);

  }

  //Porosity constraints
  {
    porosity_dof_handler.distribute_dofs (porosity_fe);
    DoFTools::extract_locally_relevant_dofs (porosity_dof_handler, porosity_locally_relevant_dofs);
    porosity_constraints.clear ();
    porosity_constraints.reinit (porosity_locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (porosity_dof_handler, porosity_constraints);
    porosity_constraints.close ();
  }

  std::vector<types::global_dof_index> darcy_dofs_per_block (2);
  DoFTools::count_dofs_per_block (darcy_dof_handler, darcy_dofs_per_block, darcy_block_component);
  const unsigned int n_u  = darcy_dofs_per_block[0],
                     n_p  = darcy_dofs_per_block[1];
                     //n_po = porosity_dof_handler.n_dofs();

  //Darcy_locally_relevant_dofs
  std::vector<IndexSet> darcy_partitioning, darcy_relevant_partitioning;
  IndexSet darcy_relevant_set;
  {
    IndexSet darcy_index_set = darcy_dof_handler.locally_owned_dofs();
    darcy_partitioning.push_back(darcy_index_set.get_view(0,n_u));
    darcy_partitioning.push_back(darcy_index_set.get_view(n_u,n_u+n_p));

    DoFTools::extract_locally_relevant_dofs (darcy_dof_handler,
                                             darcy_relevant_set);
    darcy_relevant_partitioning.push_back(darcy_relevant_set.get_view(0,n_u));
    darcy_relevant_partitioning.push_back(darcy_relevant_set.get_view(n_u,n_u+n_p));
  }


  //Darcy_constraints
  {
    darcy_constraints.clear ();
    darcy_constraints.reinit (darcy_relevant_set);
    DoFTools::make_hanging_node_constraints (darcy_dof_handler, darcy_constraints);
    darcy_constraints.close ();
  }

  //Darcy preconditioner constraints
  //Darcy_laplace_matrix_constraints
  {
    darcy_preconditioner_constraints.clear ();
    darcy_preconditioner_constraints.reinit(darcy_relevant_set);

    FEValuesExtractors::Scalar pressure(dim);

    DoFTools::make_hanging_node_constraints  (darcy_dof_handler, darcy_preconditioner_constraints);
    DoFTools::make_zero_boundary_constraints (darcy_dof_handler, darcy_preconditioner_constraints,
                                              darcy_fe.component_mask(pressure));

    darcy_preconditioner_constraints.close ();
  }

  setup_darcy_matrix         (darcy_partitioning, darcy_relevant_partitioning);
  setup_darcy_preconditioner (darcy_partitioning, darcy_relevant_partitioning);
  setup_porosity_matrix      ();


  darcy_rhs.reinit      (darcy_partitioning, MPI_COMM_WORLD);
  darcy_solution.reinit (darcy_partitioning, MPI_COMM_WORLD);

  //CHECK: What happens if we remove the fisrt argument?
  darcy_vector.reinit (darcy_partitioning, darcy_relevant_partitioning, MPI_COMM_WORLD);

  rebuild_darcy_matrix                      = true;
  rebuild_darcy_preconditioner              = true;
  rebuild_darcy_zero_bc_matrix              = true;
  rebuild_darcy_zero_bc_preconditioner      = true;
}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
template <int dim>
void SD<dim>::
setup_darcy_preconditioner ( const std::vector<IndexSet> &darcy_partitioning,
                             const std::vector<IndexSet> &darcy_relevant_partitioning)
{
  darcy_Amg_preconditioner.reset ();
  darcy_Mu_preconditioner.reset ();

  darcy_preconditioner_matrix.clear ();

  TrilinosWrappers::BlockSparsityPattern sp(darcy_partitioning, darcy_partitioning,
                                            darcy_relevant_partitioning,
                                            MPI_COMM_WORLD);

  Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
  for (unsigned int c=0; c<dim+1; ++c)
  for (unsigned int d=0; d<dim+1; ++d)
      if (c == d)
        coupling[c][d] = DoFTools::always;
      else
        coupling[c][d] = DoFTools::none;

  DoFTools::make_sparsity_pattern (darcy_dof_handler,
                                   coupling, sp,
                                   darcy_constraints, false,
                                   Utilities::MPI::
                                   this_mpi_process(MPI_COMM_WORLD));
  sp.compress();

  darcy_preconditioner_matrix.reinit (sp);
}
//---------------------------------------------------------------------

