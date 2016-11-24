//-------------------------------------------------------------
template<>
void SD<2>::create_flux_grid()
{
    TimerOutput::Scope timer_section(computing_timer, "Flux grid");

    //GridGenerator::hyper_rectangle (flux_domain, Point<2>(0,1), Point<2>(1,1));
    GridGenerator::hyper_rectangle (flux_domain, Point<2> (0,0), Point<2>(1,0));
    GridTools::shift(Point<2> (0,1), flux_domain);

    //Extra level with respect to Darcy and Stokes domains.
    flux_domain.refine_global(flux_refinements);

}
//-------------------------------------------------------------
template<>
void SD<3>::create_flux_grid()
{
    TimerOutput::Scope timer_section(computing_timer, "Flux grid");

    GridGenerator::hyper_rectangle (flux_domain, Point<3>(0,0,0), Point<3>(1,1,0));
    GridTools::shift(Point<3> (0,0,1), flux_domain);
    flux_domain.refine_global(flux_refinements);

}
//----------------------------------------------------------------
template <int dim>
void SD<dim>::setup_flux_dofs ()
{
  pcout << ">>>Setup Flux dofs" << std::endl;


  flux_dof_handler.distribute_dofs (flux_fe);

  flux_locally_owned_dofs = flux_dof_handler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (flux_dof_handler, flux_locally_relevant_dofs);

  //unsigned int flux_n_dofs = flux_dof_handler.n_dofs();
  //
  MeshData<dim-1,dim> mesh_data;
  mesh_data.compute_data(flux_dof_handler);

  //double lambda_domain_volume = GridTools::volume(lambda_domain);
  //pcout << "Lambda volume: " << lambda_domain_volume << std::endl;
  //pcout << "Flux   volume: " << flux_domain_volume << std::endl;
  flux_cell_measure           = mesh_data.volume / flux_domain.n_global_active_cells();
  //lambda_cell_measure         = lambda_domain_volume / lambda_domain.n_global_active_cells();
  sqrt_flux_cell_measure      = sqrt(flux_cell_measure);
  i_sqrt_flux_cell_measure    = 1./sqrt_flux_cell_measure;


  flux_vector.reinit (flux_locally_owned_dofs, MPI_COMM_WORLD);


}
//----------------------------------------------------------------
