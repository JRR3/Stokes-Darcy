//-------------------------------------------------------------
template<int dim>
void SD<dim>::create_flux_grid()
{
    TimerOutput::Scope timer_section(computing_timer, "Flux grid");

    //GridGenerator::hyper_rectangle (flux_domain, Point<dim>(0,1), Point<dim>(1,1));
    GridGenerator::hyper_rectangle (flux_domain, Point<dim> (0,0), Point<dim>(1,0));
    GridTools::shift(Point<dim> (0,1), flux_domain);

    typename Triangulation<dim-1,dim>::active_cell_iterator cell = flux_domain.begin(),
      ecell = flux_domain.end();
    for(; cell != ecell; ++cell)
    if(cell->at_boundary())
    for(unsigned int k = 0; k < GeometryInfo<dim-1>::faces_per_cell; ++k)
    if(cell->face(k)->at_boundary())
    {
      //pcout << "Position: " << cell->face(k)->center() 
                //<< " is @ bdry." << std::endl;
      cell->face(k)->set_boundary_id(0);
    }

    //Extra level with respect to Darcy and Stokes domains.
    flux_domain.refine_global(flux_refinements);

}
//----------------------------------------------------------------
template <int dim>
void SD<dim>::setup_flux_dofs ()
{
  pcout << ">>>Setup Flux dofs" << std::endl;

  flux_dof_handler.distribute_dofs (flux_fe);

  unsigned int flux_n_dofs = flux_dof_handler.n_dofs();

  double flux_domain_volume   = GridTools::volume(flux_domain);
  double lambda_domain_volume = GridTools::volume(lambda_domain);
  flux_cell_measure           = flux_domain_volume / flux_domain.n_active_cells();
  lambda_cell_measure         = lambda_domain_volume / lambda_domain.n_active_cells();
  sqrt_flux_cell_measure      = sqrt(flux_cell_measure);
  i_sqrt_flux_cell_measure    = 1./sqrt_flux_cell_measure;


  flux_vector.reinit(flux_n_dofs );
  flux_solution.reinit(flux_dof_handler.n_dofs());


}
//----------------------------------------------------------------
