//-------------------------------------------------------------
/** The lambda grid should match the Darcy and Stokes grids along the interface.
 * It is possible to have non matching faces between the Darcy-Lambda-Stokes regions
 * at this point it has not been necessary to use that feature.
 */
//-------------------------------------------------------------
template<int dim>
void SD<dim>::create_lambda_grid()
{
    TimerOutput::Scope timer_section(computing_timer, "Lambda grid");

    //GridGenerator::hyper_rectangle (lambda_domain, Point<dim>(0,1), Point<dim>(1,1));
    GridGenerator::hyper_rectangle (lambda_domain, Point<dim>(0,0), Point<dim>(1,0));
    GridTools::shift(Point<dim> (0,1), lambda_domain);

    typename Triangulation<dim-1,dim>::active_cell_iterator cell = lambda_domain.begin(),
      ecell = lambda_domain.end();
    for(; cell != ecell; ++cell)
    if(cell->at_boundary())
    for(unsigned int k = 0; k < GeometryInfo<dim-1>::faces_per_cell; ++k)
    if(cell->face(k)->at_boundary())
    {
      //pcout << "Position: " << cell->face(k)->center() 
                //<< " is @ bdry." << std::endl;
      cell->face(k)->set_boundary_id(0);
    }

}
//----------------------------------------------------------------
template <int dim>
void SD<dim>::setup_lambda_dofs ()
{
  pcout << ">>>Setup Lambda dofs" << std::endl;

  lambda_dof_handler.distribute_dofs (lambda_fe);

  unsigned int n_lambda_dofs = lambda_dof_handler.n_dofs();

  lambda_vector.reinit (n_lambda_dofs );
  lambda_zero.reinit   (n_lambda_dofs );
  lambda_zero = 0;


}
//----------------------------------------------------------------
