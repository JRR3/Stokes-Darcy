//-------------------------------------------------------------
/** The lambda grid should match the Darcy and Stokes grids along the interface.
 * It is possible to have non matching faces between the Darcy-Lambda-Stokes regions
 * at this point it has not been necessary to use that feature.
 */
//-------------------------------------------------------------
template<>
void SD<2>::create_lambda_grid()
{
    TimerOutput::Scope timer_section(computing_timer, "Lambda grid");

    //GridGenerator::hyper_rectangle (lambda_domain, Point<2>(0,1), Point<2>(1,1));
    GridGenerator::hyper_rectangle (lambda_domain, Point<2>(0,0), Point<2>(1,0));
    GridTools::shift(Point<2> (0,1), lambda_domain);
}
//-------------------------------------------------------------
template<>
void SD<3>::create_lambda_grid()
{
    TimerOutput::Scope timer_section(computing_timer, "Lambda grid");

    GridGenerator::hyper_rectangle (lambda_domain, Point<3>(0,0,0), Point<3>(1,1,0));
    GridTools::shift(Point<3> (0,0,1), lambda_domain);
}
//----------------------------------------------------------------
template <int dim>
void SD<dim>::setup_lambda_dofs ()
{
  pcout << ">>>Setup Lambda dofs" << std::endl;

  lambda_dof_handler.distribute_dofs (lambda_fe);

  lambda_locally_owned_dofs = lambda_dof_handler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (lambda_dof_handler, lambda_locally_relevant_dofs);

  lambda_vector.reinit (lambda_locally_owned_dofs, MPI_COMM_WORLD);


}
//----------------------------------------------------------------
