/** This function assembles the Darcy matrix and preconditioner.
 * In the Stokes setting we only assemble the preconditioner.
 * The Stokes matrix is assembled every time we call the stokes_operator()
 * function.
 */
template<int dim>
void SD<dim>::ready_matrices_and_preconditioners()
{
  rebuild_darcy_matrix = true;
  rebuild_darcy_preconditioner = true;

  assemble_darcy_system ();
  build_darcy_preconditioner ();

  build_stokes_preconditioner ();
  //stokes_operator() includes the assembly of the stokes system.
}
