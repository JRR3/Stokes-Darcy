template<int dim>
void SD<dim>::initialize_parallel_maps()
{
  TimerOutput::Scope timer_section(computing_timer, "Initialize parallel maps");

  unsigned int ref = 0;
  lambda_parallel_map.reinit( stokes_dof_handler,
                       darcy_dof_handler,
                       lambda_fe,
                       interface_id,
                       ref);

  //stokes_parallel_flux_map.create_local_mesh();

}
