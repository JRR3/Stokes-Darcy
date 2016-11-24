template<int dim>
void SD<dim>::initialize_parallel_maps()
{
  TimerOutput::Scope timer_section(computing_timer, "Initialize parallel maps");

  parallel_map.attach_relevant_objects(
                            stokes_dof_handler,
                            darcy_dof_handler,
                            interface_id,
                            flux_refinements);

  //stokes_parallel_flux_map.create_local_mesh();

}
