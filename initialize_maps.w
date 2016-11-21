template<int dim>
void SD<dim>::initialize_maps()
{
  TimerOutput::Scope timer_section(computing_timer, "Initialize maps");

  stokes_flux_map.attach_worker_id(worker_id);
  stokes_flux_map.attach_source_dof_handler(stokes_dof_handler);
  stokes_flux_map.attach_target_dof_handler(flux_dof_handler);
  stokes_flux_map.set_boundary_indicator(interface_id);
  stokes_flux_map.build_maps();
  //stokes_flux_map.test_mapping();
  
  darcy_flux_map.attach_worker_id(worker_id);
  darcy_flux_map.attach_source_dof_handler(darcy_dof_handler);
  darcy_flux_map.attach_target_dof_handler(flux_dof_handler);
  darcy_flux_map.set_boundary_indicator(interface_id);
  darcy_flux_map.build_maps();
  //darcy_flux_map.test_mapping();

  stokes_lambda_map.attach_worker_id(worker_id);
  stokes_lambda_map.attach_source_dof_handler(stokes_dof_handler);
  stokes_lambda_map.attach_target_dof_handler(lambda_dof_handler);
  stokes_lambda_map.set_boundary_indicator(interface_id);
  stokes_lambda_map.build_maps();
  //stokes_lambda_map.test_mapping();

  darcy_lambda_map.attach_worker_id(worker_id);
  darcy_lambda_map.attach_source_dof_handler(darcy_dof_handler);
  darcy_lambda_map.attach_target_dof_handler(lambda_dof_handler);
  darcy_lambda_map.set_boundary_indicator(interface_id);
  darcy_lambda_map.build_maps();
  //darcy_lambda_map.test_mapping();


  //if(!is_RT)
  //{
    //stokes_flux_map.sort_source_interface_dofs_by_x_coordinate();
    //darcy_flux_map.sort_source_interface_dofs_by_x_coordinate();
  //}

}
