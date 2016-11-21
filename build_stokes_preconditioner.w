template <int dim>
void SD<dim>::build_stokes_preconditioner ()
{

  if (rebuild_stokes_preconditioner == false)
    return;

  TimerOutput::Scope timer_section(computing_timer, "Stokes build preconditioner");
  //Timer is called in the function: assemble_stokes_preconditioner();

  pcout << "   Rebuilding Stokes preconditioner..." << std::endl;

  //Stokes_original_bc and Stokes_zero_bc
  assemble_stokes_preconditioner ();

  std::vector<std::vector<bool> > constant_modes;
  FEValuesExtractors::Vector velocity_components(0);
  DoFTools::extract_constant_modes (stokes_dof_handler,
                                    stokes_fe.component_mask(velocity_components),
                                    constant_modes);

  stokes_Mp_preconditioner.reset  (new Precondition_S_Stokes());
  stokes_Amg_preconditioner.reset (new Precondition_A_Stokes());

  //Zero bc preconditioner
  stokes_zero_bc_Mp_preconditioner.reset  (new Precondition_S_Stokes());
  stokes_zero_bc_Amg_preconditioner.reset (new Precondition_A_Stokes());

  Precondition_A_Stokes::AdditionalData Amg_data;
  Amg_data.constant_modes = constant_modes;
  Amg_data.elliptic = true;
  Amg_data.higher_order_elements = true;
  Amg_data.smoother_sweeps = 2;
  Amg_data.aggregation_threshold = 0.02;

  stokes_Mp_preconditioner->initialize (stokes_preconditioner_matrix.block(1,1));
  stokes_Amg_preconditioner->initialize (stokes_preconditioner_matrix.block(0,0), Amg_data);

  stokes_zero_bc_Mp_preconditioner->initialize (stokes_preconditioner_matrix.block(1,1));
  stokes_zero_bc_Amg_preconditioner->initialize (stokes_preconditioner_matrix.block(0,0), Amg_data);

  rebuild_stokes_preconditioner         = false;
  rebuild_stokes_zero_bc_preconditioner = false;

  pcout << std::endl;

}//end_build_stokes_preconditioner

