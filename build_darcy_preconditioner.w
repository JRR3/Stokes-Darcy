template <int dim>
void SD<dim>::build_darcy_preconditioner ()
{

  if (rebuild_darcy_preconditioner == false)
    return;

  TimerOutput::Scope timer_section(computing_timer, "Darcy build preconditioner");
  //Timer is called in the function: assemble_darcy_preconditioner();

  pcout << "   Rebuilding Darcy preconditioner..." << std::endl;

  //Darcy_original_bc and Darcy_zero_bc
  assemble_darcy_preconditioner ();

  std::vector<std::vector<bool> > constant_modes;
  FEValuesExtractors::Scalar pressure(dim);
  DoFTools::extract_constant_modes (darcy_dof_handler,
                                    darcy_fe.component_mask(pressure),
                                    constant_modes);

  darcy_Mu_preconditioner.reset  (new Precondition_A_Darcy());
  darcy_Amg_preconditioner.reset (new Precondition_S_Darcy());

  //Zero bc preconditioner
  //darcy_zero_bc_Mu_preconditioner.reset  (new TrilinosWrappers::PreconditionJacobi());
  //darcy_zero_bc_Amg_preconditioner.reset (new Precondition_S_Darcy());

  //CHECK: Do we need this in the case of scalar problems.
  Precondition_S_Darcy::AdditionalData Amg_data;
  Amg_data.constant_modes = constant_modes;
  Amg_data.elliptic = true;
  Amg_data.higher_order_elements = true;
  Amg_data.smoother_sweeps = 2;
  Amg_data.aggregation_threshold = 0.02;

  darcy_Mu_preconditioner->initialize (darcy_preconditioner_matrix.block(0,0));
  darcy_Amg_preconditioner->initialize (darcy_preconditioner_matrix.block(1,1), Amg_data);
  //darcy_Amg_preconditioner->initialize (darcy_preconditioner_matrix.block(1,1), Amg_data);

  //darcy_zero_bc_Mu_preconditioner->initialize (darcy_preconditioner_matrix.block(1,1));
  //darcy_zero_bc_Amg_preconditioner->initialize (darcy_preconditioner_matrix.block(0,0), Amg_data);

  rebuild_darcy_preconditioner         = false;
  rebuild_darcy_zero_bc_preconditioner = false;

}//end_build_darcy_preconditioner

