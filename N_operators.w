template<int dim>
void SD<dim>::N_star(Vector<double> &lambda_target, const Vector<double> &rho_source)
{
  TimerOutput::Scope timer_section(computing_timer, "N star");
  //BlockVector<double> stokes_temp (stokes_vector);
  //BlockVector<double> darcy_temp  (darcy_vector);
  stokes_star(stokes_vector, rho_source);
  darcy_star (darcy_vector,  rho_source);

  compute_flux_function(lambda_target, stokes_vector, darcy_vector);

  //MODIFIED: Only for Q2xQ1 Darcy
  //DOCUMENT
  //extract_stokes_and_darcy_interface_values();
  //lambda_target = lambda_interface;


  //Add the contribution from the lambda piece present in the
  //rho_source vector.

  //MODIFIED->REQUIRED
  //lambda_target.sadd(1,sqrt_delta, rho_source.second);
}
//---------------------------------------------------------
template<int dim>
void SD<dim>::N_prime( Vector<double> &rho_target, const Vector<double> &lambda_source)
{
  TimerOutput::Scope timer_section(computing_timer, "N prime");
  rho_target = 0;
  //BlockVector<double> stokes_temp (stokes_vector);
  //BlockVector<double> darcy_temp  (darcy_vector);
  stokes_prime     (stokes_vector, lambda_source);
  darcy_prime      (darcy_vector , lambda_source);
  compute_rho_flux (rho_target, stokes_vector, darcy_vector);

  //MODIFIED: Only for Q2xQ1 Darcy
  //DOCUMENT
  //extract_stokes_and_darcy_interface_values();
  //rho_target = flux_vector;

  //Scaling// Is it necessary?
  rho_target *= i_sqrt_flux_cell_measure;
}
//---------------------------------------------------------
template<int dim>
void SD<dim>::N_operator ( Vector<double> &rho_target, const Vector<double> &lambda_source)
{
  TimerOutput::Scope timer_section(computing_timer, "N operator");
  rho_target = 0;

  //darcy_operator  (darcy_vector , lambda_source);
  darcy_operator  ();
  //stokes_operator (stokes_vector , lambda_source);
  stokes_operator ();
  pcout << "Safe exit from Darcy and Stokes operators" << std::endl;
  return;
  compute_rho_flux(rho_target, stokes_vector , darcy_vector );

  //MODIFIED: Only for Q2xQ1 Darcy
  //extract_stokes_and_darcy_interface_values();
  //rho_target = flux_vector;


  //Scaling// Is it necessary?
  rho_target  *= i_sqrt_flux_cell_measure;
}
