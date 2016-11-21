/** This function summarizes the basic statistics of the program.
 * Number of cells
 * Degrees of freedom 
 * h parameter
 * Refine levels
 * Number of cells on the boundary
 */
template<int dim>
void SD<dim>::print_basic_stats()
{
  //Stokes data
  pcout << "*******Stokes information********" << std::endl;
  pcout << "Number of active cells: "
            << stokes_domain.n_global_active_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << stokes_dof_handler.n_dofs()
            << " (" << n_u_s << '+' << n_p_s << ')'
            << std::endl;
  pcout << "                   h           : " 
            << GridTools::maximal_cell_diameter(stokes_domain) << std::endl;

  pcout << "*******Darcy information********" << std::endl;
  pcout << "Number of active cells: "
            << darcy_domain.n_global_active_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << darcy_dof_handler.n_dofs()
            << " (" << n_u_d << '+' << n_p_d << ')'
            << std::endl;
  pcout << "                   h           : " 
            << GridTools::maximal_cell_diameter(darcy_domain) << std::endl;

  std::cout << "*******Porosity information********" << std::endl;
  std::cout << "Number of active cells: "
            << darcy_domain.n_global_active_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << porosity_dof_handler.n_dofs()
            << std::endl;

  pcout << "*******Lambda information********" << std::endl;
  pcout << "Number of active cells: "
            << lambda_domain.n_active_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << lambda_dof_handler.n_dofs()
            << std::endl;
  pcout << "                   h           : " 
            << GridTools::maximal_cell_diameter(lambda_domain) << std::endl;

  pcout << "*******Flux information********" << std::endl;
  pcout << "Number of active cells: "
            << flux_domain.n_active_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << flux_dof_handler.n_dofs()
            << std::endl;
  pcout << "   flux_cell_measure           : " 
            << flux_cell_measure  << std::endl;
  pcout << "   i_sqrt_cell_measure         : " 
            << i_sqrt_flux_cell_measure  << std::endl;


  pcout << "*******Grand total********" << std::endl;
  unsigned int total_dofs = stokes_dof_handler.n_dofs() 
                          + darcy_dof_handler.n_dofs()
                          + lambda_dof_handler.n_dofs() 
                          + flux_dof_handler.n_dofs()
                          + porosity_dof_handler.n_dofs();
  pcout << total_dofs << " dofs in the system." << std::endl;

}
