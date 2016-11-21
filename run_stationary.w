//-------------------------------------------------------------
template <int dim>
void SD<dim>::run_lambda_stokes_darcy ()
{


  cr_stokes.add_field("uh1");
  cr_stokes.add_field("ul2");
  cr_stokes.add_field("pl2");
  //--------------------------------
  cr_darcy.add_field("uhdiv");
  cr_darcy.add_field("ul2");
  cr_darcy.add_field("pl2");


  //const double pi = numbers::PI;
  std::cout << "+++++Run S+D stationary rates ++++++++++ " << std::endl;
  create_darcy_grid();
  create_stokes_grid();
  create_lambda_grid();
  create_flux_grid();
  refine_all(1);


  time       = 0;
  iteration  = 0;

  //cg_stats.attach_objects
    //(flux_dof_handler, lambda_dof_handler, flux_vector, lambda_vector, l2_norm_obj);

  for(unsigned int iter = 0; iter < 3; ++iter)
  {
    //for(unsigned int w = 0; w < cr_stokes.size(); ++w)
    //{
      //EV_stokes[w].clear();
      //EV_darcy[w].clear();
    //}

    //lambda_to_stokes.clear();
    //lambda_to_darcy.clear();
    //stokes_to_lambda.clear();
    //darcy_to_lambda.clear();

    stokes_flux_map.clear();
    darcy_flux_map.clear();
    stokes_lambda_map.clear();
    darcy_lambda_map.clear();
    //flux_to_stokes.clear();
    //flux_to_darcy.clear();
    //stokes_to_flux.clear();
    //darcy_to_flux.clear();

    refine_all(3);
    setup_stokes_dofs();
    setup_darcy_dofs();
    setup_lambda_dofs();
    setup_flux_dofs();

    //************************fe_tests
    //fe_tests();
    //return;

    //lambda_connection();
    initialize_maps();
    //connect_stokes_and_darcy_to_lambda();

    //test_flux_connection();
    //interpolate_initial_data();
    print_basic_stats();

    interpolate_interfacial_pressure();

    ++iteration;
    //Vector<double> lambda_approx ( lambda_vector.size() );
    //lambda_approx = lambda_vector;
    //std::cout << "--------------------Iteration = " << iteration << std::endl;
    //------------------------------------
    //update_extrapolation_vectors();
    //------------------------------------
    //create_stokes_constraints();
    //create_darcy_constraints();
    //------------------------------------
    //lambda_approx = -0.1;
    //compare_darcy_data();
    //return;
    //extract_stokes_and_darcy_interface_values();
    ready_matrices_and_preconditioners();
    cg_lsq();
    return;
    //inner_product_tests();
    //test_adjoint_profile();
    //test_prime_profile();
    //test_main_profile();
    //return;
    //------------------------------------

    //Computing approx SD solution
    stokes_operator();
    darcy_operator();



    ////------------------------------------
    //stokes_solution = stokes_vector;
    //darcy_solution  = darcy_vector;

    //------------------------------------Errors
    cr_stokes.compute_pressure_integral(stokes_dof_handler, q_rule, stokes_vector);
    cr_darcy.compute_pressure_integral(darcy_dof_handler, q_rule, darcy_vector);
    cr_stokes.compute_combined_mean_p(cr_darcy.p_domain, cr_darcy.measure);
    cr_darcy.compute_combined_mean_p(cr_stokes.p_domain, cr_stokes.measure);

    cr_stokes.compute_norms(stokes_dof_handler, q_rule, stokes_vector, stokes_exact_solution);
    cr_darcy.compute_norms(darcy_dof_handler, q_rule, darcy_vector, darcy_exact_solution);

    //return;



    //PlotData<dim,dim,MPI_BlockVector >::write_vtk(
                                         //darcy_dof_handler,
                                         //darcy_exact_solution,
                                         //darcy_solution,
                                         //"compare");
    //return;


    double h_param = std::max(GridTools::maximal_cell_diameter(stokes_domain),
                              GridTools::maximal_cell_diameter(darcy_domain));

    cr_stokes.save_h_value(h_param);
    cr_darcy.save_h_value(h_param);

    std::cout << "Finished one cycle" << std::endl;
    
    }//end_while

  //-----------------------------------------------------Summary
  std::cout << "Computing Stokes rates..." << std::endl;
  //cr_stokes.compute_stationary_rates(EV_stokes);
  cr_stokes.compute_stationary_rates();
  std::cout << "Computing Darcy  rates..." << std::endl;
  //cr_darcy.compute_stationary_rates(EV_darcy);
  cr_darcy.compute_stationary_rates();
  //-------------------------------------------------
  //cg_stats.print_solver_stats();
}

//-------------------------------------------------------------
template <int dim>
void SD<dim>::refine_all (const unsigned int &refinements = 1)
{
  stokes_domain.refine_global(refinements);
  darcy_domain. refine_global(refinements);
  lambda_domain.refine_global(refinements);
  flux_domain.  refine_global(refinements);

}
//-------------------------------------------------------------
//-------------------------------------------------------------
template <int dim>
void SD<dim>::update_old_vectors ()
{

}
//-------------------------------------------------------------
//-------------------------------------------------------------
template <int dim>
void SD<dim>::interpolate_initial_data ()
{

}
//-------------------------------------------------------------
template <int dim>
void SD<dim>::interpolate_interfacial_pressure ()
{
    ip_function.set_time(time);
    VectorTools::interpolate(lambda_dof_handler, ip_function, lambda_vector);  
    //VectorTools::interpolate(flux_dof_handler, ip_function, flux_solution);  
    //MINUS INCLUDED
    lambda_vector *= -1;
}
//-------------------------------------------------------------
