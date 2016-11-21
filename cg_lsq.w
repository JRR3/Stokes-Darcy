/** This function implements the CG least-squares method.
 * This function is used to solve for the interfacial pressure.
 * See the book:
 * Matrix computations by Charles F. Van Loan 
 */
template<int dim>
void SD<dim>::cg_lsq()
//void SD<dim>::cg_lsq( Vector<double> &lambda_source)
//void SD<dim>::cg_lsq(Vector<double> &lambda_target)
{
  pcout << "****** CG LSQ ******" << std::endl;
  //Initialize 
  //----------------------------------------------------(0)
  PlotIData<dim-1,dim,Vector<double> > plotter;
  plotter.attach_dof_handlers(flux_dof_handler, lambda_dof_handler);
  double sigma, tau;
  double norm_As_r, norm_As_r_sq, norm_A_p_sq, norm_As_r_new_sq, norm_As_r_new;
  lambda_vector = 0;
  flux_vector   = 0;
  Vector<double> lambda_init (lambda_vector);
  //Vector<double> true_lambda (lambda_source);
  //Vector<double> true_flux (flux_vector);
  Vector<double> p_n ( lambda_vector );
  Vector<double> x_n ( lambda_vector );
  Vector<double> As_r( lambda_vector );
  Vector<double> b   (flux_vector);
  Vector<double> r_n (flux_vector);
  Vector<double> A_p (flux_vector);
  unsigned int count    = 0;
  unsigned int max_iter = 500;
  //---
  //PlotData<dim,dim,MPI_BlockVector > darcy_plotter;
  //darcy_plotter.attach_std_dof_handler(darcy_dof_handler);
  //---
  //PlotData<dim,dim,MPI_BlockVector > stokes_plotter;
  //stokes_plotter.attach_std_dof_handler(stokes_dof_handler);
  //---
  //PlotData<dim-1,dim,Vector<double> > flux_plotter;
  //flux_plotter.attach_dof_handler(flux_dof_handler);
  //---
  //PlotData<dim-1,dim,Vector<double> > lambda_plotter;
  //lambda_plotter.attach_dof_handler(lambda_dof_handler);
  //*****************
  //CGStats<dim-1,dim,Vector<double> > cg_stats;
  //cg_stats.attach_objects(flux_dof_handler, lambda_dof_handler, flux_vector, lambda_vector, l2_norm_obj);
  //*****************

  //lambda_vector = lambda_source;

  //---------------------------------------------------Computing true_h solution
  //pcout << "Computing true solutions on current meshes" << std::endl;
  //darcy_operator  (darcy_vector, lambda_vector);
  //stokes_operator (stokes_vector, lambda_vector);
  //compute_rho_flux(flux_vector, stokes_vector, darcy_vector);
  //plotter.write_vtk(flux_vector, lambda_vector);
  //---
  //cg_stats.get_snapshot(compute_rho_flux(flux_vector, stokes_vector, darcy_vector));

  //----------------------------------------------------(1)
  pcout << "Status before excecuting the CG algorithm" << std::endl;
  lambda_init = -0.1;
  //darcy_operator  (darcy_vector , lambda_vector);
  //stokes_operator (stokes_vector, lambda_vector);
  ////---
  //cg_stats.get_snapshot(compute_rho_flux(flux_vector, stokes_vector, darcy_vector));
  ////---

  //pcout << "Flu2: " << (compute_rho_flux(flux_vector, stokes_vector, darcy_vector))
            //<< std::endl;
  //pcout << "Flux: " << sqrt(compute_rho_flux(flux_vector, stokes_vector, darcy_vector))
            //<< std::endl;
  //pcout << "Flux: " << flux_vector.norm_sqr()*flux_vector.size()
            //<< std::endl;
  //return;


  //lambda_init   = lambda_vector;
  lambda_vector = lambda_init;
  //This defines b = -(N(g_0) - [0, sqrt_delta*g_0])
  N_operator (b, lambda_vector);
  return;

  //pcout << "b flux" << std::endl;
  //for(auto &it : b)
    //pcout << it << std::endl;

  //flux_plotter.write_flux_txt(b, "cpp_vector");
  b       *= -1;
  //---------------------------------------------------Plotting_true_solution
  //stokes_plotter.write_vtk(stokes_vector, "stokes_op_");
  //darcy_plotter.write_vtk(darcy_vector, "darcy_op_");

  //sqrt(compute_product_norm_sq(b));
  
  //----------------------------------------------------(2)
  //Define initial correction direction h = x_0 = 0.01
  x_n = 0.01;
  //pcout << "x_n computation " << std::endl;
  //compute_lambda_norm(x_n);

  
  //----------------------------------------------------(3)
  //Compute initial residual
  //r_0 = b - Ax_0
  N_prime( r_n, x_n);
  //stokes_flux_map.compare_data("stokes_N_prime_rhs", stokes_rhs);
  //darcy_flux_map.compare_data("darcy_N_prime_rhs", darcy_rhs);
  //stokes_flux_map.compare_data("stokes_N_prime_sol", stokes_vector);
  //darcy_flux_map.compare_data("darcy_N_prime_sol", darcy_vector);
  //flux_plotter.write_flux_txt(r_n, "cpp_vector");
  //Slight difference of O(10^{-13})
  r_n.sadd(-1,b);

  //A.vmult(r_n,x_n);
  //r_n.sadd(-1,b);
  //flux_plotter.write_flux_txt(r_n, "cpp_vector");

    
  //----------------------------------------------------(4)
  //Compute initial p
  //p_0 = A*r_0
  //pcout << "p_n" << std::endl;
  N_star(p_n, r_n);

  //stokes_flux_map.compare_data("stokes_adj_sol", stokes_vector);
  //darcy_flux_map.compare_data("darcy_adj_sol", darcy_vector);
  //darcy_flux_map.compare_data("darcy_adj_rhs", darcy_rhs);
  //stokes_flux_map.compare_data("stokes_adj_rhs", stokes_rhs);
  //return;

  //stokes_flux_map.print_sorted_data(stokes_vector);
  //pcout << "Sorted data for Darcy vector" << std::endl;
  //darcy_flux_map.print_sorted_data(darcy_vector);
  
  //lambda_plotter.write_txt(p_n, "cpp_vector");
  //Similar profile, some overshoot on the right for C++
  //pcout << "p_0 computation " << std::endl;
  //As.vmult(p_n,r_n);
  //return;



  //----------------------------------------------------(5)
  //Compute initial norms
  //norm_As_r    = compute_lambda_norm(p_n, true);
  norm_As_r_sq = l2_norm_obj.compute_norm_sq(lambda_dof_handler, p_n);
  norm_As_r    = sqrt(norm_As_r_sq);


  //---------------Start CG-algorithm------------------
  while (count++ < max_iter)
  {
    plotter.save_res_norm(norm_As_r);
    //----------------
    if (norm_As_r < cg_eps)
    {
      //Not a real iteration of the algorithm.
      count--;
      pcout << "We are done" << std::endl;
      pcout << "------------------------Iteration: " << count << std::endl;
      pcout << "--------------------Residual norm: " << norm_As_r << std::endl;
      break;
    }
    else
    {
      pcout << "------------------------Iteration: " << count << std::endl;
      pcout << "--------------------Residual norm: " << norm_As_r << std::endl;
    }

    //---------------------------------------------------Plotting
    lambda_vector  = lambda_init;
    lambda_vector += x_n;
    darcy_operator  ();
    stokes_operator ();
    //---
    //cg_stats.get_snapshot(compute_rho_flux(flux_vector, stokes_vector, darcy_vector));
    plotter.save_flux_norm( compute_rho_flux(flux_vector, stokes_vector, darcy_vector) );
    plotter.write_vtk(flux_vector, lambda_vector);
    

    //CAREFUL WITH NORMS (ITS A VECTOR NORM)
    //----------------------------------------------------(6)
    //Compute A_p and the square of its norm
    N_prime( A_p, p_n);
    norm_A_p_sq = A_p.norm_sqr();
    pcout << "l2 Norm A_p_sq: " << norm_A_p_sq  << std::endl;

    //----------------------------------------------------(7)
    //Compute sigma
    //sigma = ||A*r_n||^2 / ||A p_n||^2
    sigma = norm_As_r_sq / norm_A_p_sq; 
    pcout << "sigma     : " << sigma << std::endl;

    //----------------------------------------------------(8)
    //Compute x_n
    //x_n = x_n + sigma p_n
    x_n.sadd(1., sigma, p_n);
    //pcout << "x_n" << std::endl;
    //compute_lambda_norm(x_n);
    //lambda_plotter.write_txt(x_n, "cpp_vector");

    //----------------------------------------------------(9)
    //Compute r_n
    //r_n = r_n - sigma Ap_n
    r_n.sadd(1.,-sigma, A_p);
    //flux_plotter.write_flux_txt(r_n, "cpp_vector");
    //pcout << "r_n computation " << std::endl;
    //compute_product_norm_sq(r_n);

    //----------------------------------------------------(10)
    //Compute As_r
    N_star(As_r, r_n);
    //lambda_plotter.write_txt(As_r, "cpp_vector");
    //As.vmult(As_r, r_n);


    //----------------------------------------------------(11)
    //Compute new norms
    //pcout << "N*(r_new) " << std::endl;
    //norm_As_r_new    = As_r.l2_norm();
    norm_As_r_new_sq = l2_norm_obj.compute_norm_sq(lambda_dof_handler, As_r);
    norm_As_r_new    = sqrt(norm_As_r_new_sq);

    //----------------------------------------------------(12)
    //Compute tau
    tau = norm_As_r_new_sq / norm_As_r_sq; 
    pcout << "tau     : " << tau << std::endl;

    //----------------------------------------------------(13)
    //Compute p_n
    p_n.sadd(tau, 1, As_r);
    //lambda_plotter.write_txt(p_n, "cpp_vector");

    //----------------------------------------------------(14)
    //Updates
    norm_As_r    = norm_As_r_new;
    norm_As_r_sq = norm_As_r_new_sq;

    pcout << "Norm As_r_sq: " << norm_As_r_sq << std::endl;


  }//end_of_CG

  //------------------------------------------------------(15)
  //Update lambda
  //lambda_vector.equ(1, lambda_source, 1, x_n);

  //---------------------------------------------------Plotting
  lambda_vector  = lambda_init;
  lambda_vector += x_n;
  darcy_operator  ();
  stokes_operator ();
  //---
  double final_flux_norm = compute_rho_flux(flux_vector, stokes_vector, darcy_vector);
  pcout << "Final flux norm: " << final_flux_norm << std::endl;
  //cg_stats.get_snapshot(final_flux_norm);
  plotter.save_flux_norm( compute_rho_flux(flux_vector, stokes_vector, darcy_vector) );
  plotter.write_vtk(flux_vector, lambda_vector);
  //---
  //cg_stats.generate_data();
  plotter.generate_data(iteration);
  plotter.save_count(count);
  //---

}//end_cg_lsq
