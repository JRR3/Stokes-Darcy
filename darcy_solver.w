template <int dim>
void SD<dim>::darcy_solver ()
{
    TimerOutput::Scope timer_section(computing_timer, "Darcy solver");

    pcout << "   Solving Darcy system... " << std::endl;

    //Optional
    darcy_solution = darcy_vector;

    //darcy_solution.block(1) /= EquationData::pressure_scaling;

    const unsigned int
    start = (darcy_solution.block(0).size() +
             darcy_solution.block(1).local_range().first),
    end   = (darcy_solution.block(0).size() +
             darcy_solution.block(1).local_range().second);

    for (unsigned int i=start; i<end; ++i)
    if (darcy_constraints.is_constrained (i))
        darcy_solution(i) = 0;

    TrilinosWrappers::MPI::Vector utmp      (darcy_rhs.block(0));
    TrilinosWrappers::MPI::Vector schur_rhs (darcy_rhs.block(1));
    const double tol = 1e-8;


    SolverControl solver_control(5000, darcy_rhs.block(0).l2_norm()*tol);
    SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
    solver.solve(darcy_matrix.block(0,0), utmp, darcy_rhs.block(0), *darcy_Mu_preconditioner);

    unsigned int n_iterations = solver_control.last_step();
    pcout << "Iterations (RHS): " << n_iterations << std::endl;


    darcy_matrix.block(1,0).vmult (schur_rhs, utmp);

    schur_rhs -= darcy_rhs.block(1);

    //BA^{-1}F - G
    
    //CHECK: Preconditioner->Jacobi
    const DarcyLinearSolvers::ApproxSchurComplement<Precondition_A_Darcy>
          approx_schur_complement (darcy_solution.block(0).locally_owned_elements(),
                                   darcy_matrix, *darcy_Mu_preconditioner);

    const DarcyLinearSolvers::InvMatrix<
          DarcyLinearSolvers::ApproxSchurComplement<Precondition_A_Darcy>,
          PreconditionIdentity >
          inv_approx_schur_complement (approx_schur_complement, identity_preconditioner);

    //CHECK: Preconditioner->Identity
    const DarcyLinearSolvers::SchurComplement<PreconditionIdentity>
          schur_complement (darcy_solution.block(0).locally_owned_elements(),
                                   darcy_matrix, identity_preconditioner);

    solver.solve(schur_complement, darcy_solution.block(1), schur_rhs, inv_approx_schur_complement);

    pcout << "*******Iterations (Inv Schur): " << solver_control.last_step() << std::endl;

    darcy_matrix.block(0,1).vmult (utmp, darcy_solution.block(1));
    utmp *= -1;
    utmp += darcy_rhs.block(0);

    solver.solve(darcy_matrix.block(0,0), darcy_solution.block(0), utmp, *darcy_Mu_preconditioner);

    pcout << "Iterations (U): " << solver_control.last_step() << std::endl;

    darcy_constraints.distribute (darcy_solution);

    darcy_vector = darcy_solution;

    //pcout << n_iterations  << " iterations." << std::endl;


}//end_darcy_solver

