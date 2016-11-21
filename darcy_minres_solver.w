template <int dim>
void SD<dim>::darcy_minres_solver ()
{
    TimerOutput::Scope timer_section(computing_timer, "Darcy solver");

    pcout << "   Solving Darcy system... " << std::endl;

    //Optional
    darcy_solution = darcy_vector;

    //darcy_solution.block(1) /= EquationData::pressure_scaling;

    //const unsigned int
    //start = (darcy_solution.block(0).size() +
             //darcy_solution.block(1).local_range().first),
    //end   = (darcy_solution.block(0).size() +
             //darcy_solution.block(1).local_range().second);
//
    //for (unsigned int i=start; i<end; ++i)
    //if (darcy_constraints.is_constrained (i))
        //darcy_solution(i) = 0;

    const DarcyLinearSolvers::ApproxSchurComplement<Precondition_A_Darcy>
          approx_schur_complement (darcy_solution.block(0).locally_owned_elements(),
                                   darcy_matrix, *darcy_Mu_preconditioner);

    const DarcyLinearSolvers::InvMatrix<
          DarcyLinearSolvers::ApproxSchurComplement<Precondition_A_Darcy>,
          PreconditionIdentity >
          inv_approx_schur_complement (approx_schur_complement, identity_preconditioner);

    const DarcyLinearSolvers::BlockDiagPreconditioner<
                Precondition_A_Darcy,
                DarcyLinearSolvers::InvMatrix<
                DarcyLinearSolvers::ApproxSchurComplement<Precondition_A_Darcy>,
                PreconditionIdentity >
                > 
          preconditioner (*darcy_Mu_preconditioner, inv_approx_schur_complement);

    SolverControl solver_control (darcy_matrix.m(), 1e-10*darcy_rhs.l2_norm());
    SolverMinRes<MPI_BlockVector> solver (solver_control);

    darcy_constraints.set_zero (darcy_solution);

    solver.solve(darcy_matrix, darcy_solution, darcy_rhs, preconditioner);
    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;

    darcy_constraints.distribute (darcy_solution);

    darcy_vector = darcy_solution;

    //pcout << n_iterations  << " iterations." << std::endl;


}//end_darcy_solver
