template <int dim>
void SD<dim>::darcy_simple_solver ()
{
  TimerOutput::Scope timer_section(computing_timer, "Darcy solver");

  {
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


    PrimitiveVectorMemory<TrilinosWrappers::MPI::BlockVector> mem;
    const DarcySimpleLinearSolvers::ApproxSchurComplement<Precondition_A_Darcy>
          approx_schur_complement (darcy_solution.block(0).locally_owned_elements(),
                                   darcy_matrix, *darcy_Mu_preconditioner);

    unsigned int n_iterations = 0;
    const double solver_tolerance = 1e-8 * darcy_rhs.l2_norm();
    SolverControl solver_control (30, solver_tolerance);


    try
      {
        const DarcySimpleLinearSolvers::BlockSimpleSchurPreconditioner<
              DarcySimpleLinearSolvers::ApproxSchurComplement<Precondition_A_Darcy>,
              Precondition_A_Darcy,
              Precondition_S_Darcy>
              preconditioner (approx_schur_complement,
                              *darcy_Mu_preconditioner, 
                              *darcy_Amg_preconditioner,
                              false);

        SolverGMRES<TrilinosWrappers::MPI::BlockVector>
        solver(solver_control, mem,
               SolverGMRES<TrilinosWrappers::MPI::BlockVector>::
               AdditionalData(30, true));
        solver.solve(darcy_matrix, darcy_solution, darcy_rhs,
                     preconditioner);

        n_iterations = solver_control.last_step();

        pcout << "The GMRES Darcy has terminated" << std::endl;
      }

    catch (SolverControl::NoConvergence)
      {
        pcout << "*********Aggresive solver***************" << std::endl;
        const DarcySimpleLinearSolvers::BlockSimpleSchurPreconditioner<
              DarcySimpleLinearSolvers::ApproxSchurComplement<Precondition_A_Darcy>,
              Precondition_A_Darcy,
              Precondition_S_Darcy>
              preconditioner (approx_schur_complement,
                              *darcy_Mu_preconditioner, 
                              *darcy_Amg_preconditioner,
                              true);

        SolverControl solver_control_refined (darcy_matrix.m(), solver_tolerance);
        SolverGMRES<TrilinosWrappers::MPI::BlockVector>
        solver(solver_control_refined, mem,
               SolverGMRES<TrilinosWrappers::MPI::BlockVector>::
               AdditionalData(50, true));
        solver.solve(darcy_matrix, darcy_solution, darcy_rhs,
                     preconditioner);

        n_iterations = (solver_control.last_step() +
                        solver_control_refined.last_step());
      }


    darcy_constraints.distribute (darcy_solution);

    //darcy_solution.block(1) *= EquationData::pressure_scaling;

    darcy_vector = darcy_solution;
    pcout << n_iterations  << " iterations."
          << std::endl;
  }//end_darcy_block


}//end_darcy_solver

