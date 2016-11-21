template <int dim>
void SD<dim>::stokes_solver ()
{
  TimerOutput::Scope timer_section(computing_timer, "Stokes solver");

  {
    pcout << "   Solving Stokes system... " << std::endl;

    //Optional
    stokes_solution = stokes_vector;

    //stokes_solution.block(1) /= EquationData::pressure_scaling;

    const unsigned int
    start = (stokes_solution.block(0).size() +
             stokes_solution.block(1).local_range().first),
    end   = (stokes_solution.block(0).size() +
             stokes_solution.block(1).local_range().second);

    for (unsigned int i=start; i<end; ++i)
    if (stokes_constraints.is_constrained (i))
        stokes_solution(i) = 0;


    PrimitiveVectorMemory<TrilinosWrappers::MPI::BlockVector> mem;

    unsigned int n_iterations = 0;
    const double solver_tolerance = 1e-8 * stokes_rhs.l2_norm();
    SolverControl solver_control (30, solver_tolerance);

    try
      {
        const StokesLinearSolvers::BlockSchurPreconditioner<Precondition_A_Stokes,
              Precondition_S_Stokes>
              preconditioner (stokes_matrix, stokes_preconditioner_matrix,
                              *stokes_Mp_preconditioner, *stokes_Amg_preconditioner,
                              false);

        SolverFGMRES<TrilinosWrappers::MPI::BlockVector>
        solver(solver_control, mem,
               SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::
               AdditionalData(30, true));
        solver.solve(stokes_matrix, stokes_solution, stokes_rhs,
                     preconditioner);

        n_iterations = solver_control.last_step();
        pcout << "The FGMRES Stokes has terminated" << std::endl;
      }

    catch (SolverControl::NoConvergence)
      {
        pcout << "*********Aggresive solver***************" << std::endl;
        const StokesLinearSolvers::BlockSchurPreconditioner<Precondition_A_Stokes,
              Precondition_S_Stokes>
              preconditioner (stokes_matrix, stokes_preconditioner_matrix,
                              *stokes_Mp_preconditioner, *stokes_Amg_preconditioner,
                              true);

        SolverControl solver_control_refined (stokes_matrix.m(), solver_tolerance);
        SolverFGMRES<TrilinosWrappers::MPI::BlockVector>
        solver(solver_control_refined, mem,
               SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::
               AdditionalData(50, true));
        solver.solve(stokes_matrix, stokes_solution, stokes_rhs,
                     preconditioner);

        n_iterations = (solver_control.last_step() +
                        solver_control_refined.last_step());
      }


    stokes_constraints.distribute (stokes_solution);

    //stokes_solution.block(1) *= EquationData::pressure_scaling;

    stokes_vector = stokes_solution;
    pcout << n_iterations  << " iterations."
          << std::endl;
  }//end_stokes_block


}//end_stokes_solver

