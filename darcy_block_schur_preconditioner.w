/** This class is used for the Darcy problem. A similar class
 * can be used for the Darcy problem. One needs to modify
 * the preconditioners since we no longer have a Laplace matrix.
 */
namespace DarcyLinearSolvers
{
  template <class PreconditionerA, class PreconditionerW>
  class BlockSchurPreconditioner : public Subscriptor
  {
  public:
    BlockSchurPreconditioner (const TrilinosWrappers::BlockSparseMatrix  &S,
                              const TrilinosWrappers::BlockSparseMatrix  &Spre,
                              const PreconditionerA                      &Apreconditioner,
                              const PreconditionerW                      &Wpreconditioner,
                              const bool                                  do_solve_W)
      :
      darcy_matrix                    (&S),
      darcy_preconditioner_matrix     (&Spre),
      a_preconditioner                (Apreconditioner),
      w_preconditioner                (Wpreconditioner),
      do_solve_W                      (do_solve_W)
    {}

    void vmult (TrilinosWrappers::MPI::BlockVector       &dst,
                const TrilinosWrappers::MPI::BlockVector &src) const
    {
      TrilinosWrappers::MPI::Vector utmp(src.block(0));

      if (do_solve_W == true)
        {
          SolverControl solver_control(5000, 1e-2 * src.block(1).l2_norm());
          TrilinosWrappers::SolverCG solver(solver_control);
          solver.solve(darcy_preconditioner_matrix->block(1,1),
                       dst.block(1), src.block(1),
                       w_preconditioner);
        }
      else
        w_preconditioner.vmult (dst.block(1), src.block(1));

      dst.block(1) *= -1.0;


      {
        darcy_matrix->block(0,1).vmult(utmp, dst.block(1));
        utmp*=-1.0;
        utmp.add(src.block(0));
      }

      {
        SolverControl solver_control(5000, utmp.l2_norm()*1e-6);

        SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);

        solver.solve(darcy_matrix->block(0,0), dst.block(0), utmp,
                     a_preconditioner);

      }
    }//end_v_mult

  private:
    const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> darcy_matrix;
    const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> darcy_preconditioner_matrix;
    const PreconditionerA  &a_preconditioner;
    const PreconditionerW  &w_preconditioner;
    const bool do_solve_W;
  };
}

