/** This class is used for the Darcy problem. A similar class
 * can be used for the Darcy problem. One needs to modify
 * the preconditioners since we no longer have a Laplace matrix.
 */
namespace DarcySimpleLinearSolvers
{
  template <class PreconditionerA>
  class ApproxSchurComplement : public Subscriptor
  {
  public:
    ApproxSchurComplement ( const IndexSet                             &partitioning,
                            const TrilinosWrappers::BlockSparseMatrix  &M,
                            const PreconditionerA                      &Apreconditioner)
      :
      tmp1                            (partitioning),
      tmp2                            (partitioning),
      darcy_matrix                    (&M),
      a_preconditioner                (Apreconditioner)
    {}

    void vmult        (TrilinosWrappers::MPI::Vector       &dst,
                const TrilinosWrappers::MPI::Vector        &src) const
    {
      //TrilinosWrappers::MPI::Vector tmp1, tmp2;

      darcy_matrix->block(0,1).vmult (tmp1, src);
      a_preconditioner.vmult (tmp2, tmp1);
      darcy_matrix->block(1,0).vmult (dst, tmp2);
    }//end_v_mult

  private:
    mutable TrilinosWrappers::MPI::Vector tmp1, tmp2;
    const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> darcy_matrix;
    const PreconditionerA  &a_preconditioner;
  };

  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------

  template <class PreconditionerA>
  class SchurComplement : public Subscriptor
  {
  public:
    SchurComplement       ( const IndexSet                             &partitioning,
                            const TrilinosWrappers::BlockSparseMatrix  &M,
                            const PreconditionerA                      &Apreconditioner)
      :
      tmp1                            (partitioning),
      tmp2                            (partitioning),
      darcy_matrix                    (&M),
      a_preconditioner                (Apreconditioner)
    {}

    void vmult        (TrilinosWrappers::MPI::Vector       &dst,
                const TrilinosWrappers::MPI::Vector        &src) const
    {
      //TrilinosWrappers::MPI::Vector tmp1, tmp2;

      darcy_matrix->block(0,1).vmult (tmp1, src);

      SolverControl solver_control(5000, src.l2_norm()*1e-6);
      SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
      solver.solve(darcy_matrix->block(0,0), tmp2, tmp1,
                   a_preconditioner);
      darcy_matrix->block(1,0).vmult (dst, tmp2);
    }//end_v_mult

  private:
    mutable TrilinosWrappers::MPI::Vector tmp1, tmp2;
    const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> darcy_matrix;
    const PreconditionerA  &a_preconditioner;
  };
  //-----------------------------------------------------------------------------

  template <class PreconditionerS, class PreconditionerA, class PreconditionerW>
  class BlockSimpleSchurPreconditioner  : public Subscriptor
  {
  public:
    BlockSimpleSchurPreconditioner (  const PreconditionerS &Spreconditioner,
                                      const PreconditionerA &Apreconditioner,
                                      const PreconditionerW &Wpreconditioner,
                                      const bool            do_solve_W)
      :
      approx_schur_complement         (&Spreconditioner),
      a_preconditioner                (Apreconditioner),
      w_preconditioner                (Wpreconditioner),
      do_solve_W                      (do_solve_W)
    {}

    void vmult (TrilinosWrappers::MPI::BlockVector       &dst,
                const TrilinosWrappers::MPI::BlockVector &src) const
    {
      TrilinosWrappers::MPI::Vector utmp(src.block(0));


      a_preconditioner.vmult (dst.block(0), src.block(0));

      if (do_solve_W == true)
      {
        SolverControl solver_control(5000, src.block(1).l2_norm()*1e-6);

        SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
        //TrilinosWrappers::SolverCG solver(solver_control);

        solver.solve(*approx_schur_complement, dst.block(1), src.block(1),
                     w_preconditioner);
      }
      else
      {
        w_preconditioner.vmult (dst.block(1), src.block(1));
      }


    }//end_v_mult

  private:
    const SmartPointer<const PreconditionerS> approx_schur_complement;
    const PreconditionerA &a_preconditioner;
    const PreconditionerW &w_preconditioner;
    const bool do_solve_W;
  };
}
