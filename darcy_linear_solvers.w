/** This class is used for the Darcy problem. A similar class
 * can be used for the Darcy problem. One needs to modify
 * the preconditioners since we no longer have a Laplace matrix.
 */
namespace DarcyLinearSolvers
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

  template <class MATRIX, class PreconditionerA>
  class InvMatrix : public Subscriptor
  {
  public:
    InvMatrix ( 
                            const MATRIX            &M,
                            const PreconditionerA   &Apreconditioner)
      :
      matrix                          (&M),
      a_preconditioner                (Apreconditioner)
    {}

    void print_summary() const
    {
      std::cout << "Inv Approx Schur Complement was called" << std::endl;
      std::cout << ">>> itertations: " << n_iterations << std::endl;
    }

    void vmult        (TrilinosWrappers::MPI::Vector       &dst,
                 const TrilinosWrappers::MPI::Vector        &src) const
    {
      SolverControl solver_control(5000, src.l2_norm()*1e-2);
      SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
      solver.solve(*matrix, dst, src, a_preconditioner);

      n_iterations = solver_control.last_step();
      print_summary();
    }//end_v_mult


  private:
    const SmartPointer<const MATRIX> matrix;
    mutable unsigned int           n_iterations;
    //const MATRIX           &matrix;
    const PreconditionerA  &a_preconditioner;
  };
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

    void print_summary() const
    {
      std::cout << "Schur Complement was called (A^{-1})" << std::endl;
      std::cout << ">>> itertations: " << n_iterations << std::endl;
    }

    void vmult        (TrilinosWrappers::MPI::Vector       &dst,
                const TrilinosWrappers::MPI::Vector        &src) const
    {
      //TrilinosWrappers::MPI::Vector tmp1, tmp2;

      darcy_matrix->block(0,1).vmult (tmp1, src);

      SolverControl solver_control(5000, tmp1.l2_norm()*1e-10);
      SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
      solver.solve(darcy_matrix->block(0,0), tmp2, tmp1, a_preconditioner);
      darcy_matrix->block(1,0).vmult (dst, tmp2);

      n_iterations = solver_control.last_step();
      print_summary();
    }//end_v_mult


  private:
    mutable unsigned int                  n_iterations;
    mutable TrilinosWrappers::MPI::Vector tmp1, tmp2;
    const SmartPointer<const TrilinosWrappers::BlockSparseMatrix>     
                                          darcy_matrix;
    const PreconditionerA                 &a_preconditioner;
  };
  //-----------------------------------------------------------------------------

  template <class PreconditionerA, class PreconditionerS>
  class BlockDiagPreconditioner  : public Subscriptor
  {
  public:
    BlockDiagPreconditioner (         const PreconditionerA &Apreconditioner,
                                      const PreconditionerS &Spreconditioner,
                                      const bool            do_solve_W = false)
      :
      a_preconditioner                (Apreconditioner),
      s_preconditioner                (Spreconditioner),
      do_solve_W                      (do_solve_W)
    {}

    void vmult (TrilinosWrappers::MPI::BlockVector       &dst,
                const TrilinosWrappers::MPI::BlockVector &src) const
    {
      a_preconditioner.vmult (dst.block(0), src.block(0));
      s_preconditioner.vmult (dst.block(1), src.block(1));

    }//end_v_mult

  private:
    const PreconditionerA &a_preconditioner;
    const PreconditionerS &s_preconditioner;
    const bool do_solve_W;
  };
}
