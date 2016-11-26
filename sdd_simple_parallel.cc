/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2008 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Javier Ruiz Ram\'{i}rez, Clemson University, Dec 2015
 */


// @sect3{Include files}

// As usual, we start by including some well-known files:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/function_bessel.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
//#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_raviart_thomas.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
//#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/fe_field_function.h>

//Use of the timer data class
#include <deal.II/base/timer.h>
//Conditional output stream
#include <deal.II/base/conditional_ostream.h>
// Then we need to include the header file for the sparse direct solver
// UMFPACK,
#include <deal.II/lac/sparse_direct.h>
//and a few iterative solvers.
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/solver_gmres.h>

// This includes the library for the incomplete LU factorization that will be
// used as a preconditioner in 3D:
#include <deal.II/lac/sparse_ilu.h>

//SparsityTools
#include <deal.II/lac/sparsity_tools.h>

//Libraries to work in parallel
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/base/index_set.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

// This is C++:
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <deal.II/base/std_cxx11/shared_ptr.h>

#include "convergence_rates.h"
//#include "darcy_flux_poly.h"
//#include "stokes_exact_solution.h"
//#include "darcy_exact_solution.h"
#include "stokes_exact_solution_2.h"
#include "darcy_exact_solution_2.h"
//#include "flux_function.h"
#include "darcy_K_inverse.h"
#include "flux_function_2.h"
//#include "interfacial_pressure_function.h"
//#include "test_polynomial_function.h"
#include "interfacial_pressure_function_2.h"
//#include "product_space.h"
//#include "plot_data.h"
#include "plot_i_data.h"
#include "map_linker.h"
#include "parallel_map_linker.h"
//#include "generate_template_matrices.h"
//#include "local_averaging.h"
#include "my_l2_norm.h"
#include "mesh_data.h"
//#include "cg_lsq_stats.h"

// As in all programs, the namespace dealii is included:
//namespace SDSPACE
//{
  //namespace SDSPACE initiated
  using namespace dealii;

  //------------------------------CLASS----------------------------------

  template <int dim>
  class SD
  {
  public:
    SD ();
    void run_lambda_stokes_darcy();

  private:
    //-------------------------------------------------Typedefs
    typedef TrilinosWrappers::MPI::BlockVector   MPI_BlockVector;
    //-------------------------------------------------Global_data
    const unsigned int   n_workers;
    const unsigned int   worker_id;
    ConditionalOStream   pcout;
    TimerOutput          computing_timer;
    const unsigned int   wall_id;
    const unsigned int   interface_id;
    const unsigned int   q_rule_increment;
    double               time;
    double               T_initial;
    double               DT;
    unsigned int         iteration;

    //-------------------------------------------------NP_functions
    void refine_all (const unsigned int &);
    void compute_all_norms(const DoFHandler<dim> &dof_handler, 
                           const   Function<dim> &exact_solution,
                           const   MPI_BlockVector &solution,
                           const   unsigned int &degree);
    void update_extrapolation_vectors();
    void update_old_vectors();
    void interpolate_initial_data();
    void initialize_maps();
    void initialize_parallel_maps();
    void print_basic_stats();

    //------------------------------------------N_operators
    void ready_matrices_and_preconditioners();
    void cg_lsq();
    //void cg_lsq( Vector<double> &lambda_source);
    //void cg_lsq(Vector<double> &lambda_target, const Vector<double> &lambda_source);
          //const Vector<double> &lambda_source);
    void N_operator ( Vector<double> &rho_target, const Vector<double> &lambda_source);
    void N_prime( Vector<double> &rho_target, const Vector<double> &lambda_source);
    void N_star(Vector<double> &lambda_target, const Vector<double> &rho_source);

    //-------------------------------------------------Stokes_functions
    void setup_stokes_dofs ();
    void setup_stokes_matrix (const std::vector<IndexSet> &stokes_partitioning,
                              const std::vector<IndexSet> &stokes_relevant_partitioning);
    void setup_stokes_preconditioner (const std::vector<IndexSet> &stokes_partitioning,
                                      const std::vector<IndexSet> &stokes_relevant_partitioning);
    void setup_stokes_zero_bc_matrix (const std::vector<IndexSet> &stokes_partitioning,
                              const std::vector<IndexSet> &stokes_relevant_partitioning);
    void setup_stokes_zero_bc_preconditioner (const std::vector<IndexSet> &stokes_partitioning,
                                      const std::vector<IndexSet> &stokes_relevant_partitioning);

    void build_stokes_preconditioner ();
    void assemble_stokes_preconditioner ();

    void stokes_star (MPI_BlockVector &stokes_target,
                     const Vector<double> &flux_source);
    void stokes_prime (MPI_BlockVector &stokes_target,
                      const Vector<double> &lambda_source);
    void stokes_operator ();
    //void stokes_operator (MPI_BlockVector &stokes_target,
                         //const Vector<double> &lambda_source);
    void stokes_solver ();
    void output_results (const unsigned int refinement_cycle) const;
    void create_stokes_grid();
    void create_stokes_constraints        (IndexSet &stokes_relevant_set);
    void create_stokes_zero_bc_constraints(IndexSet &stokes_relevant_set);
    void compute_all_stokes_norms();

    //-------------------------------------------------Stokes data
    const double         viscosity;
    const unsigned int   stokes_degree;
    bool                 rebuild_stokes_matrix;
    bool                 rebuild_stokes_preconditioner;
    bool                 rebuild_stokes_zero_bc_matrix;
    bool                 rebuild_stokes_zero_bc_preconditioner;

    parallel::distributed::Triangulation<dim>   stokes_domain;
    FESystem<dim>                               stokes_fe;
    DoFHandler<dim>                             stokes_dof_handler;

    ConstraintMatrix     stokes_constraints;
    ConstraintMatrix     stokes_zero_bc_constraints;

    TrilinosWrappers::BlockSparseMatrix       stokes_matrix;
    TrilinosWrappers::BlockSparseMatrix       stokes_zero_bc_matrix;
    TrilinosWrappers::BlockSparseMatrix       stokes_preconditioner_matrix;
    TrilinosWrappers::BlockSparseMatrix       stokes_zero_bc_preconditioner_matrix;
    TrilinosWrappers::MPI::BlockVector        stokes_vector;
    TrilinosWrappers::MPI::BlockVector        stokes_solution;
    TrilinosWrappers::MPI::BlockVector        old_stokes_vector;
    TrilinosWrappers::MPI::BlockVector        stokes_rhs;

    MapLinker<dim> stokes_flux_map;
    MapLinker<dim> stokes_lambda_map;

    ParallelMapLinker<dim-1,dim> parallel_map;

    unsigned int n_u_s;
    unsigned int n_p_s;

    //-------------------------------------------------Darcy functions
    void setup_darcy_dofs ();
    void setup_darcy_matrix ( const std::vector<IndexSet> &darcy_partitioning,
                              const std::vector<IndexSet> &darcy_relevant_partitioning);
    void setup_darcy_preconditioner ( const std::vector<IndexSet> &darcy_partitioning,
                                      const std::vector<IndexSet> &darcy_relevant_partitioning);
    void assemble_darcy_system ();
    void assemble_darcy_preconditioner ();
    void build_darcy_preconditioner ();
    void darcy_star (MPI_BlockVector &darcy_target,
                    const Vector<double> &flux_source);
    void darcy_prime (MPI_BlockVector &darcy_target,
                     const Vector<double> &lambda_source);
    void darcy_operator ();
    //void darcy_operator (MPI_BlockVector &darcy_target,
                        //const Vector<double> &lambda_source);
    void darcy_solver ();
    void darcy_minres_solver ();
    void darcy_fgmres_solver ();
    //void darcy_simple_solver ();
    void compute_darcy_errors () const;
    void create_darcy_grid();
    void compute_all_darcy_norms();
    void create_darcy_constraints ();
    void create_darcy_zero_constraints ();
    void compare_darcy_data();
    //-------------------------------------------------Darcy data
    const unsigned int    darcy_degree;
    bool                  rebuild_darcy_matrix;
    bool                  rebuild_darcy_preconditioner;
    bool                  rebuild_darcy_zero_bc_matrix;
    bool                  rebuild_darcy_zero_bc_preconditioner;

    parallel::distributed::Triangulation<dim>   darcy_domain;
    FESystem<dim>        darcy_fe;
    const double         grad_div;
    DoFHandler<dim>      darcy_dof_handler;


    TrilinosWrappers::BlockSparseMatrix  darcy_matrix;
    TrilinosWrappers::BlockSparseMatrix  darcy_preconditioner_matrix;

    TrilinosWrappers::MPI::BlockVector   darcy_solution;
    TrilinosWrappers::MPI::BlockVector   darcy_rhs;
    TrilinosWrappers::MPI::BlockVector   darcy_vector;

    ConstraintMatrix          darcy_constraints;
    ConstraintMatrix          darcy_preconditioner_constraints;

    MapLinker<dim> darcy_flux_map;
    MapLinker<dim> darcy_lambda_map;

    unsigned int n_u_d;
    unsigned int n_p_d;

    KInverse<dim>            k_inverse;

    //-------------------------------------------------Q_rules
    const unsigned int    q_rule_degree;
    const QGauss<dim>     q_rule;
    const QGauss<dim-1>   face_q_rule;
    const unsigned int    n_q_points;
    const unsigned int    n_face_q_points;
    //-------------------------------------------------L2_norm
    L2_norm<dim-1, dim> l2_norm_obj;
    //-------------------------------------------------CGStats
    //CGStats<dim-1,dim,Vector<double> > cg_stats;
    //-------------------------------------------------Lambda functions
    void create_lambda_grid();
    void create_lambda_constraints();
    void setup_lambda_dofs ();
    void assemble_lambda_system();
    //void lambda_connection();
    void interpolate_interfacial_pressure();
    void test_lambda_connection();
    void test_flux_connection();
    void solve_lambda();
    //void connect_stokes_and_darcy_to_lambda();
    void extract_stokes_and_darcy_interface_values();

    //-------------------------------------------------Lambda data
    Triangulation<dim-1,dim>   lambda_domain;

    const unsigned int       lambda_degree;
    //double                   lambda_cell_measure;

    FE_Q<dim-1,dim>          lambda_fe;
    DoFHandler<dim-1,dim>    lambda_dof_handler;

    TrilinosWrappers::MPI::Vector  lambda_vector;

    IndexSet                 lambda_locally_owned_dofs; 
    IndexSet                 lambda_locally_relevant_dofs;


    //-------------------------------------------------Flux_functions
    double compute_rho_flux_norm(const Vector<double> &rho_source, bool print);
    void   create_flux_grid();
    void   setup_flux_dofs();
    void   flux_connection();
    double compute_rho_flux(        
                                    Vector<double> &flux_target,
                         const MPI_BlockVector &stokes_source, 
                         const MPI_BlockVector &darcy_source,
                         bool  print = false );
    void   compute_flux_function(       
                                      Vector<double> &lambda_target,
                           const MPI_BlockVector &stokes_source, 
                           const MPI_BlockVector &darcy_source);

    //-------------------------------------------------Flux_data
    Triangulation<dim-1,dim>   flux_domain;

    const unsigned int    flux_refinements;
    FE_DGQ<dim-1,dim>     flux_fe;
    DoFHandler<dim-1,dim> flux_dof_handler;

    TrilinosWrappers::MPI::Vector  flux_vector;

    IndexSet              flux_locally_owned_dofs; 
    IndexSet              flux_locally_relevant_dofs;

    const double          sqrt_delta;
    const double          cg_eps;
    const bool            darcy_bc_are_enforced_weakly;
    double                flux_cell_measure;
    double                sqrt_flux_cell_measure;
    double                i_sqrt_flux_cell_measure;

    //------------------------------------------------Porosity data
    const unsigned int                   porosity_degree;
    FE_Q<dim>                            porosity_fe;
    DoFHandler<dim>                      porosity_dof_handler;
    ConstraintMatrix                     porosity_constraints;

    TrilinosWrappers::SparseMatrix       porosity_matrix;

    TrilinosWrappers::MPI::Vector        porosity_solution;
    TrilinosWrappers::MPI::Vector        old_porosity_solution;
    TrilinosWrappers::MPI::Vector        old_old_porosity_solution;
    TrilinosWrappers::MPI::Vector        porosity_rhs;

    IndexSet                             porosity_locally_owned_dofs; 
    IndexSet                             porosity_locally_relevant_dofs;

    //------------------------------------------------Porosity functions
    void setup_porosity_matrix();

    //------------------------------------------------Typedef
    typedef TrilinosWrappers::PreconditionJacobi Precondition_A_Darcy;
    typedef TrilinosWrappers::PreconditionAMG    Precondition_S_Darcy;
    //typedef PreconditionIdentity                 Precondition_S_Darcy;
    typedef TrilinosWrappers::PreconditionAMG    Precondition_A_Stokes;
    typedef TrilinosWrappers::PreconditionJacobi Precondition_S_Stokes;

    //------------------------------------------------Preconditioners
    std_cxx11::shared_ptr<Precondition_S_Darcy>    darcy_Amg_preconditioner;
    std_cxx11::shared_ptr<Precondition_A_Darcy>    darcy_Mu_preconditioner;

    std_cxx11::shared_ptr<Precondition_A_Stokes>   stokes_Amg_preconditioner;
    std_cxx11::shared_ptr<Precondition_S_Stokes>   stokes_Mp_preconditioner;

    std_cxx11::shared_ptr<Precondition_A_Stokes>   stokes_zero_bc_Amg_preconditioner;
    std_cxx11::shared_ptr<Precondition_S_Stokes>   stokes_zero_bc_Mp_preconditioner;

    PreconditionIdentity                           identity_preconditioner;
    //------------------------------------------------External functions
    InterfacialPressureFunction2<dim> ip_function;
    FluxFunction2<dim>                flux_function;
    StokesExactSolution2<dim>         stokes_exact_solution;
    DarcyExactSolution2<dim>          darcy_exact_solution;
    //------------------------------------------------
    //InterfacialPressureFunction<dim> ip_function;
    //FluxFunction<dim>                flux_function;
    //StokesExactSolution<dim>         stokes_exact_solution;
    //DarcyExactSolution<dim>          darcy_exact_solution;
    //------------------------------------------------Convergence rates
    ConvergenceRates<dim, StokesExactSolution2<dim>, MPI_BlockVector> cr_stokes;
    ConvergenceRates<dim, DarcyExactSolution2<dim>, MPI_BlockVector>  cr_darcy;

  };

  template <int dim>
  SD<dim>::SD ()
    :
    //------------------------------------------Global
    n_workers (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
    worker_id (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
    pcout (std::cout, worker_id == 0),
    computing_timer(pcout, TimerOutput::summary, TimerOutput::wall_times),
    wall_id (1),
    interface_id (0),
    q_rule_increment (3),
    //------------------------------------------Stokes
    viscosity          (0.5),
    stokes_degree      (1),
    stokes_domain      (MPI_COMM_WORLD, 
                       typename Triangulation<dim>::MeshSmoothing
                       (Triangulation<dim>::smoothing_on_refinement |
                        Triangulation<dim>::smoothing_on_coarsening)),
        //Triangulation<dim>::maximum_smoothing),
    stokes_fe          (FE_Q<dim>(stokes_degree+1), dim,
                        FE_Q<dim>(stokes_degree), 1),
    stokes_dof_handler (stokes_domain),
    //------------------------------------------Darcy
    darcy_degree      (1),
    darcy_domain      (MPI_COMM_WORLD, 
                       typename Triangulation<dim>::MeshSmoothing
                       (Triangulation<dim>::smoothing_on_refinement |
                       Triangulation<dim>::smoothing_on_coarsening)),
        //Triangulation<dim>::maximum_smoothing),
    darcy_fe           (FE_Q<dim>(darcy_degree+1), dim,
                        FE_Q<dim>(darcy_degree), 1),
    grad_div          (0),
    //grad_div          (500),
    darcy_dof_handler (darcy_domain),
    //------------------------------------------Q_rule
    q_rule_degree     (std::max(stokes_degree, darcy_degree) + q_rule_increment),
    q_rule            (q_rule_degree),
    face_q_rule       (q_rule_degree),
    n_q_points        (q_rule.size()),
    n_face_q_points   (face_q_rule.size()),
    //------------------------------------------L2_norm
    l2_norm_obj       (q_rule_degree),
    //------------------------------------------Lambda
    //lambda_domain      (MPI_COMM_WORLD),
    lambda_degree      (2),
    lambda_fe          (lambda_degree),
    //lambda_dof_handler (lambda_domain),
    //------------------------------------------Flux
    //flux_domain      (MPI_COMM_WORLD),
    flux_refinements (2),
    flux_fe          (0),
    //flux_dof_handler (flux_domain),
    //------------------------------------------
    sqrt_delta       (0),
    cg_eps           (1e-6),
    darcy_bc_are_enforced_weakly (true),
    //------------------------------------------Porosity
    porosity_degree      (1),
    porosity_fe          (porosity_degree),
    porosity_dof_handler (darcy_domain)
    //------------------------------------------
    //End of initialization
    //------------------------------------------
  {
    //-------------------------------------------------
    pcout << "Darcy  fe: " << darcy_fe.get_name() << std::endl;
    pcout << "Stokes fe: " << stokes_fe.get_name() << std::endl;
    //pcout << "Lambda fe: " << lambda_fe.get_name() << std::endl;
    //pcout << "Flux   fe: " << flux_fe.get_name() << std::endl;
    //-------------------------------------------------
    pcout << "Grad div         : " << grad_div << std::endl;
    pcout << "CG LSQ eps       : " << cg_eps << std::endl;
    pcout << "Darcy BC weak    : " << (darcy_bc_are_enforced_weakly ? "True" : "False") << std::endl;
    pcout << "Flux refinements : " << flux_refinements << std::endl;
    //-------------------------------------------------
    //std::cout << "Worker id: " << worker_id << std::endl;
    pcout << "# workers        : " << n_workers << std::endl;
  }
  //--------------------------------------------- 
  //#include "enforce_flux_continuity.w"
  //#include "darcy_playground.w"
  //#include "test_adjoint.w"
  //#include "darcy_simple_linear_solvers.w"
  //#include "darcy_simple_solver.w"
  #include "ready_matrices_and_preconditioners.w"

  #include "darcy_linear_solvers.w"
  #include "darcy_block_schur_preconditioner.w"
  #include "darcy_solver.w"
  #include "darcy_minres_solver.w"
  #include "darcy_fgmres_solver.w"

  #include "stokes_linear_solvers.w"
  #include "stokes_solver.w"
  //---------------------------------------------
  //#include "domain_connection.w"
  //#include "flux_connection.w"
  //#include "compute_flux.w"
  //#include "flux_initialization.w"
  //---------------------------------------------
  //#include "N_operators.w"
  //#include "cg_lsq.w"
  //---------------------------------------------
  #include "stokes_initialization.w"
  //#include "stokes_assembly.w"
  //#include "stokes_operators.w"
  //#include "stokes_stationary.w"
  #include "assemble_stokes_preconditioner.w"
  #include "build_stokes_preconditioner.w"
  #include "stokes_operator.w"
  //#include "stokes_stationary2.w"
  //---------------------------------------------
  #include "darcy_initialization.w"
  //#include "darcy_assembly.w"
  //#include "darcy_operators.w"
  #include "assemble_darcy_preconditioner.w"
  #include "assemble_darcy_system.w"
  #include "build_darcy_preconditioner.w"
  #include "darcy_stationary.w"
  //-------------------------------------------------------------
  //#include "lambda_initialization.w"
  //-------------------------------------------------------------
  //#include "error_norms.w"
  #include "print_basic_stats.w"
  //#include "plot_data_dealii.w"
  //-------------------------------------------------------------
  //#include "test_flux_communication.w"
  //#include "initialize_maps.w"
  #include "initialize_parallel_maps.w"
  //-------------------------------------------------------------
  //#include "run_lambda_stokes_darcy.w"
  #include "run_stationary.w"
  //-------------------------------------------------------------
  //#include "output_results.w"

//}//end_namespace

// @sect3{The <code>main</code> function}

// The main function is the same as in step-20. We pass the element degree as
// a parameter and choose the space dimension at the well-known template slot.
int main (int argc, char* argv[] )
{
  using namespace dealii;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);
  // This program can only be run in serial. Otherwise, throw an exception.
  //AssertThrow( Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1,
               //ExcMessage("This program can only be run in serial, use ./filename"));

  try
    {
      //using namespace SDSPACE;

      deallog.depth_console (0);

      SD<3> flow_problem;
      //flow_problem.run_simple ();
      //flow_problem.run ();
      flow_problem.run_lambda_stokes_darcy ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
