#ifndef CONVERGENCE_RATES_H
#define CONVERGENCE_RATES_H
//------------------------------------------------------------Dependencies
#include <map>
#include <deal.II/lac/block_vector.h>
#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/base/quadrature_lib.h>
//------------------------------------------------------------
DEAL_II_NAMESPACE_OPEN
//------------------------------------------------------------
/**This class computes the convergence rates for stationary and
 * time dependent problems. The functionality for each case is different.
 */
//------------------------------------------------------------
template<int dim, class FF, typename V = BlockVector<double> >
class ConvergenceRates
{
  public:
   ConvergenceRates(); 

   double p_domain;
   double measure;

   void add_field(const std::string &);
   void add_time_space_error(const std::vector<std::vector<double> > &, 
                             const double &, const double &);
   //void compute_stationary_rates (std::vector<std::vector<double> > &iEV);
   void compute_stationary_rates ();
   void compute_rates();
   void save_h_value (const double &h_param);
   unsigned int size();
   void clean();
   void compute_combined_mean_p(const double &n_p_domain, const double &n_measure);
   void compute_norms(const DoFHandler<dim> &dh,
                      const QGauss<dim>     &q_rule,
                      const V               &fem_solution,
                      const FF              &exact_solution);
   void compute_pressure_integral
                     (const DoFHandler<dim> &dh,
                      const QGauss<dim>     &q_rule,
                      const V               &fem_solution);
  private:
   std::vector<std::vector<double> > EV;
   std::vector<std::vector<double> > NV;
   double p_true_mean;
   double p_mean;
   enum Norms {gul2, ul2, pl2, udiv, uhdiv, uh1};
   std::map<std::string, Norms> norm_map;
   std::vector<std::string>          fields;
   std::vector<std::vector<double> > errors;
   std::vector<double>               dt_vec;
   std::vector<double>               h_vec;
   std::vector<unsigned int>         iter_vec;
   void update_norms();
};
//-------------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
//------------------------------------------------------------
#endif
