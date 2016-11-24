#include "convergence_rates.h"
#include "stokes_exact_solution_2.h"
#include "darcy_exact_solution_2.h"
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
//------------------------------------------------------------
DEAL_II_NAMESPACE_OPEN
//------------------------------------------------------------
template<int dim, class FF, typename V>
ConvergenceRates<dim, FF, V>::ConvergenceRates():
  EV         (4),
  p_true_mean( exp(1)*(cos(1)-1)/2 + sin(1)*(exp(2)-1)/2 + 1./3 ),
  p_mean     (0),
  norm_map   ( {{"uh1",uh1}, {"ul2",ul2}, {"pl2",pl2},
               {"uhdiv",uhdiv}, {"udiv",udiv}, {"gul2",gul2}})
{
}
//-------------------------------------------------------------------
template<int dim, class FF, typename V>
unsigned int ConvergenceRates<dim, FF, V>::size() 
{
  return fields.size();
}
//-------------------------------------------------------------------
template<int dim, class FF, typename V>
void ConvergenceRates<dim, FF, V>::clean() 
{
  for(unsigned int i = 0; i < fields.size(); ++i)
  {
    errors[i].clear();
  }
}
//-------------------------------------------------------------------
template<int dim, class FF, typename V>
void ConvergenceRates<dim, FF, V>::add_field(const std::string &name)
{
  fields.push_back(name);
  errors.push_back(std::vector<double>());
  NV.push_back(std::vector<double>());
}
//-------------------------------------------------------------------
template<int dim, class FF, typename V>
void ConvergenceRates<dim, FF, V>::add_time_space_error
                       (const std::vector<std::vector<double> > &vec, 
                        const double &DT,
                        const double &h_param)
{
  for(unsigned int f = 0; f < fields.size(); ++f)
  {
    double s = 0;
    for(unsigned int i = 0; i < vec[f].size(); ++i)
    {
      s += vec[f][i] * vec[f][i];
    }
    errors[f].push_back(sqrt(s*DT));
  }
    dt_vec.push_back(DT);
    h_vec.push_back(h_param);
}
//-------------------------------------------------------------------
template<int dim, class FF, typename V>
void ConvergenceRates<dim, FF, V>::save_h_value (const double &h_param)
{
  h_vec.push_back(h_param);
}
//-------------------------------------------------------------------
template<int dim, class FF, typename V>
void ConvergenceRates<dim, FF, V>::compute_stationary_rates()
//void ConvergenceRates<dim, FF, V>::compute_stationary_rates
                      //(std::vector<std::vector<double> > &EV)
{
  update_norms();

  std::vector<std::vector<double> > log_ratio_vec (fields.size(), std::vector<double>());

  for(unsigned int f = 0; f < fields.size(); ++f)
    log_ratio_vec[f].push_back(-100);

  std::cout << std::endl;
  std::cout << "Condensed view ..." << std::endl;

  for(unsigned int f = 0; f < fields.size(); ++f)
  {
    std::cout << "Rates for : " << fields[f] << std::endl;
    for(unsigned int k = 1; k < NV[f].size(); ++k)
    {
      double ratio     = NV[f][k-1]/NV[f][k];
      double log_ratio = log2(ratio);
      log_ratio_vec[f].push_back(log_ratio);
      std::cout << log_ratio << std::endl;
    }
  }//end_for_condensed

  //---------------------------------------

  std::cout << std::endl;
  std::cout << "Latex view ..." << std::endl;
  std::cout << "\\hline" << std::endl;
  for(unsigned int k = 0; k < h_vec.size(); ++k)
  {
    std::cout << std::fixed;
    std::cout.precision(3);
    std::cout << h_vec[k] ;

    //std::cout << std::scientific;
    //std::cout.precision(2);
    //std::cout <<  " & " << dt_vec[k] ;

    for(unsigned int f = 0; f < fields.size(); ++f)
    {
      std::cout << std::scientific;
      std::cout.precision(2);
      std::cout << " & " << NV[f][k]  ;
      std::cout << std::fixed;
      if(k > 0)
        std::cout << " & " << log_ratio_vec[f][k]  ;
      else
        std::cout <<  " & - " ;
    }
        std::cout <<  " \\\\ " << std::endl;
        std::cout <<  "\\hline" << std::endl;
  }
}//end_compute_stat_rates
//-------------------------------------------------------------------
template<int dim, class FF, typename V>
void ConvergenceRates<dim, FF, V>::compute_rates()
{
  std::vector<std::vector<double> > log_ratio_vec (fields.size(), std::vector<double>());

  for(unsigned int f = 0; f < fields.size(); ++f)
    log_ratio_vec[f].push_back(-100);

  std::cout << std::endl;
  std::cout << "Condensed view ..." << std::endl;

  for(unsigned int f = 0; f < fields.size(); ++f)
  {
    std::cout << "Rates for : " << fields[f] << std::endl;
    for(unsigned int k = 1; k < errors[f].size(); ++k)
    {
      double ratio     = errors[f][k-1]/errors[f][k];
      double log_ratio = log2(ratio);
      log_ratio_vec[f].push_back(log_ratio);
      std::cout << log_ratio << std::endl;
    }
  }//end_for_condensed

  //---------------------------------------

  std::cout << std::endl;
  std::cout << "Latex view ..." << std::endl;
  std::cout << "\\hline" << std::endl;
  for(unsigned int k = 0; k < h_vec.size(); ++k)
  {
    std::cout << std::fixed;
    std::cout.precision(3);
    std::cout << h_vec[k] ;

    std::cout << std::scientific;
    std::cout.precision(2);
    std::cout <<  " & " << dt_vec[k] ;

    for(unsigned int f = 0; f < fields.size(); ++f)
    {
      std::cout << std::scientific;
      std::cout.precision(2);
      std::cout << " & " << errors[f][k]  ;
      std::cout << std::fixed;
      if(k > 0)
        std::cout << " & " << log_ratio_vec[f][k]  ;
      else
        std::cout <<  " & - " ;
    }
        std::cout <<  " \\\\ " << std::endl;
        std::cout <<  "\\hline" << std::endl;
  }
}
//-------------------------------------------------------------------------
template<int dim, class FF, typename V>
void ConvergenceRates<dim, FF, V>::compute_pressure_integral
                                            (const DoFHandler<dim> &dh,
                                             const QGauss<dim>     &q_rule,
                                             const V               &fem_solution)
{
  typename DoFHandler<dim>::active_cell_iterator 
    cell = dh.begin_active(),
   ecell = dh.end();
  unsigned int n_q_points = q_rule.size();
  FEValuesExtractors::Scalar pressure(dim);
  //-------------------------------------------------------
  FEValues<dim> fe_values(dh.get_fe(), q_rule,
                          update_values    |
                          update_JxW_values);

  std::vector<double> p_values (n_q_points);
  p_domain   = 0;
  measure    = 0;

  for(; cell != ecell; ++cell)
  {
    fe_values.reinit(cell);  
    fe_values[pressure].get_function_values     (fem_solution, p_values);
    for(unsigned int q = 0; q < n_q_points; ++q)
    {
      p_domain    += p_values[q] * fe_values.JxW(q);
      measure     += fe_values.JxW(q);
    }
  }

  p_mean = p_domain / measure;

}
//-------------------------------------------------------------------------
template<int dim, class FF, typename V>
void ConvergenceRates<dim, FF, V>::compute_combined_mean_p(const double &n_p_domain, const double &n_measure)
{
  p_mean = (p_domain + n_p_domain) / (measure + n_measure);
}
//-------------------------------------------------------------------------
template<int dim, class FF, typename V>
void ConvergenceRates<dim, FF, V>::compute_norms(const DoFHandler<dim> &dh,
                                             const QGauss<dim>     &q_rule,
                                             const V               &fem_solution,
                                             const FF              &exact_solution)
{
  std::cout << ">>>Compute  errors " << std::endl;

  typename DoFHandler<dim>::active_cell_iterator 
    cell = dh.begin_active(),
   ecell = dh.end();
  unsigned int n_q_points = q_rule.size();
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure(dim);
  //exact_solution.set_time(time);
  //-------------------------------------------------------
  FEValues<dim> fe_values(dh.get_fe(), q_rule,
                          update_values    |
                          update_quadrature_points |
                          update_gradients |
                          update_JxW_values);
  std::vector<double>             p_values          (n_q_points);
  std::vector<double>             p_true_values     (n_q_points);
  std::vector<Tensor<1,dim> >     u_values          (n_q_points);
  std::vector<Tensor<1,dim> >     u_true_values     (n_q_points);
  std::vector<Tensor<2,dim> >     grad_u_values     (n_q_points);
  std::vector<Tensor<2,dim> >     grad_u_true_values(n_q_points);
  std::vector<Point<dim> >        q_points          (n_q_points);
  double                          grad_error = 0;
  double                          u_error    = 0;
  double                          div_u      = 0;
  double                          p_error    = 0;
  //double                          p_mean     = 0;
  //double                          p_true_mean= 0;
  //double                          measure    = 0;
  //for(; cell != ecell; ++cell)
  //{
    //fe_values.reinit(cell);  
    //fe_values[pressure].get_function_values     (fem_solution, p_values);
    //exact_solution.get_p_list (q_points, p_true_values);
    //for(unsigned int q = 0; q < n_q_points; ++q)
    //{
      //p_mean      += p_values[q] * fe_values.JxW(q);
      //p_true_mean += p_true_values[q] * fe_values.JxW(q);
      //measure     += fe_values.JxW(q);
    //}
  //}
  //p_mean /= measure;
  //p_true_mean /= measure;
  ////-------------------------------------
  //cell = dh.begin_active();
  ////-------------------------------------
  for(; cell != ecell; ++cell)
  {
    fe_values.reinit(cell);  
    fe_values[velocities].get_function_gradients(fem_solution, grad_u_values);
    fe_values[velocities].get_function_values   (fem_solution, u_values);
    fe_values[pressure].get_function_values     (fem_solution, p_values);
    q_points = fe_values.get_quadrature_points();
    exact_solution.get_gradient_list (q_points, grad_u_true_values);
    exact_solution.get_u_list (q_points, u_true_values);
    exact_solution.get_p_list (q_points, p_true_values);
    for(unsigned int q = 0; q < n_q_points; ++q)
    {
      grad_u_values[q] -= grad_u_true_values[q];
      u_values[q]      -= u_true_values[q];
      grad_error       += grad_u_values[q].norm_square() * fe_values.JxW(q);
      u_error          += u_values[q].norm_square()      * fe_values.JxW(q);
      //p_error          += pow(p_values[q] - p_mean*(1-fix_pressure) - p_true_values[q],2) 
                          //* fe_values.JxW(q);
      div_u            += pow(trace(grad_u_values[q]),2) * fe_values.JxW(q);
      p_error          += pow( (p_values[q] - 1*p_mean          )
                           -   (p_true_values[q] - 1*p_true_mean),2)
                                                         * fe_values.JxW(q);
    }
  }//end_for_cell
  grad_error = sqrt(grad_error);
  u_error    = sqrt(u_error);
  p_error    = sqrt(p_error);
  div_u      = sqrt(div_u);
  std::cout << "H10_norm = " << grad_error << std::endl;
  std::cout << "UL2_norm = " << u_error << std::endl;
  std::cout << "PL2_norm = " << p_error << std::endl;
  std::cout << "DIV_norm = " << div_u << std::endl;
  std::cout << "P_avg.   = " << p_mean << std::endl;

  EV[gul2].push_back(grad_error);
  EV[ul2].push_back(u_error);
  EV[pl2].push_back(p_error);
  EV[udiv].push_back(div_u);

}
//-------------------------------------------------------------------------
template<int dim, class FF, typename V>
void ConvergenceRates<dim, FF, V>::update_norms()
{
  double val;
  for(unsigned int f = 0; f < fields.size(); ++f)
  {
    Norms norm_type = norm_map[fields[f]];
    for(unsigned int k = 0; k < EV[0].size(); ++k)
    switch(norm_type)
    {
      case gul2:
        NV[f].push_back(EV[gul2][k]);
        break;
      case ul2:
        NV[f].push_back(EV[ul2][k]);
        break;
      case pl2:
        NV[f].push_back(EV[pl2][k]);
        break;
      case udiv:
        NV[f].push_back(EV[udiv][k]);
        break;
      case uhdiv:
        val = sqrt( EV[ul2][k]*EV[ul2][k] + EV[udiv][k]*EV[udiv][k]);
        NV[f].push_back(val);
        break;
      case uh1:
        val = sqrt( EV[ul2][k]*EV[ul2][k] + EV[gul2][k]*EV[gul2][k]);
        NV[f].push_back(val);
        break;
    }//end_switch
  }//end_fields

}
//-------------------------------------------------------------------
template class ConvergenceRates<2, StokesExactSolution2<2>, TrilinosWrappers::MPI::BlockVector >;
template class ConvergenceRates<2, DarcyExactSolution2<2>,  TrilinosWrappers::MPI::BlockVector >;
template class ConvergenceRates<3, StokesExactSolution2<3>, TrilinosWrappers::MPI::BlockVector >;
template class ConvergenceRates<3, DarcyExactSolution2<3>,  TrilinosWrappers::MPI::BlockVector >;
//-------------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
//------------------------------------------------------------
