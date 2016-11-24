#include "my_l2_norm.h"
#include <deal.II/fe/fe_values.h>
//---------------------------
DEAL_II_NAMESPACE_OPEN
//---------------------------
template<int dim, int spacedim, typename VEC>
L2_norm<dim, spacedim, VEC>::L2_norm(unsigned int degree):
  q_rule(degree),
  n_q_points( q_rule.size() )
{}
//---------------------------
template<int dim, int spacedim, typename VEC>
double L2_norm<dim, spacedim, VEC>::
compute_norm_sq(const DoFHandler<dim,spacedim> &dof_handler, const VEC &source, bool print) const
{
  FEValues<dim, spacedim> fe_values (dof_handler.get_fe(), q_rule, 
                                        update_values | 
                                        update_JxW_values);

  std::vector<double>                values (n_q_points);
  double                             norm_sq = 0;

  typename DoFHandler<dim, spacedim>::active_cell_iterator 
    cell  = dof_handler.begin_active(),
    ecell = dof_handler.end();

  for(; cell != ecell; ++cell)
  {
    fe_values.reinit(cell);
    fe_values.get_function_values(source, values);

    for(unsigned int q = 0; q < n_q_points; ++q)
    {
      norm_sq += values[q] * values[q] * fe_values.JxW(q);
    }//end_q_point

  }//end_for_each_cell

  if(print)
  std::cout << "The L2 norm squared of the given function   : "
            << norm_sq
            << std::endl;
  
  return norm_sq;
}
//---------------------------
template<int dim, int spacedim, typename VEC>
double L2_norm<dim, spacedim, VEC>::
compute_norm(const DoFHandler<dim,spacedim> &dof_handler, const VEC &source, bool print) const
{
  double norm = sqrt(compute_norm_sq(dof_handler, source, print));

  if(print)
  std::cout << "The L2 norm of the given function   : "
            << norm
            << std::endl;

  return norm;
}
//---------------------------
template class L2_norm <1,2>;
template class L2_norm <2,2,BlockVector<double> >;
template class L2_norm <2,3>;
template class L2_norm <2,3,BlockVector<double> >;
//---------------------------
DEAL_II_NAMESPACE_CLOSE
//---------------------------

