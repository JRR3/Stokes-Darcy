#ifndef MY_L2_NORM_H
#define MY_L2_NORM_H
#include <deal.II/lac/block_vector.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/base/quadrature_lib.h>
//-----------------------------------------
DEAL_II_NAMESPACE_OPEN
//------------------------------------------------------------
/** User defined L2_norm for a generic vector VEC. The user
 * must specify a dof_handler object and a vector.
 */
//------------------------------------------------------------
template<int dim, int spacedim, typename VEC = Vector<double> >
class L2_norm
{
  private:
    const QGauss<dim> q_rule;
    const unsigned int n_q_points;
  public:
    L2_norm(unsigned int degree);
    double compute_norm_sq(const DoFHandler<dim,spacedim> &dof_handler, 
                           const VEC &function, 
                           bool print = false) const;
    double compute_norm(const DoFHandler<dim,spacedim> &dof_handler, 
                        const VEC &function,
                        bool print = false) const;
  //private:
    //QGauss<dim> q_rule;
    
    
};
DEAL_II_NAMESPACE_CLOSE
//-----------------------------------------
#endif
