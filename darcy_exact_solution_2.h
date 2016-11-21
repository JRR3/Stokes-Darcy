#ifndef DARCY_EXACT_SOLUTION_2_H
#define DARCY_EXACT_SOLUTION_2_H
//------------------------------------------------------------Dependencies
#include <deal.II/base/function.h>
#include <deal.II/lac/block_vector.h>
DEAL_II_NAMESPACE_OPEN
//------------------------------------------------------------
/** Darcy exact solution from the paper:
 * Approximation of the Stokes–Darcy System by Optimization
 * See page 14.
 */
//------------------------------------------------------------Darcy
template <int dim>
class DarcyExactSolution2 : public Function<dim>
{
public:
  DarcyExactSolution2 () : Function<dim>(dim+1) {}


  virtual void vector_value (const Point<dim> &p,
                             Vector<double>   &value) const;

  void get_div_list(const std::vector<Point<dim> > &p,
                       std::vector<double> &values) const;

  void get_gradient_list (const std::vector<Point<dim> > &points,
                                std::vector<Tensor<2,dim> > &gradients) const;

  void get_u_list (const std::vector<Point<dim> > &p,
                        std::vector<Tensor<1,dim> > &values) const;

  void get_rhs_list (const std::vector<Point<dim> > &p,
                        std::vector<Tensor<1,dim> > &values) const;

  void get_p_list (const std::vector<Point<dim> > &p,
                        std::vector<double> &values) const;
};
//}//end_namespace
DEAL_II_NAMESPACE_CLOSE
#endif
