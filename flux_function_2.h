#ifndef FLUX_FUNCTION_2_H
#define FLUX_FUNCTION_2_H
//------------------------------------------------------------Dependencies
#include <deal.II/base/function.h>
#include <deal.II/lac/block_vector.h>
//------------------------------------------------------------
DEAL_II_NAMESPACE_OPEN
//------------------------------------------------------------
/** Flux exact solution from the paper:
 * Approximation of the Stokes–Darcy System by Optimization
 * See page 14.
 */
//------------------------------------------------------------
template<int dim>
class FluxFunction2 : public Function<dim>
{
  public:
  FluxFunction2() : Function<dim>(dim) {}

  virtual void vector_value (const Point<dim> &p,
                              Vector<double> &value) const;
  
};
//----------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
#endif
