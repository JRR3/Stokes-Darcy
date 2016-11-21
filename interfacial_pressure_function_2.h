#ifndef INTERFACIAL_PRESSURE_FUNCTION_2_H
#define INTERFACIAL_PRESSURE_FUNCTION_2_H
//------------------------------------------------------------Dependencies
#include <deal.II/base/function.h>
#include <deal.II/lac/block_vector.h>
//------------------------------------------------------------
/** Darcy pressure solution from the paper:
 * Approximation of the Stokes–Darcy System by Optimization
 * See page 14.
 * See the function:
 * void SD<dim>::interpolate_interfacial_pressure ()
 */
//------------------------------------------------------------
DEAL_II_NAMESPACE_OPEN
template<int dim>
class InterfacialPressureFunction2 : public Function<dim>
{
  public:
  InterfacialPressureFunction2() : Function<dim>(1) {}

  virtual void value_list (const std::vector<Point<dim> > &points,
                                 std::vector<double> &values) const;
  virtual double value    (const Point<dim> &point, const unsigned int component = 0) const;
  
};
//----------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
#endif
