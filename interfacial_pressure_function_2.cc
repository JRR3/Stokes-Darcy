#include "interfacial_pressure_function_2.h"
//----------------------------------------------------------------
DEAL_II_NAMESPACE_OPEN
//----------------------------------------------------------------
template <int dim>
void
InterfacialPressureFunction2<dim>::value_list (const std::vector<Point<dim> > &points,
                                               std::vector<double> &values) const
{
  double x;
  double y;

  for(unsigned int k = 0; k < points.size(); ++k)
  {
    x = points[k][0];
    y = points[k][1];
    //--------------------------
    values[k] = -sin(y)*exp(1) + cos(x)*exp(y) + y*y - 2*y + 1;
  }
}
//----------------------------------------------------------------
template <int dim>
double
InterfacialPressureFunction2<dim>::value(const Point<dim> &point, const unsigned int ) const
{
  double x = point[0];
  double y = point[1];

  return (-sin(y)*exp(1) + cos(x)*exp(y) + y*y - 2*y + 1);
}

template class InterfacialPressureFunction2<2>;
//----------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
//----------------------------------------------------------------
