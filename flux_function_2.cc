#include "flux_function_2.h"
//----------------------------------------------------------------
DEAL_II_NAMESPACE_OPEN
//----------------------------------------------------------------
template <int dim>
void
FluxFunction2<dim>:: vector_value (const Point<dim> &p,
                             Vector<double> &value) const
{
  double x = p[0];
  double y = p[1];
  //--------------------------
  value[0] = -x*(sin(y)*exp(1) + 2*(y-1));
  value[1] = -cos(y)*exp(1) + (y-1)*(y-1);
}

template class FluxFunction2<2>;
//----------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
//----------------------------------------------------------------
