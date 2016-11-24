#include "darcy_K_inverse.h"
//----------------------------------------------------------------
DEAL_II_NAMESPACE_OPEN
//-----------------------------------------------------------------
template <int dim>
void
KInverse<dim>::value_list (const std::vector<Point<dim> > &points,
                           std::vector<Tensor<2,dim> >    &values) const
{
  //Assert (points.size() == values.size(),
          //ExcDimensionMismatch (points.size(), values.size()));

  for (unsigned int p=0; p<points.size(); ++p)
  {
    values[p].clear ();

    for (unsigned int d=0; d<dim; ++d)
      values[p][d][d] = 1;
  }
}
//----------------------------------------------------------------
template class KInverse<2>;
template class KInverse<3>;
//----------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
//----------------------------------------------------------------
