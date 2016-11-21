#ifndef DARCY_K_INVERSE
#define DARCY_K_INVERSE
//------------------------------------------------------------Dependencies
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
DEAL_II_NAMESPACE_OPEN
//------------------------------------------------------------
/** Darcy inverse permeability tensor.
 */
//------------------------------------------------------------Darcy
template <int dim>
class KInverse : public TensorFunction<2,dim>
{
public:
  KInverse ()
    :
    TensorFunction<2,dim> ()
  {}

  virtual void value_list (const std::vector<Point<dim> > &points,
                           std::vector<Tensor<2,dim> >    &values) const;
};
//-------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
//-------------------------------------------------------------
#endif
