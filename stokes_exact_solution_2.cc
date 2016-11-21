#include "stokes_exact_solution_2.h"
////------------------------------------------------------------Stokes
DEAL_II_NAMESPACE_OPEN
//----------------------------------------------------------------
template <int dim>
void
StokesExactSolution2<dim>::get_rhs_list (const std::vector<Point<dim> > &p,
                                           std::vector<Tensor<1,dim> > &values) const
{
  //double t = this->get_time();
  double x;
  double y;
  Tensor<1,dim> rhs;
  //--------------------------
  for(unsigned int k = 0; k < p.size(); ++k)
  {
    x = p[k][0];
    y = p[k][1];
    rhs[0] = -x*x*x - 6*x*(y - 1)*(y - 1) - exp(y)*sin(x);
    rhs[1] = -3*x*x*(y-1) - cos(y)*exp(1) + cos(x)*exp(y) + 2*y - 2;
    values[k] = rhs;

  }
}
//----------------------------------------------------------------
template <int dim>
void
StokesExactSolution2<dim>::get_u_list (const std::vector<Point<dim> > &p,
                                         std::vector<Tensor<1,dim> > &values) const
{
  //double t = this->get_time();
  double x;
  double y;
  Tensor<1,dim> u;
  //--------------------------
  for(unsigned int k = 0; k < p.size(); ++k)
  {
    x = p[k][0];
    y = p[k][1];
    //--------------------------
    u[0] = x*x*x*(y-1)*(y-1);
    u[1] = -exp(1)*cos(y);
    values[k] = u;

  }
}
//----------------------------------------------------------------
template <int dim>
void
StokesExactSolution2<dim>::get_div_list (const std::vector<Point<dim> > &p,
                                         std::vector<double> &values) const
{
  //double t = this->get_time();
  double x;
  double y;
  //--------------------------
  for(unsigned int k = 0; k < p.size(); ++k)
  {
    x = p[k][0];
    y = p[k][1];
    //--------------------------
    values[k] = 3*x*x*(y - 1)*(y - 1) + exp(1)*sin(y);
  }
}
//----------------------------------------------------------------
template <int dim>
void
StokesExactSolution2<dim>::get_p_list (const std::vector<Point<dim> > &p,
                                         std::vector<double> &values) const
{
  //double t = this->get_time();
  double x;
  double y;
  //--------------------------
  for(unsigned int k = 0; k < p.size(); ++k)
  {
    x = p[k][0];
    y = p[k][1];
    //--------------------------
    values[k] = cos(x)*exp(y) + y*y - 2*y + 1;
  }
}
//----------------------------------------------------------------
template <int dim>
void
StokesExactSolution2<dim>::get_stat_rhs_list (const std::vector<Point<dim> > &p,
                                          std::vector<Tensor<1,dim> > &values) const
{
  //double t = this->get_time();
  double x;
  double y;
  Tensor<1,dim> rhs;
  //--------------------------
  for(unsigned int k = 0; k < p.size(); ++k)
  {
    x = p[k][0];
    y = p[k][1];
    rhs[0] = -2*x*x*x - 6*x*(y - 1)*(y - 1) - exp(y)*sin(x);
    rhs[1] = -cos(y)*exp(1) + cos(x)*exp(y) + 2*y - 2;
    values[k] = rhs;

  }
}
//----------------------------------------------------------------
template <int dim>
void StokesExactSolution2<dim>::vector_value (const Point<dim> &p,
                                   Vector<double>   &values) const
{
  //double t = this->get_time();
  double x = p[0];
  double y = p[1];
  //--------------------------
  values[0]  = x*x*x*(y-1)*(y-1);
  values[1]  = -exp(1)*cos(y);
  values[2]  = cos(x)*exp(y) + y*y + -2*y + 1;
}

//----------------------------------------------------------------
//-----------------------------------------------------------
template<int dim>
void StokesExactSolution2<dim>::get_gradient_list (const std::vector<Point<dim> > &points,
                               std::vector<Tensor<2,dim> > &gradients) const
{
  //Using jacobian notation
  //double t  = this->get_time();
  double x,y;
  Tensor<2,dim> grad;
  grad      = 0;
  for(unsigned int k = 0; k < points.size(); ++k)
  {
    x = points[k][0];
    y = points[k][1];
    //---------------------
    grad[0][0] = 3*x*x*(y - 1)*(y - 1);
    grad[0][1] = 2*x*x*x*(y - 1);
    grad[1][0] = 0;
    grad[1][1] = exp(1)*sin(y);
    //---------------------
    gradients[k] = grad;
  }
}
////------------------------------------------------------------Stokes
template class StokesExactSolution2<2>;
////------------------------------------------------------------Stokes
DEAL_II_NAMESPACE_CLOSE
//------------------------------------------------------------Stokes
