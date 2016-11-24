#include "darcy_exact_solution_2.h"
//----------------------------------------------------------------
//namespace SDSPACE
//{
//using namespace dealii;

DEAL_II_NAMESPACE_OPEN
//-----------------------------------------------------------------
template <int dim>
void
DarcyExactSolution2<dim>::get_u_list (const std::vector<Point<dim> > &p,
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
    u[0] = -x*(sin(y)*exp(1) + 2*(y-1));
    u[1] = -cos(y)*exp(1) + (y-1)*(y-1);
    values[k] = u;

  }
}
//-----------------------------------------------------------------
template <int dim>
void
DarcyExactSolution2<dim>::get_rhs_list (const std::vector<Point<dim> > &p,
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
    u[0] = -(exp(1)*sin(y) + 2*y - 2)*x - exp(y)*sin(x);
    u[1] = (y - 1)*(y - 1) - 2*cos(y)*exp(1) + cos(x)*exp(y) + 2*y - 2;
    values[k] = u;

  }
}
//----------------------------------------------------------------
template <int dim>
void
DarcyExactSolution2<dim>::get_p_list (const std::vector<Point<dim> > &p,
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
    values[k] = -sin(y)*exp(1) + cos(x)*exp(y) + y*y - 2*y + 1;
  }
}
//----------------------------------------------------------------
template<int dim>
void DarcyExactSolution2<dim>::get_gradient_list (const std::vector<Point<dim> > &points,
                               std::vector<Tensor<2,dim> > &gradients) const
{
  double x,y;
  Tensor<2,dim> grad;
  grad      = 0;
  for(unsigned int k = 0; k < points.size(); ++k)
  {
    x = points[k][0];
    y = points[k][1];
    //---------------------
    grad[0][0] = -(sin(y)*exp(1) + 2*(y-1));
    grad[0][1] = -x*(cos(y)*exp(1) + 2);
    grad[1][0] = 0;
    grad[1][1] = sin(y)*exp(1) + 2*(y-1);
    //---------------------
    gradients[k] = grad;
  }
}
//-------------------------------------------------------------------
//----------------------------------------------------------------
template <int dim>
void
DarcyExactSolution2<dim>::get_div_list (const std::vector<Point<dim> > &p,
                                          std::vector<double> &values) const
{
  //double t = this->get_time();
  //double x;
  //double y;
  //--------------------------
  for(unsigned int k = 0; k < p.size(); ++k)
  {
    //x = p[k][0];
    //y = p[k][1];
    values[k] = 0;
  }
}

//----------------------------------------------------------------
template <int dim>
void
DarcyExactSolution2<dim>::vector_value (const Point<dim> &p,
                                   Vector<double>   &values) const
{
  //double t = this->get_time();
  double x = p[0];
  double y = p[1];
  //--------------------------
  values[0] = -x*(sin(y)*exp(1) + 2*(y-1));
  values[1] = -cos(y)*exp(1) + (y-1)*(y-1);
  values[2] = -sin(y)*exp(1) + cos(x)*exp(y) + y*y - 2*y + 1;
}
//----------------------------------------------------------------

//----------------------------------------------------------------
template class DarcyExactSolution2<2>;
template class DarcyExactSolution2<3>;
//----------------------------------------------------------------

DEAL_II_NAMESPACE_CLOSE
//}//end_namespace
