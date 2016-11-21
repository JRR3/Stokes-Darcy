#include "map_linker.h"
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
//---------------------------
DEAL_II_NAMESPACE_OPEN
//---------------------------
template<int dim, int dim_T, int spacedim>
MapLinker<dim, dim_T, spacedim>::MapLinker():
owns_cells_on_the_interface( false )
{
}
//---------------------------
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::attach_worker_id(const unsigned int &id)
{
  worker_id = id;
}
//---------------------------
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::attach_source_dof_handler(const DoFHandler<dim, spacedim> &source)
{
  source_dof_handler   = &source;
  source_dofs_per_cell = source_dof_handler->get_fe().dofs_per_cell;
  source_dofs_per_face = source_dof_handler->get_fe().dofs_per_face;
}
//---------------------------
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::
        determine_if_we_own_cells_on_the_interface()
{
  typename DoFHandler<dim, spacedim>::active_cell_iterator 
    source_cell  = source_dof_handler->begin_active(),
    source_ecell = source_dof_handler->end();

  const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;

  for(;source_cell != source_ecell; ++source_cell)
  if(source_cell->is_locally_owned())
  if(source_cell->at_boundary())
  for(unsigned int f = 0; f < faces_per_cell; ++f)
  if(source_cell->face(f)->boundary_id() == interface_id)
  {
    owns_cells_on_the_interface = true;
  }
}
//---------------------------
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::generate_local_domain()
{
  if( owns_cells_on_the_interface == false )
    return;

  typename DoFHandler<dim, spacedim>::active_cell_iterator 
    source_cell  = source_dof_handler->begin_active(),
    source_ecell = source_dof_handler->end();

  std::ofstream out ("local_triangulation.inp");

  const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;

  for(;source_cell != source_ecell; ++source_cell)
  if(source_cell->is_locally_owned())
  if(source_cell->at_boundary())
  for(unsigned int f = 0; f < faces_per_cell; ++f)
  if(source_cell->face(f)->boundary_id() == interface_id)
  {
    return;
  }
}
//---------------------------
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::attach_target_dof_handler(const DoFHandler<dim_T, spacedim> &target)
{
  target_dof_handler = &target;
  target_dofs_per_cell = target_dof_handler->get_fe().dofs_per_cell;
}
//---------------------------
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::set_boundary_indicator(const unsigned int &id)
{
  interface_id = id;
}
//---------------------------
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::is_ready()
{
  if(target_dof_handler == NULL || source_dof_handler == NULL)
    throw std::invalid_argument("Some DoFHandler pointer is NULL in MapLinker");
}
//---------------------------
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::clear()
{
  target_dof_handler = NULL;
  source_dof_handler = NULL;
  source_to_target.clear();
  target_to_source.clear();
}
//---------------------------
template<int dim, int dim_T, int spacedim>
const typename MapLinker<dim, dim_T, spacedim>::P_cell_face& MapLinker<dim, dim_T, spacedim>::
get_source_pair(const DHIt &target) 
{
  return target_to_source[target];
}
//---------------------------
template<int dim, int dim_T, int spacedim>
typename MapLinker<dim, dim_T, spacedim>::DHIt MapLinker<dim, dim_T, spacedim>::
get_target(const DHIs &source) 
{
  std::vector<DHIt> target_vec = source_to_target[source];
  //std::cout << "Size of target_vec: " << target_vec.size() << std::endl;
  //std::cout << "ITC " << target_vec[0]->center() << std::endl;
  return target_vec[0];
}
//------------------------------------------------------------
/** This function projects a point inside the unit cell to the
 * face corresponding to face_no under the lexicographic order.
 * E.g. in 2D 
 * 0 -> left
 * 1 -> right
 * 2 -> bottom
 * 3 -> top
 */
//------------------------------------------------------------
template<int dim, int dim_T, int spacedim>
Point<dim_T> MapLinker<dim, dim_T, spacedim>::
project_unit_cell_point_to_face(const Point<dim> &p, const unsigned int &face_no)
{
  if(dim == 2)
  {
    if(GeometryInfo<dim>::unit_normal_direction[face_no] == 0)
      return Point<dim-1> (p(1));
    else
      return Point<dim-1> (p(0));
  }
  else if (dim == 3)
  {
    if(GeometryInfo<dim>::unit_normal_direction[face_no] == 0)
      return Point<dim-1> (p(1), p(2));
    if(GeometryInfo<dim>::unit_normal_direction[face_no] == 1)
      return Point<dim-1> (p(0), p(2));
    else
      return Point<dim-1> (p(0), p(1));
  }
}
//------------------------------------------------------------
/** This function computes quadrature points on a macrocell.
 * First it maps the quadrature points of all subcells to the unit macrocell
 * and then projects them to the corresponding face.
 */
//------------------------------------------------------------
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::
compute_quadrature_points( FEValues<dim_T, spacedim> &target_fe_values)
{
  DHIs source_cell = target_to_source.begin()->second.first;
  unsigned int source_face_no = target_to_source.begin()->second.second;
  compute_quadrature_points(source_cell, 
                            source_face_no,
                            target_fe_values);
}
//------------------------------------------------------------
/** This function computes quadrature points on a macrocell.
 * First it maps the quadrature points of all subcells to the unit macrocell
 * and then projects them to the corresponding face.
 */
//------------------------------------------------------------
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::
compute_quadrature_points(
                              const DHIs                      &source_cell, 
                              const unsigned int              &source_face_no,
                              FEValues<dim_T,spacedim>        &target_fe_values)
{
  unsigned int n_q_points = target_fe_values.get_quadrature().size();
  unsigned int n_subfaces = source_to_target[source_cell].size();
  n_source_unit_face_q_points = n_subfaces * n_q_points;
  std::vector<Point<spacedim> > q_points (n_subfaces * n_q_points);
  source_unit_face_q_points.resize(n_source_unit_face_q_points);
  auto it = q_points.begin();
  //auto lamb = [&](const DHIs &cell, const Point<spacedim> &point) -> Point<spacedim>
  //{return this->source_map.transform_real_to_unit_cell(cell, point);};
  for(auto &target_cell : source_to_target[source_cell])
  {
    target_fe_values.reinit(target_cell);
    std::copy(target_fe_values.get_quadrature_points().begin(), 
              target_fe_values.get_quadrature_points().end(), it);
    it += n_q_points;
  }
  unsigned int q = 0;
  for(auto &point: q_points)
  {
    point = source_map.transform_real_to_unit_cell(source_cell, point);
    source_unit_face_q_points[++q] = project_unit_cell_point_to_face(point, source_face_no);
  }
}
//------------------------------------------------------------
//------------------------------------------------------------
/** This function allocates a source_fe_face_values object using
 * dynamic memory. Careful with memory leaks.
 */
//------------------------------------------------------------
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::create_source_fe_face_values_object()
{
  source_face_q_rule    = new Quadrature<dim_T> (source_unit_face_q_points);
  source_fe_face_values = new 
      FEFaceValues<dim, spacedim> (source_dof_handler->get_fe(), 
                                   *source_face_q_rule, 
                                   update_values |
                                   update_normal_vectors);
}
//------------------------------------------------------------
//------------------------------------------------------------
/** This function tries to improve the function:
 * get_source_basis_functions_values 
 * by collecting all the quadrature points on the cell subfaces
 * and create only one FEFaceValues. In exchange, we have to loop
 * over the subfaces twice. An acceptable tradeoff.
 * */
//------------------------------------------------------------
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::
flux_integration_on_cell_subface(
                              const DHIs                      &source_cell, 
                              const unsigned int              &source_face_no,
                              const FEValues<dim_T, spacedim> &target_fe_values,
                              const double                    &scaling,
                              const Vector<double>            &target_function,
                              Vector<double>                  &target_rhs)

{

}//end_get_source_basis_functions
//------------------------------------------------------------
/** This function takes a target cell (a face F), a FEValues object,
 * a scaling factor (c) and a function supported on the target cell (f). 
 * The function computes the integral
 * \int_{F} c (phi_u \cdot \bfn) f \, dF
 * CHECK: Why target_fe_values is passed instead of having one stored in the class?
 * Because we alredy have one FEValues object created in the function calling
 * this function. See the functions: stokes_star and darcy_star.
 */
//------------------------------------------------------------
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::
get_source_basis_functions_values(
                              const FEValues<dim_T, spacedim> &target_fe_values,
                              const double                    &scaling,
                              const DHIt                      &cell, 
                              const Vector<double>            &target_function,
                              std::vector<unsigned int>       &source_dof_indices,
                              Vector<double>                  &target_rhs)

{
  //double scaling = sqrt(32.);
  P_cell_face source_pair = get_source_pair(cell);
  DHIs source_cell        = source_pair.first;
  unsigned int source_face= source_pair.second;
  source_cell->get_dof_indices(source_dof_indices);
 //std::cout << "----------------------------" << std::endl;
 //std::cout << "Source cell id: " << source_cell->id() << std::endl;

  //unsigned int source_dofs_per_face = source_dof_handler->get_fe().source_dofs_per_face;
  unsigned int n_q_points = target_fe_values.get_quadrature().size();

  std::vector<Point<spacedim> > q_points =
      target_fe_values.get_quadrature_points();
  std::vector<double> values (n_q_points);
      target_fe_values.get_function_values(target_function, values);
  //std::vector<types::global_dof_index> local_dof_indices (source_dofs_per_cell);
  std::vector<Point<dim_T> > face_q_points (n_q_points);

  const FEValuesExtractors::Vector velocities (0);
  //std::vector<unsigned int> mask = {1,4,7};

  //This function can be made more general.
  //We are getting the x component of the 
  for(unsigned int q = 0; q < q_points.size(); ++q)
  {
      Point<spacedim> point_in_unit_cell = 
        source_map.transform_real_to_unit_cell(source_cell, q_points[q]);

      face_q_points[q] = 
      project_unit_cell_point_to_face(point_in_unit_cell, source_face);
  }

  Quadrature<dim_T> q_rule (face_q_points);
  FEFaceValues<dim, spacedim> fe_face_values(source_dof_handler->get_fe(), 
                                   q_rule, 
                                   update_values |
                                   update_normal_vectors);

  fe_face_values.reinit(source_cell, source_face);

  double temp;

  for(unsigned int q = 0; q < n_q_points; ++q)
  for(unsigned int i = 0; i < source_dofs_per_cell; ++i)
  {
     temp = fe_face_values[velocities].value(i,q)           *
                      fe_face_values.normal_vector(q)       *
                      values[q]                             *
                      scaling                               *
                      target_fe_values.JxW(q);
     //std::cout << temp << std::endl;
     target_rhs(i) += temp;
     //if(fabs(temp) > 1e-10)
     //{
       //std::cout.precision(0);
       //std::cout << "Contribution from q_point " << q << " to dof at ";
       //std::cout.precision(4);
       //std::cout << map_dof_to_x[source_dof_indices[i]] << std::endl;
       //std::cout.precision(6);
       //std::cout << "Contribution: " << temp << std::endl;
       //std::cout << "++++++++++++++" << std::endl;
       //target_rhs(i) += temp;
     //}
  }//end_q_and_i

}//end_get_source_basis_functions
//---------------------------
//---------------------------------------
/**
  This is a simplified version of the function get_function_values.
  By using piecewise constant elements, we do not have to create
  q_rules, nor FEValues objects. 
  This function is used in the N'* operator.
  */
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::
get_function_values_pwc(const DHIs &cell, 
                        const Vector<double> &function_vector, 
                        const std::vector<Point<spacedim> > &q_points,
                              std::vector<double> &values)
{
  std::vector<types::global_dof_index> local_dof_indices 
                                       (target_dofs_per_cell);
  std::vector<DHIt> cell_vec = source_to_target[cell];
  for(unsigned int q = 0; q < q_points.size(); ++q)
  for(typename std::vector<DHIt>::const_iterator it  = cell_vec.begin(); 
                                                  it != cell_vec.end(); ++it)
  if((*it)->point_inside(q_points[q]))
  {
    //std::cout << "----------------------------" << std::endl;
    //std::cout << "Flux_cell center: " << (*it)->center() << std::endl;
    //std::cout << "Qpoint          : " << q_points[q] << std::endl;
    //std::cout << "Distance        : " << (*it)->center().distance(q_points[q]) << std::endl;
    (*it)->get_dof_indices(local_dof_indices);
    values[q] = function_vector(local_dof_indices[0]);
    break;
  }
}
//--------------------------------------------
/*!
  Evaluate the target function at the q_points living
  on the source mesh. We assume that the map between the
  source cell and the target cell is injective, i.e., each
  source cell has only one target cell.
  */
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::
get_function_values(const DHIs                 &cell, 
                    const Vector<double>       &function_vector, 
                    FEValues<dim_T,spacedim>   &target_fe_values,
                          std::vector<double>  &values)
{
  target_fe_values.reinit(source_to_target[cell][0]);
  target_fe_values.get_function_values(function_vector, values);
}//end_get_function_values
//--------------------------------------------
/*!
  Evaluate the target function at the q_points living
  on the source mesh. Note that for each q_point we
  create a q_rule and a FEValues object. At this moment
  this is unavoidable since the q_points might belong
  to different target cells. This function is used in the
  N'* operator when we integrate along the interface.
  For Q0 elements in the target, see get_function_values_pwc.
  */
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::
get_function_values(const DHIs &cell, 
                    const Vector<double> &function_vector, 
                    const std::vector<Point<spacedim> > &q_points,
                          std::vector<double> &values)
{
  //std::vector<DHIt> cell_vec = source_to_target[cell];
  Point<dim_T> point_in_unit_cell;
  std::vector<double> target_value (1);

  for(unsigned int q = 0; q < q_points.size(); ++q)
  for(auto &target_cell : source_to_target[cell])
  {
    //DHIt target_cell = *it;
    if(target_cell->point_inside(q_points[q]))
    {
      point_in_unit_cell = target_map.transform_real_to_unit_cell(target_cell, q_points[q]);
      Quadrature<dim_T> q_rule (point_in_unit_cell);
      FEValues<dim_T, spacedim> fe_values(target_dof_handler->get_fe(), q_rule,
                                    update_values);
      fe_values.reinit(target_cell);
      fe_values.get_function_values(function_vector, target_value);
      values[q] = target_value[0];
      break;
    }
  }//end_for_q_points

}//end_get_function_values
//---------------------------
template<int dim, int dim_T, int spacedim>
const typename MapLinker<dim, dim_T, spacedim>::P_cell_face& MapLinker<dim, dim_T, spacedim>::
operator[](const DHIt &target)
{
  return target_to_source[target];
}
//--------------------------------------------------------------
/**
  Connect two different meshes using a std::map object.
  We assume that the coarse object is the source and the fine one is the target.
  The inverse map is also created. Note that when the meshes have a different 
  level of refinement, the map source-->target is a one-to-many map. 
  The inverse map is always injective.
  This function works in 3D.
  */
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::build_maps()
{

  is_ready();

  //std::cout << ">>>Connect domains through target" << std::endl;

  const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;

  typename DoFHandler<dim, spacedim>::active_cell_iterator 
    source_cell  = source_dof_handler->begin_active(),
    source_ecell = source_dof_handler->end();
  typename DoFHandler<dim_T, spacedim>::active_cell_iterator 
    target_cell   = target_dof_handler->begin_active(),
    target_ecell  = target_dof_handler->end();
  //------------------------------------------Data
  //Compute the maps (Cell,id)
  //------------------------------------------Stokes

  unsigned int counter = 0;
  M_source_cell_id          source_cell_num;
  std::vector<Point<spacedim> > source_center_vec;
  std::vector<unsigned int> source_face_vec;
  

  for(;source_cell != source_ecell; ++source_cell)
  if(source_cell->is_locally_owned())
  if(source_cell->at_boundary())
  for(unsigned int f = 0; f < faces_per_cell; ++f)
  if(source_cell->face(f)->boundary_id() == interface_id)
  {
    Point<spacedim> center = source_cell->face(f)->center();
    source_center_vec.push_back(center);
    source_face_vec.push_back(f);
    source_cell_num[source_cell] = counter;
    ++counter;
  }

  n_cells_on_the_boundary = counter;
  //------------------------------------------Flux
  counter = 0;
  M_target_cell_id          target_cell_num;
  std::vector<Point<spacedim> >  target_center_vec;
  
  for(;target_cell != target_ecell; ++target_cell)
  {
    Point<spacedim> center = target_cell->center();
    target_center_vec.push_back(center);
    target_cell_num[target_cell] = counter;
    ++counter;
  }
  //------------------------------------------
  //Here we build the maps (Source <-> Target)
  //-----------------------------------------------------------------Source
  for(typename M_source_cell_id::const_iterator source_it = source_cell_num.begin();
      source_it != source_cell_num.end(); ++source_it)
  {
    std::vector<DHIt> target_cell_vec;
    unsigned int source_id = source_it->second;
    source_cell  = source_it->first;
    //Point<spacedim> source_center = source_center_vec[source_id];

    for(typename M_target_cell_id::const_iterator target_it = target_cell_num.begin();
        target_it != target_cell_num.end(); ++target_it)
    {
      unsigned int target_id = target_it->second;
      Point<spacedim> target_center = target_center_vec[target_id];
      //bool is_match = target_center.distance(source_center) < target_cell_h;
      bool is_match = source_cell->point_inside(target_center);

      if(is_match)
      {
        unsigned int source_face= source_face_vec[source_id];
        P_cell_face source_pair = std::make_pair(source_it->first, source_face);
        target_cell_vec.push_back( target_it->first );
        //source_to_target[source_it->first] = target_it->first;
        target_to_source[target_it->first] = source_pair;
        //std::cout << "We have a match in cell " << source_it->first->id() << std::endl;
            //std::cout << "SC " 
                      //<< source_pair.first->face(source_pair.second)->center() 
                      //<< std::endl;
            //std::cout << "TC " << target_it->first->center() << std::endl;
      }

    }//end_for_target_cell

    //source_to_target[source_it->first] = target_cell_vec;
    source_to_target[source_cell] = target_cell_vec;

  }//end_for_source_cell

}//end_build_map
//--------------------------------------------------------------
/** Following the Test-Driven-Design paradigm, this test
 * verifies that the connection between the domains is correct.
 * We consider
 */
template<int dim, int dim_T, int spacedim>
void MapLinker<dim, dim_T, spacedim>::test_mapping()
{
  for(auto &source_map : source_to_target)
  {
    std::cout << " I am worker : " << worker_id << std::endl;
    std::cout << " I am inside source cell: " << source_map.first->id() << std::endl;
    std::cout << " The center of this cell: " << source_map.first->center() << std::endl;
    for(auto &target_cell    : source_map.second)
    {
      std::cout << ">>> Child cell center: " << target_cell->center() << std::endl;
    }
  }

}
//------------------------------------------------------------
template class MapLinker<2>;
//template class MapLinker<3>;
//---------------------------
DEAL_II_NAMESPACE_CLOSE
//---------------------------
