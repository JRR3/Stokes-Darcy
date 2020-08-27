#include "parallel_map_linker.h"
//#include <deal.II/lac/block_vector.h>
//#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/geometry_info.h>
//--------------------------------------------
//--------------------------------------------
DEAL_II_NAMESPACE_OPEN
//--------------------------------------------
//--------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::set_solution_vector( double &val)
{
  solution_vector = val;
}
//--------------------------------------------
//--------------------------------------------
/**
* Evaluate the target function at the q_points living
* on the source mesh. We assume that the map between the
* source cell and the target cell is injective, i.e., each
* source cell has only one target cell.
*/
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::
get_function_values(const DHIs                 &cell, 
                    FEValues<dim,spacedim>     &target_fe_values,
                    const Vector<double>       &function_vector, 
                          std::vector<double>  &values)
{
  target_fe_values.reinit(source_to_target[cell][0]);
  target_fe_values.get_function_values(function_vector, values);
}//end_get_function_values
//--------------------------------------------
//--------------------------------------------
/**
 * Compute the flux on each individual cell.
*/
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::compute_flux(
                                       Vector<double> &flux_target,
                                const MPI_BlockVector &stokes_source, 
                                const MPI_BlockVector &darcy_source)
{
  if(owns_cells_on_the_interface == false)
    return;

  const MPI_BlockVector             * source_vector;

  if(domain_type == Stokes)
    source_vector = &stokes_source;
  else
    source_vector = &darcy_source;

  //if(worker_id == 0)
    //std::cout << "Size of vec: " << source_vector->size() << std::endl;

  typename DoFHandler<spacedim>::active_cell_iterator 
    source_cell     = source_dof_handler->begin_active(),
    end_source_cell = source_dof_handler->end();

  QGauss<dim> face_q_rule (3);
  unsigned int n_face_q_points = face_q_rule.size();

  FEValues<dim,spacedim> flux_fe_values (*target_fe, face_q_rule, 
                                        //update_values |
                                        update_quadrature_points |
                                        update_JxW_values);

  //std::vector<Point<spacedim> > q_points;
  std::vector<Point<dim> >      source_unit_face_q_points;
  std::vector<double>           weights;
  std::vector<Tensor<1,spacedim> > source_values;
  std::vector<Point<spacedim> > normal_vectors;
  const FEValuesExtractors::Vector velocity (0);
  std::vector<types::global_dof_index> local_dof_indices (1);
  double worker_flux = 0;

  for(;source_cell != end_source_cell; ++source_cell)
  if(source_cell->is_locally_owned())
  if(source_cell->at_boundary())
  for(unsigned int fn = 0; fn < GeometryInfo<spacedim>::faces_per_cell; ++fn)
  if(source_cell->face(fn)->boundary_id() == interface_id)
  {
    compute_quadrature_points_and_weights(source_cell, 
                                                   fn, 
                                       flux_fe_values,
                                             //q_points,
                            source_unit_face_q_points,
                                             weights);

    source_values.resize(source_unit_face_q_points.size());
    //normal_vectors.resize(q_points.size());

    //unsigned int n_subcells = q_points.size() / n_face_q_points;

    Quadrature<dim> cell_quadrature (source_unit_face_q_points);

    //std::cout << "Size of quad: " << cell_quadrature.size() << std::endl;

    FEFaceValues<spacedim> source_fe_face_values (
                                        source_dof_handler->get_fe(),
                                        cell_quadrature,
                                        update_values |
                                        update_normal_vectors);

    source_fe_face_values.reinit( source_cell, fn );

    source_fe_face_values[velocity].get_function_values(*source_vector, 
                                                        source_values);
    normal_vectors = source_fe_face_values.get_normal_vectors();
    
    unsigned int q = 0;

    for(auto const &target_cell : source_to_target[source_cell])
    {
      target_cell->get_dof_indices( local_dof_indices );
      
      double flux = 0;
      for(unsigned int j = 0; j < n_face_q_points; ++j)
      {
        flux += source_values[q] * normal_vectors[q] * weights[q];
        ++q;
      }//end_for_j

      worker_flux += flux;
      solution_vector(local_dof_indices[0]) = flux;

    }//end_for_target_cell

  }//end_for_source_cell

  //std::cout << "***Worker " 
            //<< worker_id 
            //<< " of type " 
            //<< domain_type_str
            //<< " has local flux: " 
            //<< worker_flux
            //<< std::endl;

  transfer_data();
  print_flux();
}//end_compute_flux
//------------------------------------------------------------
//------------------------------------------------------------
/** This function shares the data in the solution_vector with the 
 * unique neighbor processor.
 */
template<int dim, int spacedim>
void ParallelMapLinker<dim,spacedim>::transfer_data()
{
  for(unsigned int i = 0; i < n_elements_to_send; ++i) 
    send_dof_value_vec[i] = solution_vector(send_dof_index_vec[i]);

  std::vector<double> temp_recv_data (n_elements_to_recv, 0);

  MPI_Alltoallv(&send_dof_value_vec[0], &send_size_indexed_by_foreign_worker[0],
                &send_disp_indexed_by_foreign_worker[0], MPI_DOUBLE,
                &temp_recv_data[0], &recv_size_indexed_by_foreign_worker[0],
                &recv_disp_indexed_by_foreign_worker[0], MPI_DOUBLE,
                intercomm);

  for(unsigned int i = 0; i < n_elements_to_recv; ++i) 
    recv_dof_value_vec[recv_dof_index_vec[i]] = temp_recv_data[i];

  for(unsigned int i = 0; i < solution_vector.size(); ++i) 
    solution_vector(i) += recv_dof_value_vec[i];

}
//------------------------------------------------------------
//------------------------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim,spacedim>::print_flux()
{

  typename DoFHandler<dim,spacedim>::active_cell_iterator
    cell = target_dof_handler.begin_active(),
   ecell = target_dof_handler.end();

  std::vector<types::global_dof_index> local_dof_indices (1);

  double flux = 0;

  for(; cell != ecell; ++cell)
  {
    cell->get_dof_indices( local_dof_indices );
    flux += solution_vector( local_dof_indices [0] );
  }

  double total_flux = 0;
  MPI_Reduce( &flux, &total_flux, 1, MPI_DOUBLE, MPI_SUM, 0, sd_comm);
  //MPI_Gather( &flux, 1 , MPI_DOUBLE, &total_flux, 1, MPI_DOUBLE, 0, sd_comm);


  std::cout << "***Worker " 
            << worker_id 
            << " of type " 
            << domain_type_str
            << " has local flux: " 
            << flux
            << " and total flux: "
            << total_flux
            << std::endl;

}
//------------------------------------------------------------
//------------------------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim,spacedim>::
compute_quadrature_points_and_weights(
                              const DHIs                      &source_cell, 
                              const unsigned int              &source_face_no,
                              FEValues<dim,spacedim>          &target_fe_values,
                              //std::vector<Point<spacedim> >   &q_points,
                              std::vector<Point<dim> >        &source_unit_face_q_points,
                              std::vector<double>             &weights)
{
  unsigned int n_q_points = target_fe_values.get_quadrature().size();
  unsigned int n_subfaces = source_to_target[source_cell].size();
  unsigned int n_source_unit_face_q_points = n_subfaces * n_q_points;
  std::vector<Point<spacedim> > q_points ( n_source_unit_face_q_points );
  weights.resize  ( n_source_unit_face_q_points );
  source_unit_face_q_points.resize  ( n_source_unit_face_q_points );
  auto q_points_it = q_points.begin();
  auto weights_it  = weights.begin();
  //auto lamb = [&](const DHIs &cell, const Point<spacedim> &point) -> Point<spacedim>
  //{return this->source_map.transform_real_to_unit_cell(cell, point);};
  for(auto &target_cell : source_to_target[source_cell])
  {
    target_fe_values.reinit(target_cell);
    std::copy(target_fe_values.get_quadrature_points().begin(), 
              target_fe_values.get_quadrature_points().end(), q_points_it );
    std::copy(target_fe_values.get_JxW_values().begin(), 
              target_fe_values.get_JxW_values().end(), weights_it );
    q_points_it  += n_q_points;
    weights_it   += n_q_points;
  }
  unsigned int q = 0;
  for(auto &point: q_points)
  {
    point = source_map.transform_real_to_unit_cell(source_cell, point);
    source_unit_face_q_points[q++] = 
                project_unit_cell_point_to_face(point, source_face_no);
  }
}//end_compute_q_points
//------------------------------------------------------------
//------------------------------------------------------------
template<int dim, int spacedim>
Point<dim> ParallelMapLinker<dim, spacedim>::
project_unit_cell_point_to_face(const Point<spacedim> &p, const unsigned int &face_no)
{
  if(spacedim == 2)
  {
    if(GeometryInfo<spacedim>::unit_normal_direction[face_no] == 0)
      return Point<spacedim-1> (p(1));
    else
      return Point<spacedim-1> (p(0));
  }
  else if (spacedim == 3)
  {
    if(GeometryInfo<spacedim>::unit_normal_direction[face_no] == 0)
      return Point<spacedim-1> (p(1), p(2));
    if(GeometryInfo<spacedim>::unit_normal_direction[face_no] == 1)
      return Point<spacedim-1> (p(0), p(2));
    else
      return Point<spacedim-1> (p(0), p(1));
  }
}
//------------------------------------------------------------
//--------------------------------------------
//--------------------------------------------
/**
 * Copy the solution vector contained in this object to the given
 * vector.
*/
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::copy_solution_vector( Vector<double>  &vec)
{
  vec = solution_vector;
}//end_get_function_values
//------------------------------------------------------------
template class ParallelMapLinker<1,2>;
template class ParallelMapLinker<2,3>;
//--------------------------------------------
DEAL_II_NAMESPACE_CLOSE
//--------------------------------------------
