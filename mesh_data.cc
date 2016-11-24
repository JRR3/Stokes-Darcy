//---------------------------
#include "mesh_data.h"
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/utilities.h>
//---------------------------
DEAL_II_NAMESPACE_OPEN
//---------------------------
template<int dim, int spacedim>
MeshData<dim, spacedim>::MeshData():
  worker_id (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
{}
//---------------------------
template<int dim, int spacedim>
void MeshData<dim, spacedim>::
compute_data(DoFHandler<dim,spacedim> &dof_handler)
{
  typename DoFHandler<dim,spacedim>::active_cell_iterator cell  = dof_handler.begin_active(),
                                                          ecell = dof_handler.end();
  QGauss<dim> q_rule (1);
  unsigned int n_q_points = q_rule.size();
  FEValues<dim, spacedim> fe_values (dof_handler.get_fe(), q_rule, update_JxW_values);
  volume = 0;
  max_h  = 0;
  min_h  = std::numeric_limits<double>::infinity();
  double diam;

  for(; cell != ecell; ++cell)
  if(cell->is_locally_owned())
  {
    fe_values.reinit(cell);

    diam = cell->diameter();

    if(diam < min_h)
      min_h = diam;

    if(max_h < diam)
      max_h = diam;

    for(unsigned int q = 0; q < n_q_points; ++q)
      volume += fe_values.JxW(q); 
  }

  min_h  = Utilities::MPI::min(min_h, MPI_COMM_WORLD);
  max_h  = Utilities::MPI::max(max_h, MPI_COMM_WORLD);
  volume = Utilities::MPI::sum(volume, MPI_COMM_WORLD);

}
//---------------------------
template class MeshData<3,3>;
template class MeshData<2,2>;
template class MeshData<2,3>;
//---------------------------
DEAL_II_NAMESPACE_CLOSE
//---------------------------
