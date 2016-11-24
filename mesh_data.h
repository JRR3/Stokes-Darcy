#ifndef MESH_DATA_H
#define MESH_DATA_H
//---------------------------
#include <deal.II/dofs/dof_handler.h>
//---------------------------
DEAL_II_NAMESPACE_OPEN
//---------------------------
template<int dim, int spacedim>
class MeshData
{
  public:
    double min_h;
    double max_h;
    double volume;

  private:
    const unsigned int worker_id;

  public: 
    MeshData();
    void compute_data(DoFHandler<dim,spacedim> &dof_handler);
    
};
//---------------------------
DEAL_II_NAMESPACE_CLOSE
//---------------------------
#endif
