#ifndef PARALLEL_MAP_LINKER_H
#define PARALLEL_MAP_LINKER_H
//---------------------------
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/component_mask.h>
//---------------------------
DEAL_II_NAMESPACE_OPEN
//------------------------------------------------------------
/** This class connects a mesh in R^{p} with a mesh in R^{q}.
 * Both meshes live in the same ambient space R^{spacedim}.
 */
//------------------------------------------------------------
template<int dim, int spacedim = dim+1>
class ParallelMapLinker
{
  private:
    enum                              CellType {Stokes, Darcy};
    MPI_Comm                          interface_comm;
    MPI_Comm                          sd_comm;
    MPI_Comm                          intercomm;
    const DoFHandler<spacedim>        * stokes_dof_handler;
    const DoFHandler<spacedim>        * darcy_dof_handler;
    unsigned int                      interface_id;
    unsigned int                      stokes_dofs_per_cell;
    unsigned int                      stokes_dofs_per_face;
    unsigned int                      darcy_dofs_per_cell;
    unsigned int                      darcy_dofs_per_face;
    unsigned int                      n_stokes_cells;
    unsigned int                      n_darcy_cells;
    unsigned int                      n_interface_workers;
    unsigned int                      n_sd_workers;
    unsigned int                      flux_refinements;
    unsigned int                      i_worker_id;
    unsigned int                      sd_worker_id;
    bool                              owns_cells_on_the_interface;
    bool                              owns_cells_on_both_domains;
    bool                              owns_cells_on_stokes;
    bool                              owns_cells_on_darcy;
    const unsigned int                worker_id;
    const unsigned int                n_workers;
    const unsigned int                faces_per_cell;
    Triangulation<dim,spacedim>       stokes_lambda_triangulation;
    Triangulation<dim,spacedim>       stokes_flux_triangulation;
    Triangulation<dim,spacedim>       darcy_lambda_triangulation;
    Triangulation<dim,spacedim>       darcy_flux_triangulation;
    std::vector<double>               stokes_center_vector;
    std::vector<double>               darcy_center_vector;
    struct Comparator
    {
      bool operator()(const Point<spacedim> &p1, const Point<spacedim> &p2) const;
    };
    //std::map<Point<spacedim>, unsigned int, Comparator> map_point_to_num;

  public:
    ParallelMapLinker();
    void clear();
    void attach_relevant_objects(  const DoFHandler<spacedim> &stokes,
                                   const DoFHandler<spacedim> &darcy,
                                   const unsigned int         &id,
                                   const unsigned int         &refinements);
    //void attach_source_dof_handler(const DoFHandler<spacedim, spacedim> &source);
    //void attach_target_dof_handler(const DoFHandler<dim, spacedim> &target);
    //void set_boundary_indicator(const unsigned int &id);

  private:
    void find_cells_on_the_interface();
    void create_meshes();
    void create_local_mesh( 
          const DoFHandler<spacedim> * source_dof_handler, 
          Triangulation<dim, spacedim> &coarse_triangulation,
          Triangulation<dim, spacedim> &fine_triangulation,
          std::vector<double>          &center_vector,
          const unsigned int           &n_cells);
    //void generate_local_domain();
    void plot_triangulation();
    void split_communication();
    void send_and_compare_center_vector
                      (CellType &cell_type,
                       std::vector<double> &center_vector);


};
//---------------------------
DEAL_II_NAMESPACE_CLOSE
//---------------------------
#endif
