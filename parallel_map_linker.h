#ifndef PARALLEL_MAP_LINKER_H
#define PARALLEL_MAP_LINKER_H
//---------------------------
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
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
    enum                              DomainType {Darcy, Stokes};
    DomainType                        domain_type;
    MPI_Comm                          interface_comm;
    MPI_Comm                          sd_comm;
    MPI_Comm                          intercomm;
    const DoFHandler<spacedim>        * stokes_dof_handler;
    const DoFHandler<spacedim>        * darcy_dof_handler;
    const DoFHandler<spacedim>        * source_dof_handler;
    //const DoFHandler<spacedim>        * target_dof_handler;
    unsigned int                      interface_id;
    unsigned int                      dofs_per_cell;
    unsigned int                      dofs_per_face;
    unsigned int                      stokes_dofs_per_cell;
    unsigned int                      stokes_dofs_per_face;
    unsigned int                      darcy_dofs_per_cell;
    unsigned int                      darcy_dofs_per_face;
    unsigned int                      n_cells;
    unsigned int                      n_stokes_cells;
    unsigned int                      n_darcy_cells;
    unsigned int                      n_interface_workers;
    unsigned int                      flux_refinements;
    unsigned int                      i_worker_id;
    unsigned int                      sd_worker_id;
    unsigned int                      n_sd_workers;
    unsigned int                      n_foreign_workers;
    bool                              owns_cells_on_the_interface;
    bool                              owns_cells_on_both_domains;
    bool                              owns_cells_on_stokes;
    bool                              owns_cells_on_darcy;
    const unsigned int                worker_id;
    const unsigned int                n_workers;
    const unsigned int                faces_per_cell;
    Triangulation<dim,spacedim>       lambda_triangulation;
    Triangulation<dim,spacedim>       flux_triangulation;
    struct Comparator
    {
      bool operator()(const Point<spacedim> &p1, const Point<spacedim> &p2) const;
    };

  public:
    ParallelMapLinker();
    void clear();
    void attach_relevant_objects(  const DoFHandler<spacedim> &stokes,
                                   const DoFHandler<spacedim> &darcy,
                                   const unsigned int         &id,
                                   const unsigned int         &refinements);
  private:
    FE_Q<dim, spacedim>               lambda_fe;
    DoFHandler<dim,spacedim>          lambda_dof_handler;
    Vector<double>                    lambda_vector;

    FE_DGQ<dim, spacedim>             flux_fe;
    DoFHandler<dim,spacedim>          flux_dof_handler;
    Vector<double>                    flux_vector;

    std::vector<unsigned int>         source_face_vec;

    typedef typename DoFHandler<spacedim>::active_cell_iterator      DHIs;
    typedef typename DoFHandler<dim,spacedim>::active_cell_iterator  DHIt;
    typedef std::map<DHIs,unsigned int>     M_source_cell_id;
    typedef std::map<DHIt,unsigned int>     M_target_cell_id;
    typedef std::pair<DHIs, unsigned int >  P_cell_face;
    typedef std::map<Point<spacedim>, DHIt, Comparator> M_target_center_cell;
    typedef std::map<Point<spacedim>, types::global_dof_index, Comparator> M_point_dof;
    std::map<DHIs, std::vector<DHIt> >      source_to_lambda;
    std::map<DHIs, std::vector<DHIt> >      source_to_flux;
    std::map<DHIt, P_cell_face>             lambda_to_source;
    std::map<DHIt, P_cell_face>             flux_to_source;
    M_point_dof                             lambda_point_to_dof;
    M_point_dof                             flux_point_to_dof;
    M_source_cell_id                        source_cell_num;
    std::vector<double>                     expanded_lambda_dof_coord_vec;
    std::vector<double>                     expanded_flux_dof_coord_vec;
    //std::vector<std::vector<DHIt> >         foreign_worker_owns_lambda_cells_vec;
    //std::vector<std::vector<DHIt> >         foreign_worker_owns_flux_cells_vec;


  private:
    void find_cells_on_the_interface();
    void split_communication();
    void initialize_source_dof_handler();
    void create_local_mesh();
    void setup_dofs();
    void create_maps();
    void plot_triangulation();
    void build_source_target_map(
      const DoFHandler<dim, spacedim>         &target_dof_handler,
      std::map<DHIs, std::vector<DHIt> >      &source_to_target,
      std::map<DHIt, P_cell_face>             &target_to_source);
    void build_maps();
    void relate_target_foreign_dofs(const DoFHandler<dim, spacedim>  &target_dof_handler);
    void relate_foreign_dofs();


};
//---------------------------
DEAL_II_NAMESPACE_CLOSE
//---------------------------
#endif
