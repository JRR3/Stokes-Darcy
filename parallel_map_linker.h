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
    Triangulation<dim,spacedim>       triangulation;
    MPI_Comm                          interface_comm;
    MPI_Comm                          sd_comm;
    MPI_Comm                          intercomm;
    const FiniteElement<dim,spacedim> * target_fe;
    const DoFHandler<spacedim>        * stokes_dof_handler;
    const DoFHandler<spacedim>        * darcy_dof_handler;
    const DoFHandler<spacedim>        * source_dof_handler;
    DoFHandler<dim,spacedim>          target_dof_handler;
    Vector<double>                    solution_vector;
    unsigned int                      interface_id;
    unsigned int                      n_cells;
    unsigned int                      n_elements_to_send;
    unsigned int                      n_elements_to_recv;
    unsigned int                      n_stokes_cells;
    unsigned int                      n_darcy_cells;
    unsigned int                      n_interface_workers;
    unsigned int                      i_worker_id;
    unsigned int                      sd_worker_id;
    unsigned int                      n_sd_workers;
    unsigned int                      n_foreign_workers;
    bool                              owns_cells_on_the_interface;
    bool                              owns_cells_on_both_domains;
    bool                              owns_cells_on_stokes;
    bool                              owns_cells_on_darcy;
    unsigned int                      worker_id;
    unsigned int                      n_workers;
    unsigned int                      n_refinements;
    struct Comparator
    {
      bool operator()(const Point<spacedim> &p1, const Point<spacedim> &p2) const;
    };
    typedef std::pair<Point<spacedim>,types::global_dof_index> P_point_dof;

  public:
    ParallelMapLinker();
    void reinit(
      const DoFHandler<spacedim> &stokes,
      const DoFHandler<spacedim> &darcy,
      const FiniteElement<dim,spacedim> &fe,
      const unsigned int         &id,
      const unsigned int         &refinements);
    void clear();

  private:
    std::vector<unsigned int>               source_face_vec;
    typedef typename DoFHandler<spacedim>::active_cell_iterator      DHIs;
    typedef typename DoFHandler<dim,spacedim>::active_cell_iterator  DHIt;
    typedef std::map<DHIs,unsigned int>     M_source_cell_id;
    typedef std::map<DHIt,unsigned int>     M_target_cell_id;
    typedef std::pair<DHIs, unsigned int >  P_cell_face;
    typedef std::map<Point<spacedim>, DHIt, Comparator> M_target_center_cell;
    typedef std::map<Point<spacedim>, types::global_dof_index, Comparator> M_point_dof;
    std::map<DHIs, std::vector<DHIt> >      source_to_target;
    std::map<DHIt, P_cell_face>             target_to_source;
    M_source_cell_id                        source_cell_num;
    std::vector<int>                        send_size_indexed_by_foreign_worker;
    std::vector<int>                        recv_size_indexed_by_foreign_worker;
    std::vector<int>                        send_disp_indexed_by_foreign_worker;
    std::vector<int>                        recv_disp_indexed_by_foreign_worker;
    std::vector<std::size_t>                request_hash_vec;
    std::vector<unsigned int>               recv_dof_index_vec;
    std::vector<unsigned int>               send_dof_index_vec;
    std::vector<double>                     send_dof_value_vec;
    std::vector<double>                     recv_dof_value_vec;
    std::vector<unsigned int>               foreign_ids;
    std::vector<std::size_t>                local_hash_vec;
    std::map<std::size_t, unsigned int>     local_hash_map;


  private:
    void find_cells_on_the_interface();
    void split_communication();
    void initialize_source_dof_handler();
    void create_local_mesh();
    void setup_dofs();
    void plot_triangulation();
    void build_source_target_map();
    double point_to_double(const Point<spacedim> &p);
    std::size_t point_hasher(const Point<spacedim> &p);
    double zero_to_pi(const double &num);
    void verify_domains_match();
    void relate_foreign_dofs();
    void test_communication();


};
//---------------------------
DEAL_II_NAMESPACE_CLOSE
//---------------------------
#endif
