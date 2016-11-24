#ifndef MAP_LINKER_H
#define MAP_LINKER_H
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
template<int dim, int dim_T = dim-1, int spacedim = dim>
class MapLinker
{
  private:
    const DoFHandler<dim, spacedim>   * source_dof_handler;
    const DoFHandler<dim_T, spacedim> * target_dof_handler;
    FEFaceValues<dim, spacedim>       * source_fe_face_values;
    const Quadrature<dim_T>           * source_face_q_rule;
    const MappingQ1<dim_T, spacedim>  target_map;
    const MappingQ1<dim, spacedim>    source_map;
    typedef typename DoFHandler<dim, spacedim>::active_cell_iterator    DHIs;
    typedef typename DoFHandler<dim_T, spacedim>::active_cell_iterator  DHIt;
    typedef std::map<DHIs,unsigned int>     M_source_cell_id;
    typedef std::map<DHIt,unsigned int>     M_target_cell_id;
    typedef std::pair<unsigned int, double> pair_t;
    typedef std::pair<DHIs, unsigned int >  P_cell_face;
    std::map<DHIs, std::vector<DHIt> >      source_to_target;
    std::map<DHIt, P_cell_face>             target_to_source;
    unsigned int                            interface_id;
    unsigned int                            source_dofs_per_cell;
    unsigned int                            source_dofs_per_face;
    unsigned int                            target_dofs_per_cell;
    std::vector<Point<dim_T> >              source_unit_face_q_points;
    unsigned int                            n_source_unit_face_q_points;
    unsigned int                            worker_id;

  public:
    MapLinker();
    unsigned int                            n_cells_on_the_boundary;
    void attach_worker_id(const unsigned int &id);
    void attach_source_dof_handler(const DoFHandler<dim, spacedim> &source);
    void attach_target_dof_handler(const DoFHandler<dim_T, spacedim> &target);
    Point<dim_T> project_unit_cell_point_to_face(const Point<dim> &p, const unsigned int &face_no);
    void set_boundary_indicator(const unsigned int &id);
    void is_ready();
    void clear();
    void build_maps();
    void get_function_values_pwc(const DHIs &cell, 
                                 const Vector<double> &function_vector, 
                                 const std::vector<Point<spacedim> > &q_points,
                                       std::vector<double> &values);
    void get_function_values(const DHIs &cell, 
                             const Vector<double> &function_vector, 
                             const std::vector<Point<spacedim> > &q_points,
                                   std::vector<double> &values);
    void get_function_values( const DHIs                 &cell, 
                              const Vector<double>       &function_vector, 
                              FEValues<dim_T,spacedim>   &target_fe_values,
                                    std::vector<double>  &values);
    void compute_quadrature_points(
                              const DHIs                      &source_cell, 
                              const unsigned int              &source_face_no,
                              FEValues<dim_T, spacedim> &target_fe_values);

    void compute_quadrature_points(FEValues<dim_T, spacedim> &target_fe_values);

    void flux_integration_on_cell_subface(
                              const DHIs                      &source_cell, 
                              const unsigned int              &source_face_no,
                              const FEValues<dim_T, spacedim> &target_fe_values,
                              const double                    &scaling,
                              const Vector<double>            &target_function,
                              Vector<double>                  &target_rhs);

    void create_source_fe_face_values_object();

    void get_source_basis_functions_values(
                                  const FEValues<dim_T, spacedim> &target_fe_values,
                                  const double                    &scaling,
                                  const DHIt                      &cell, 
                                  const Vector<double>            &target_function,
                                  std::vector<unsigned int>       &source_dof_indices,
                                  Vector<double>                  &target_rhs);
    const P_cell_face& get_source_pair(const DHIt &target) ;
    DHIt get_target(const DHIs &source) ;
    const P_cell_face& operator[](const DHIt &target) ;
    void test_mapping();



};
//---------------------------
DEAL_II_NAMESPACE_CLOSE
//---------------------------
#endif
