#ifndef PLOT_I_DATA_H
#define PLOT_I_DATA_H
#include <deal.II/base/function.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <iostream>
#include <fstream>
#include <sstream>
//------------------------------------
DEAL_II_NAMESPACE_OPEN
//------------------------------------
/** This class takes care of the plotting of the evolution
 * of the flux along the interface.
 */
//------------------------------------------------------------
template<int dim, int spacedim, typename V>
class PlotIData
{
  private:
    const DoFHandler<dim,spacedim> *dh1;
    const DoFHandler<dim,spacedim> *dh2;
    unsigned int                    counter;
    std::string                     counter_str;
    std::vector<double>             res_vec;
    std::vector<double>             flux_vec;
  public:
    //------------------------------------------------
    PlotIData();
    void attach_dof_handlers(const DoFHandler<dim,spacedim> &dh_1, 
                             const DoFHandler<dim,spacedim> &dh_2);
    void save_count(const unsigned int &val);
    void save_flux_norm(const double &val);
    void save_res_norm(const double &val);
    void generate_data(const unsigned int &cycle);
    void reset_counter();
    void update_counter();
    void write_vtk( const V &v1, const V &v2);
};
//------------------------------------
DEAL_II_NAMESPACE_CLOSE
//------------------------------------
#endif
