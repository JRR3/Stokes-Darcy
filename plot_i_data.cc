#include "plot_i_data.h"
//------------------------------------
DEAL_II_NAMESPACE_OPEN
//------------------------------------
template<int dim, int spacedim, typename V>
PlotIData<dim,spacedim,V>::PlotIData():
  counter (0)
{
}
//------------------------------------
template<int dim, int spacedim, typename V>
void PlotIData<dim,spacedim,V>::reset_counter()
{
  counter = 0;
}
//------------------------------------
template<int dim, int spacedim, typename V>
void PlotIData<dim,spacedim,V>::update_counter()
{
  //counter_str = std::static_cast<std::ostringstream*>(&(std::ostringstream() << counter))->str();
  counter_str = std::to_string(counter);
  ++counter;
}
//------------------------------------
template<int dim, int spacedim, typename V>
void PlotIData<dim,spacedim,V>::attach_dof_handlers(const DoFHandler<dim,spacedim> &dh_1, 
                                                    const DoFHandler<dim,spacedim> &dh_2)
{
  dh1 = &dh_1;
  dh2 = &dh_2;
}
//------------------------------------
template<int dim, int spacedim, typename V>
void PlotIData<dim,spacedim,V>::write_vtk (const V &v1, const V &v2)
{
  update_counter();
  DataOut<dim,DoFHandler<dim,spacedim> > data_out;
  data_out.attach_dof_handler (*dh1);
  data_out.add_data_vector (v1, "flux", DataOut<dim,DoFHandler<dim,spacedim> >::type_dof_data);
  data_out.build_patches ();
  std::string fname = "flux_" + counter_str + ".vtk";
  std::ofstream output1 ( fname.c_str() );
  data_out.write_vtk (output1);
  //-----------------------------------
  data_out.clear();
  data_out.attach_dof_handler (*dh2);
  data_out.add_data_vector (v2, "lambda");
  data_out.build_patches ();
  fname = "lambda_" + counter_str + ".vtk";
  std::ofstream output2 ( fname.c_str() );
  data_out.write_vtk (output2);
}
//------------------------------------
template<int dim, int spacedim, typename V>
void PlotIData<dim,spacedim,V>::save_count(const unsigned int &val)
{
  std::ofstream out("iterations.dat", std::ios_base::app);
  out << val << std::endl;
}
//------------------------------------
template<int dim, int spacedim, typename V>
void PlotIData<dim,spacedim,V>::save_flux_norm(const double &val)
{
  flux_vec.push_back(val);
}
//------------------------------------
template<int dim, int spacedim, typename V>
void PlotIData<dim,spacedim,V>::save_res_norm(const double &val)
{
  res_vec.push_back(val);
}
//------------------------------------
template<int dim, int spacedim, typename V>
void PlotIData<dim, spacedim, V>::generate_data(const unsigned int &cycle)
{
  std::string fname;
  fname = "o_flux_" + std::to_string(cycle) + ".dat";
  std::ofstream o_flux (fname.c_str());
  fname = "o_res_" + std::to_string(cycle) + ".dat";
  std::ofstream o_res (fname.c_str());

  for(unsigned int k = 0; k < res_vec.size(); ++k)
  {
    o_res << k << "\t" << res_vec[k] << std::endl;
  }

  for(unsigned int k = 0; k < flux_vec.size(); ++k)
  {
    o_flux << k << "\t" << flux_vec[k] << std::endl;
  }
}
//------------------------------------
template class PlotIData<1,2,Vector<double> >;
//------------------------------------
DEAL_II_NAMESPACE_CLOSE
//------------------------------------
