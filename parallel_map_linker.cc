#include "parallel_map_linker.h"
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/grid_out.h>
//---------------------------
DEAL_II_NAMESPACE_OPEN
//---------------------------
template<int dim, int spacedim>
ParallelMapLinker<dim, spacedim>::ParallelMapLinker():
worker_id ( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ),
n_workers ( Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) ),
faces_per_cell ( GeometryInfo<spacedim>::faces_per_cell )
{
  //Object requires no specific initialization
}
//---------------------------
template<int dim, int spacedim>
bool
ParallelMapLinker<dim, spacedim>::Comparator::
        operator()(const Point<spacedim> &p1, const Point<spacedim> &p2) const
{
  for(unsigned int i = 0; i < spacedim; ++i)
  {
    double diff = p2[i] - p1[i];
    if(std::abs(diff) <  1e-9)
      continue;
    else if(diff > 0)
      return true;
    else //(p1[i] > p2[i])
      return false;
  }

  return false;
}
//---------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::
     attach_relevant_objects(  const DoFHandler<spacedim> &stokes,
                               const DoFHandler<spacedim> &darcy,
                               const unsigned int         &id,
                               const unsigned int         &refinements)
{
  stokes_dof_handler   = &stokes;
  stokes_dofs_per_cell = stokes_dof_handler->get_fe().dofs_per_cell;
  stokes_dofs_per_face = stokes_dof_handler->get_fe().dofs_per_face;
  //target_dof_handler   = &target;
  darcy_dof_handler   = &darcy;
  darcy_dofs_per_cell = darcy_dof_handler->get_fe().dofs_per_cell;
  darcy_dofs_per_face = darcy_dof_handler->get_fe().dofs_per_face;
  //target_dofs_per_cell = target_dof_handler->get_fe().dofs_per_cell;
  interface_id         = id;
  flux_refinements     = refinements;

  find_cells_on_the_interface();
  split_communication();
  create_meshes();
  //plot_triangulation();
}
//---------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::split_communication()
{
  int color = 0;

  if(owns_cells_on_the_interface == false)
    color = MPI_UNDEFINED;

  MPI_Comm_split(MPI_COMM_WORLD, color, worker_id, &interface_comm);

  if(owns_cells_on_darcy)
    color = 0;
  else
    color = MPI_UNDEFINED;

  MPI_Comm_split(interface_comm, color, worker_id, &darcy_comm);


  if(owns_cells_on_stokes)
    color = 0;
  else
    color = MPI_UNDEFINED;

  MPI_Comm_split(interface_comm, color, worker_id, &stokes_comm);

  n_interface_workers = 0;

  if(owns_cells_on_the_interface == true)
  {
    int size;
    MPI_Comm_size(interface_comm, &size);
    n_interface_workers = size;

    MPI_Comm_rank(interface_comm, &size);
    i_worker_id = size;
  }
  
  //if(owns_cells_on_darcy)
  //{
    //std::cout << "Worker " << worker_id << " is Darcy" << std::endl;
    //int darcy_rank;
    //MPI_Comm_size(darcy_comm, &darcy_rank);
    //std::cout << "# of Darcy members: " << darcy_rank << std::endl;
  //}
  //
  //if(owns_cells_on_stokes)
  //{
    //std::cout << "Worker " << worker_id << " is Stokes" << std::endl;
    //int stokes_rank;
    //MPI_Comm_size(stokes_comm, &stokes_rank);
    //std::cout << "# of stokes members: " << stokes_rank << std::endl;
  //}

  if(worker_id == 0)
  {
    std::cout << "# of Interface members: " 
      << Utilities::MPI::n_mpi_processes(interface_comm) << std::endl;
  }
}
//---------------------------
//---------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::find_cells_on_the_interface()
{
  owns_cells_on_the_interface = false;

  n_stokes_cells = 0;

  typename DoFHandler<spacedim>::active_cell_iterator 
    stokes_cell  = stokes_dof_handler->begin_active(),
    stokes_ecell = stokes_dof_handler->end();

  for(;stokes_cell != stokes_ecell; ++stokes_cell)
  if(stokes_cell->is_locally_owned())
  if(stokes_cell->at_boundary())
  for(unsigned int f = 0; f < faces_per_cell; ++f)
  if(stokes_cell->face(f)->boundary_id() == interface_id)
  {
    ++n_stokes_cells;
  }

  //----------------------------------------------------------
  //----------------------------------------------------------
  n_darcy_cells = 0;

  typename DoFHandler<spacedim>::active_cell_iterator 
    darcy_cell  = darcy_dof_handler->begin_active(),
    darcy_ecell = darcy_dof_handler->end();

  for(;darcy_cell != darcy_ecell; ++darcy_cell)
  if(darcy_cell->is_locally_owned())
  if(darcy_cell->at_boundary())
  for(unsigned int f = 0; f < faces_per_cell; ++f)
  if(darcy_cell->face(f)->boundary_id() == interface_id)
  {
    ++n_darcy_cells;
  }

  owns_cells_on_stokes        = (n_stokes_cells > 0);
  owns_cells_on_darcy         = (n_darcy_cells > 0);
  owns_cells_on_the_interface = owns_cells_on_stokes || owns_cells_on_darcy;
  owns_cells_on_both_domains  = owns_cells_on_stokes && owns_cells_on_darcy;
  

  std::cout << "Worker " << worker_id << " owns " << n_stokes_cells << " stokes cells " << std::endl;
  std::cout << "Worker " << worker_id << " owns " << n_darcy_cells << " darcy cells " << std::endl;

}
//---------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::create_meshes()
{
  if(owns_cells_on_stokes)
    create_local_mesh(stokes_dof_handler, 
                      stokes_lambda_triangulation,
                      stokes_flux_triangulation,
                      stokes_center_vector,
                      n_stokes_cells);
  if(owns_cells_on_darcy)
    create_local_mesh(darcy_dof_handler, 
                      darcy_lambda_triangulation,
                      darcy_flux_triangulation,
                      darcy_center_vector,
                      n_darcy_cells);
}
//--------------------------------------------------------
//--------------------------------------------------------
//--------------------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::create_local_mesh( 
    const DoFHandler<spacedim> * source_dof_handler, 
    Triangulation<dim, spacedim> &coarse_triangulation,
    Triangulation<dim, spacedim> &fine_triangulation,
    std::vector<double> &center_vector,
    const unsigned int           &n_cells)
{

  typename DoFHandler<spacedim>::active_cell_iterator 
    source_cell  = source_dof_handler->begin_active(),
    source_ecell = source_dof_handler->end();

  const unsigned int vertices_per_face = GeometryInfo<spacedim>::vertices_per_face;
  const bool is_new = true;

  std::map<Point<spacedim>, unsigned int, Comparator> map_point_to_num;

  std::vector<CellData<dim> > cell_vector (n_cells, CellData<dim>());

  unsigned int vertex_counter = 0;
  unsigned int cell_counter   = 0;
  unsigned int index;
  typedef std::map<Point<spacedim>, unsigned int, Comparator> MAP;

  for(;source_cell != source_ecell; ++source_cell)
  if(source_cell->is_locally_owned())
  if(source_cell->at_boundary())
  for(unsigned int f = 0; f < faces_per_cell; ++f)
  if(source_cell->face(f)->boundary_id() == interface_id)
  {
    for(unsigned int i = 0; i <  vertices_per_face; ++i)
    {
      //Add elements only if they do not exists
      Point<spacedim> p = source_cell->face(f)->vertex(i);
      std::pair<typename MAP::iterator, bool> status = 
            map_point_to_num.insert({p,vertex_counter});
      index = status.first->second;
      if(status.second == is_new)
        ++vertex_counter;
      cell_vector[cell_counter].vertices[i] = index;
    }

    //Extract the center from each face on the interface
    //Note that we pass the point as contiguous elements
    Point<spacedim> center = source_cell->face(f)->center();
    for(unsigned int i = 0; i < spacedim; ++i)
      center_vector.push_back(center[i]);

    //Specify the default material id
    cell_vector[cell_counter].material_id = 0;
    ++cell_counter;
  }

  std::vector<Point<spacedim> > point_vector (map_point_to_num.size());

  //Convert the keys to elements of a std::vector<Point<dim> >
  for( auto &it : map_point_to_num)
    point_vector[it.second] = it.first;

  //std::cout << "Worker " << worker_id << " says vertex counter: " << vertex_counter << std::endl;
  //std::cout << "Worker " << worker_id << " says cell   counter: " << cell_counter << std::endl;

  coarse_triangulation.create_triangulation (point_vector, cell_vector, SubCellData());
  fine_triangulation.copy_triangulation( coarse_triangulation );
  fine_triangulation.refine_global( flux_refinements );
}
//--------------------------------------------
//--------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::send_and_compare_center_vector
                      (CellType &cell_type,
                       std::vector<double> &center_vector)
                       
{
  //if(owns_cells_on_the_interface == false)
    //return;

  std::vector<int> vector_size (n_interface_workers);
  unsigned int size = center_vector.size();
  //std::cout << " worker: " << worker_id << " has " << size << std::endl;
  MPI_Allgather(&size, 1, MPI_UNSIGNED, &vector_size[0], 1, MPI_INT, interface_comm);
  
  std::vector<int> disp (n_interface_workers, 0);
  for(unsigned int i = 1; i < n_interface_workers; ++i)
    disp[i] = disp[i-1] + vector_size[i-1];

  unsigned int total_size = std::accumulate(vector_size.begin(), vector_size.end(), 0);
  std::vector<double> center_data (total_size);

  MPI_Allgatherv(&center_vector[0], vector_size[i_worker_id], MPI_DOUBLE, 
      &center_data[0], &vector_size[0], &disp[0], MPI_DOUBLE, interface_comm);

  std::vector<Point<spacedim> > 
  for(unsigned int i = 0; i < center_data.size(); i+=spacedim)


}
//--------------------------------------------
//--------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::plot_triangulation()
{
  if(owns_cells_on_the_interface == false)
    return;

  GridOut plotter;
  unsigned int out_index = 0;
  const std::string filename = ("grid-" +
                                Utilities::int_to_string (out_index, 2) +
                                "." +
                                Utilities::int_to_string
                                (worker_id, 2) +
                                ".vtu");
  std::ofstream output (filename.c_str());
  //plotter.write_vtu(lambda_triangulation, output);

  if(worker_id != 0)
    return;

  //std::vector<std::string> filenames;
  //for(unsigned int i = 0; i < n_workers; ++i)
  //filenames.push_back (std::string("grid-") +
                       //Utilities::int_to_string (out_index, 2) +
                       //"." +
                       //Utilities::int_to_string(i, 2) +
                       //".vtu");

}
//--------------------------------------------
//--------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::clear()
{
  stokes_dof_handler = NULL;
  darcy_dof_handler = NULL;
}
//--------------------------------------------------------------
//------------------------------------------------------------
template class ParallelMapLinker<1,2>;
template class ParallelMapLinker<2,3>;
//---------------------------
DEAL_II_NAMESPACE_CLOSE
//---------------------------
