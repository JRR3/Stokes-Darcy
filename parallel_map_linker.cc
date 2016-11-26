#include "parallel_map_linker.h"
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <deal.II/lac/block_vector.h>
#include <deal.II/fe/mapping_q1.h>
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
faces_per_cell ( GeometryInfo<spacedim>::faces_per_cell ),
//-----------------------------------------------------
lambda_fe (2),
lambda_dof_handler ( lambda_triangulation),
flux_fe (0),
flux_dof_handler ( flux_triangulation)
{
  //Object requires no specific initialization
}
//---------------------------
template<int dim, int spacedim>
bool ParallelMapLinker<dim, spacedim>::Comparator::
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
  initialize_source_dof_handler();
  create_local_mesh();
  setup_dofs();
  build_maps();
  relate_foreign_dofs();
  //plot_triangulation();
}
//---------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::split_communication()
{
  int color;

  if(owns_cells_on_the_interface == false)
    color = MPI_UNDEFINED;
  else
    color = 0;

  MPI_Comm_split(MPI_COMM_WORLD, color, worker_id, &interface_comm);

  if(owns_cells_on_the_interface == false)
    color = MPI_UNDEFINED;
  else if(owns_cells_on_darcy)
    color = 0;
  else if(owns_cells_on_stokes)
    color = 1;

  MPI_Comm_split(interface_comm, color, worker_id, &sd_comm);


  if(owns_cells_on_the_interface == false)
    return;

  int size,rank;
  MPI_Comm_size(interface_comm, &size);
  n_interface_workers = size;
  MPI_Comm_rank(interface_comm, &rank);
  i_worker_id = rank;

  MPI_Comm_size(sd_comm, &size);
  n_sd_workers = size;
  MPI_Comm_rank(sd_comm, &rank);
  sd_worker_id = rank;

  //Number of workers in the other group
  n_foreign_workers = n_interface_workers - n_sd_workers;

  //unsigned int min_rank;
  //unsigned int id = i_worker_id;
  //MPI_Allreduce(&id, &min_rank, 1, MPI_UNSIGNED, MPI_MIN, sd_comm);

  if(i_worker_id < n_sd_workers)
    MPI_Intercomm_create(sd_comm, 0, interface_comm, n_sd_workers, 111, &intercomm);
  else
    MPI_Intercomm_create(sd_comm, 0, interface_comm, 0, 111, &intercomm);

  int flag;
  MPI_Comm_test_inter(intercomm, &flag);
  if (flag == false)
    Assert (false, ExcInternalError());

  //if(i_worker_id == 0)
    //std::cout << " Worker " << worker_id << " says: Flag is " << flag << std::endl;

  int val=worker_id;
  std::vector<int> reci (2);
  MPI_Allgather(&val, 1, MPI_INT, &reci[0], 1, MPI_INT, intercomm);

  //for(unsigned int i = 0; i < reci.size(); ++i)
    //std::cout << "I am worker " << worker_id << " and reci["
      //<< i << "]: " << reci[i] << std::endl;
}
//---------------------------
//---------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::find_cells_on_the_interface()
{
  owns_cells_on_the_interface = false;

  unsigned int counter = 0;

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
    //Point<spacedim> center = stokes_cell->face(f)->center();
    //source_face_center_vec.push_back(center);
    source_face_vec.push_back(f);
    source_cell_num[stokes_cell] = counter;
    ++counter;

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
    //Point<spacedim> center = darcy_cell->face(f)->center();
    //source_face_center_vec.push_back(center);
    source_face_vec.push_back(f);
    source_cell_num[darcy_cell] = counter;
    ++counter;

    ++n_darcy_cells;
  }

  //Initialize boolean variables
  owns_cells_on_stokes        = (n_stokes_cells > 0);
  owns_cells_on_darcy         = (n_darcy_cells > 0);
  owns_cells_on_the_interface = owns_cells_on_stokes || owns_cells_on_darcy;
  owns_cells_on_both_domains  = owns_cells_on_stokes && owns_cells_on_darcy;
  

  //std::cout << "Worker " << worker_id << " owns " << n_stokes_cells << " stokes cells " << std::endl;
  //std::cout << "Worker " << worker_id << " owns " << n_darcy_cells << " darcy cells " << std::endl;

  if(owns_cells_on_both_domains)
  {
    std::cout << "Owns cells on both domains" << std::endl;
    Assert (false, ExcNotImplemented());
  }


}
//--------------------------------------------------------
//--------------------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::initialize_source_dof_handler()
{
  if(owns_cells_on_the_interface == false)
    return;

  if(owns_cells_on_stokes)
  {
    source_dof_handler = stokes_dof_handler;
    n_cells = n_stokes_cells;
    domain_type = Stokes;
  } 

  if(owns_cells_on_darcy)
  {
    source_dof_handler = darcy_dof_handler;
    n_cells = n_darcy_cells;
    domain_type = Darcy;
  } 
  
}
//--------------------------------------------------------
//--------------------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::setup_dofs()
{
  if(owns_cells_on_the_interface == false)
    return;

  lambda_dof_handler.distribute_dofs(lambda_fe);
  lambda_vector.reinit (lambda_dof_handler.n_dofs());

  flux_dof_handler.distribute_dofs(flux_fe);
  flux_vector.reinit (flux_dof_handler.n_dofs());
}
//--------------------------------------------------------
//--------------------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::create_local_mesh()
{

  if(owns_cells_on_the_interface == false)
    return;

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

  lambda_triangulation.create_triangulation (point_vector, cell_vector, SubCellData());
  double h_lambda = GridTools::maximal_cell_diameter(lambda_triangulation);
  if(worker_id == 0)
  {
    std::cout << "*******Lambda information********" << std::endl;
    std::cout << "                   h           : " 
              << h_lambda  << std::endl;
  }
  flux_triangulation.copy_triangulation( lambda_triangulation );
  flux_triangulation.refine_global( flux_refinements );
}//end_create_local_mesh
//--------------------------------------------
//--------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::relate_foreign_dofs()
{ 
  relate_target_foreign_dofs (lambda_dof_handler);
}
//--------------------------------------------
//--------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::relate_target_foreign_dofs (
                       const DoFHandler<dim, spacedim>  &target_dof_handler)
{
  if(owns_cells_on_the_interface == false)
    return;

  M_point_dof  target_point_to_dof;
  DoFTools::map_support_points_to_dofs(MappingQ1<dim,spacedim>(),
                         target_dof_handler, target_point_to_dof);

  //if(worker_id == 0)
    //std::cout << "Size of xpand: " << expanded_dof_coord_vec.size() << std::endl;

  std::vector<double> expanded_dof_coord_vec (spacedim * target_point_to_dof.size());
  unsigned int count = 0;

  //Expand vector
  for(auto const &it : target_point_to_dof)
  for(unsigned int i = 0; i < spacedim; ++i)
    expanded_dof_coord_vec[count++] = it.first[i];

  std::vector<int> in_size_vec (n_foreign_workers);
  int my_vec_size = expanded_dof_coord_vec.size();
  MPI_Allgather(&my_vec_size, 1, MPI_INT, &in_size_vec[0], 1, MPI_INT, intercomm);

  //for(unsigned int i = 0; i < in_size_vec.size(); ++i)
  //std::cout << " Worker: " << worker_id << " has " << in_size_vec[i] << std::endl;

  //Displacement vector
  std::vector<int> disp (n_foreign_workers, 0);
  for(unsigned int i = 1; i < n_interface_workers; ++i)
    disp[i] = disp[i-1] + in_size_vec[i-1];

  //Incoming vector
  unsigned int total_size = std::accumulate(in_size_vec.begin(), in_size_vec.end(), 0);
  std::vector<double> incoming_data (total_size);

  //Transmission
  MPI_Allgatherv(&expanded_dof_coord_vec[0], expanded_dof_coord_vec.size(), MPI_DOUBLE, 
     &incoming_data[0], &in_size_vec[0], &disp[0], MPI_DOUBLE, intercomm);

  std::vector<bool> is_taken (target_dof_handler.n_dofs(), false);
  std::vector<std::vector<unsigned int> > data_received_ordered_by_worker (n_foreign_workers);

  unsigned int counter = 0;
  Point<spacedim> temp;

  //From the point of view of a processor on a certain domain,
  //he owns some cells on the boundary that have a corresponding
  //neighboring cell owned by a different processor belonging
  //to the other domain. Once we get the list of the degrees of 
  //freedom owned by each foreign processor, we determine from
  //which processor we expect to receive data. The list is ordered
  //by the rank of each processor on the foreign group. Afterwards,
  //we relate every dof that we own with a dof from the list of foreign
  //processors. Since the dof numbering in one processor might be different
  //to the numbering in the foreign processor, instead of sending 
  //the dof numbering, we send the coordinates of the support points.
  //This idea only works when the finite element space does indeed have
  //the concept of support points. (This holds for the Lagrangian elements)
  //The next step is to sort the support points so that both sender and
  //receiver agree on the order of the degrees of freedom.
  //
  std::vector<Point<spacedim> > support_points ( target_dof_handler.n_dofs() );
  DoFTools::map_dofs_to_support_points(MappingQ1<dim,spacedim>(), 
                                       target_dof_handler, 
                                       support_points);

  for(unsigned int i = 0; i < n_foreign_workers; ++i)
  {
    std::vector<unsigned int> temp_dof_vec;
    std::vector<P_point_dof> temp_pairs;

    for(unsigned int j = 0; j < in_size_vec[i]; ++j)
    {
      for(unsigned int k = 0; k < spacedim; ++k)
      {
        temp[k] = incoming_data[counter++]; 
      }

      auto it = target_point_to_dof.find(temp);
      if( it != target_point_to_dof.end())
      if(is_taken[it->second] == false)
      {
        P_point_dof pair (it->first, it->second);
        temp_pairs.push_back(pair);
        is_taken[it->second] = true;
        //std::cout << "Worker " << worker_id << " says center " 
          //<< temp << " is his." << std::endl;
      }
    }//end_elements_of_i-th_process

    std::vector<P_point_dof> clone_temp_pairs ( temp_pairs );

    //for(unsigned int p = 0; p < temp_pairs.size(); ++p)
    //{
      //if(worker_id == 0)
        //std::cout << "--- Worker " << worker_id << " owns p("  << p << "): " 
                  //<< "( " 
                  //<< temp_pairs[p].second << " , "
                  //<< temp_pairs[p].first 
                  //<< ") " << std::endl;
    //}

    //We sort the elements of the temp_pairs vector by coordinate
    //using a lambda expression. Note that we pass by reference.
    std::sort(temp_pairs.begin(), temp_pairs.end(), 
          [&](const P_point_dof &left, const P_point_dof &right) 
          { return Comparator()(left.first, right.first);});

    data_received_ordered_by_worker[i].resize(temp_pairs.size());
    //The ordering becomes important in 3D.

    //for(unsigned int p = 0; p < temp_pairs.size(); ++p)
    //{
      //if(worker_id == 0)
      //if(temp_pairs[p].second != clone_temp_pairs[p].second)
        //std::cout << "*** Worker " << worker_id << " says they are not equal " << std::endl;

        //std::cout << "*** Worker " << worker_id << " owns p("  << p << "): " 
                  //<< "( " 
                  //<< temp_pairs[p].second << " , "
                  //<< temp_pairs[p].first 
                  //<< ") " << std::endl;
    //}

  }//end_for_each_foreign_worker



}//end_compare_target_centers
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
  darcy_dof_handler  = NULL;
  source_dof_handler = NULL;
}
//--------------------------------------------------------------
//------------------------------------------------------------
template class ParallelMapLinker<1,2>;
template class ParallelMapLinker<2,3>;
//---------------------------
DEAL_II_NAMESPACE_CLOSE
//---------------------------
