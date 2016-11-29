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
ParallelMapLinker<dim, spacedim>::ParallelMapLinker()
{}
//---------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::reinit(
    const DoFHandler<spacedim> &stokes,
    const DoFHandler<spacedim> &darcy,
    const FiniteElement<dim,spacedim> &fe,
    const unsigned int         &id,
    const unsigned int         &refinements)
{
  stokes_dof_handler = &stokes;
  darcy_dof_handler  = &darcy;
  target_fe          = &fe;
  interface_id       = id;
  n_refinements      = refinements;

  find_cells_on_the_interface();
  initialize_source_dof_handler();
  split_communication();
  create_local_mesh();
  setup_dofs();
  build_source_target_map();
  verify_domains_match();
  relate_foreign_dofs();
  test_communication();
}
//---------------------------
template<int dim, int spacedim>
bool ParallelMapLinker<dim, spacedim>::Comparator::
        operator()(const Point<spacedim> &p1, const Point<spacedim> &p2) const
{
  for(unsigned int i = 0; i < spacedim; ++i)
  {
    //double diff = p2[i] - p1[i];
    //if(std::abs(diff) <  1e-9)
    if(p1[i] == p2[i])
      continue;
    else if(p1[i] < p2[i])
      return true;
    else //(p1[i] > p2[i])
      return false;
  }

  return false;
}
//---------------------------
//---------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::split_communication()
{
  int size,rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &size);

  worker_id = rank;
  n_workers = size;

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

  if(i_worker_id < n_sd_workers)
    MPI_Intercomm_create(sd_comm, 0, interface_comm, n_sd_workers, 111, &intercomm);
  else
    MPI_Intercomm_create(sd_comm, 0, interface_comm, 0, 111, &intercomm);

  int flag;
  MPI_Comm_test_inter(intercomm, &flag);
  if (flag == false)
    Assert (false, ExcInternalError());

  //Get foreign ids
  foreign_ids.resize (n_foreign_workers);
  MPI_Allgather(&worker_id, 1, MPI_UNSIGNED, 
                &foreign_ids[0], 1, MPI_UNSIGNED, intercomm);

}
//---------------------------
//---------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::find_cells_on_the_interface()
{
  owns_cells_on_the_interface = false;

  unsigned int faces_per_cell = GeometryInfo<spacedim>::faces_per_cell;
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

  if(domain_type == Stokes)
    std::cout << "Worker " << worker_id << " has domain type Stokes" << std::endl;
  if(domain_type == Darcy)
    std::cout << "Worker " << worker_id << " has domain type Darcy" << std::endl;

  target_dof_handler.initialize(triangulation, *target_fe);
  //target_dof_handler.distribute_dofs(*target_fe);
  solution_vector.reinit (target_dof_handler.n_dofs());

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

  unsigned int faces_per_cell = GeometryInfo<spacedim>::faces_per_cell;
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

  triangulation.create_triangulation (point_vector, cell_vector, SubCellData());
  triangulation.refine_global( n_refinements );
  double h_param = GridTools::maximal_cell_diameter(triangulation);
  if(worker_id == 0)
  {
    std::cout << "*******Embedded triangulation information********" << std::endl;
    std::cout << "                   h           : " 
              << h_param  << std::endl;
  }
}//end_create_local_mesh
//--------------------------------------------
//--------------------------------------------
template<>
double ParallelMapLinker<1, 2>::point_to_double (const Point<2> &p)
{
  double x = p[0];
  double y = p[1];
  return exp(x) + sin(y) + 0.5;
}
//--------------------------------------------
//--------------------------------------------
template<>
double ParallelMapLinker<2, 3>::point_to_double (const Point<3> &p)
{
  double x = p[0];
  double y = p[1];
  double z = p[2];
  return exp(x) + sin(y) + sqrt(z) - 0.5;
}
//--------------------------------------------
//--------------------------------------------
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
//the dof numbering, we send a hash correspondig to coordinates 
//of the support points. This idea only works when the finite 
//element space does indeed have the concept of support points. 
//(This holds for Lagrangian elements)
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::relate_foreign_dofs () 
{
  if(owns_cells_on_the_interface == false)
    return;

  unsigned int n_local_dofs = target_dof_handler.n_dofs();

  //Idea: Each worker compares his hash map with the full hash map of the other group

  //Exchange sizes
  int size = local_hash_vec.size();
  std::vector<int> in_size_vec (n_foreign_workers);
  MPI_Allgather(&size, 1, MPI_INT, &in_size_vec[0], 1, MPI_INT, intercomm);

  //for(unsigned int i = 0; i < in_size_vec.size(); ++i)
  //std::cout << " Worker: " << worker_id << " has " << in_size_vec[i] << std::endl;

  //Displacement vector
  std::vector<int> disp (n_foreign_workers, 0);
  for(unsigned int i = 1; i < n_foreign_workers; ++i)
    disp[i] = disp[i-1] + in_size_vec[i-1];

  //Incoming vector
  unsigned int total_size = std::accumulate(in_size_vec.begin(), in_size_vec.end(), 0);
  //std::cout << "The total size is: " << total_size << std::endl;
  std::vector<std::size_t> foreign_hash_vec (total_size);

  //Transmission
  MPI_Allgatherv(&local_hash_vec[0], local_hash_vec.size(), MPI_UNSIGNED_LONG, 
     &foreign_hash_vec[0], &in_size_vec[0], &disp[0], MPI_UNSIGNED_LONG, intercomm);

  //for(unsigned int i = 0; i < in_size_vec.size(); ++i)
  //std::cout << "Worker " << worker_id << " says his in_size_vec[" << i << "]: "
    //<< in_size_vec[i] << std::endl;

  //for(unsigned int i = 0; i <  foreign_hash_vec.size(); ++i)
    //std::cout << "Worker " << worker_id << " fvec[" << i << "]: " 
      //<< foreign_hash_vec[i] << std::endl;

  std::vector<bool> is_taken (n_local_dofs, false);
  send_size_indexed_by_foreign_worker.resize ( n_foreign_workers, 0 );
  recv_size_indexed_by_foreign_worker.resize ( n_foreign_workers, 0 );
  send_disp_indexed_by_foreign_worker.resize ( n_foreign_workers, 0 );
  recv_disp_indexed_by_foreign_worker.resize ( n_foreign_workers, 0 );
  recv_dof_index_vec.resize                  ( n_local_dofs );
  request_hash_vec.resize                    ( n_local_dofs );
  
  unsigned int hash_counter = 0;
  unsigned int dof_counter = 0;

  for(unsigned int i = 0; i < n_foreign_workers; ++i)
  {
    unsigned int n_elements = 0;

    for(unsigned int j = 0; j < in_size_vec[i]; ++j)
    {
      std::size_t hash = foreign_hash_vec[hash_counter++];
      auto it = local_hash_map.find(hash);
      if( it != local_hash_map.end())
      if(is_taken[it->second] == false)
      {
        is_taken[it->second] = true;
        request_hash_vec  [dof_counter] = it->first;
        recv_dof_index_vec[dof_counter] = it->second;
        ++dof_counter;
        ++n_elements;
        //std::cout << "I am worker " << worker_id 
        //<< " and I am updating " << n_elements << std::endl;
      }
    }

    recv_size_indexed_by_foreign_worker[i] = n_elements;

    if(n_local_dofs == dof_counter)
      break;

  }//end_for_each_foreign_worker

  //Total size
  n_elements_to_recv = 
  std::accumulate(recv_size_indexed_by_foreign_worker.begin(),
                  recv_size_indexed_by_foreign_worker.end(), 0);

  //Displacement vector for receive
  for(unsigned int i = 1; i < n_foreign_workers; ++i)
    recv_disp_indexed_by_foreign_worker[i] = 
        recv_disp_indexed_by_foreign_worker[i-1] 
      + recv_size_indexed_by_foreign_worker[i-1];


  for(unsigned int i = 0; i < n_foreign_workers; ++i)
  std::cout << "I am worker " << worker_id 
            << " and I expect to receive " 
            << recv_size_indexed_by_foreign_worker[i]
            << " units from worker " 
            << foreign_ids[i]
            << " with disp " 
            << recv_disp_indexed_by_foreign_worker[i]
            <<  std::endl;

  //Request information

  //Get send size
  MPI_Alltoall(&recv_size_indexed_by_foreign_worker[0], 1, MPI_INT, 
               &send_size_indexed_by_foreign_worker[0], 1, MPI_INT, 
               intercomm);

  //Total size
  n_elements_to_send = 
  std::accumulate(send_size_indexed_by_foreign_worker.begin(),
                  send_size_indexed_by_foreign_worker.end(), 0);

  std::vector<std::size_t> hash_data_to_send ( n_elements_to_send );

  //Displacement vector for send
  for(unsigned int i = 1; i < n_foreign_workers; ++i)
    send_disp_indexed_by_foreign_worker[i] = 
        send_disp_indexed_by_foreign_worker[i-1] 
      + send_size_indexed_by_foreign_worker[i-1];

  //Receive the request from each worker
  MPI_Alltoallv( &request_hash_vec[0], &recv_size_indexed_by_foreign_worker[0],
                 &recv_disp_indexed_by_foreign_worker[0], MPI_UNSIGNED_LONG,
                 &hash_data_to_send[0], &send_size_indexed_by_foreign_worker[0],
                 &send_disp_indexed_by_foreign_worker[0], MPI_UNSIGNED_LONG, intercomm);

  for(unsigned int i = 0; i < n_foreign_workers; ++i)
  std::cout << "I am worker " << worker_id 
            << " and I expect to send " 
            << send_size_indexed_by_foreign_worker[i]
            << " units to worker " 
            << foreign_ids[i]
            << " with disp " 
            << send_disp_indexed_by_foreign_worker[i]
            <<  std::endl;

  //Populate the send vector of indices 
  unsigned int counter = 0;
  send_dof_index_vec.resize( n_elements_to_send );
  for(const auto &hash : hash_data_to_send)
  {
    auto it = local_hash_map.find(hash);
    send_dof_index_vec[counter++] = it->second;
  }

  send_dof_value_vec.resize (n_elements_to_send, 0);
  recv_dof_value_vec.resize (n_elements_to_recv, 0);
}//relate_foreign_dofs 
//--------------------------------------------
//--------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::test_communication()
{
  std::vector<Point<spacedim> > support_points (target_dof_handler.n_dofs());
  DoFTools::map_dofs_to_support_points(MappingQ1<dim, spacedim>(), 
                                       target_dof_handler, support_points);
  std::vector<double> solution (target_dof_handler.n_dofs(), 0);
  std::transform(support_points.begin(), support_points.end(), 
                 solution.begin(),
                 [&](Point<spacedim> &p)->double 
                 { return this->point_to_double(p);});

  for(unsigned int i = 0; i < n_elements_to_send; ++i) 
    send_dof_value_vec[i] = solution[send_dof_index_vec[i]];

  std::vector<double> temp_recv_data (n_elements_to_recv, 0);

  MPI_Alltoallv(&send_dof_value_vec[0], &send_size_indexed_by_foreign_worker[0],
                &send_disp_indexed_by_foreign_worker[0], MPI_DOUBLE,
                &temp_recv_data[0], &recv_size_indexed_by_foreign_worker[0],
                &recv_disp_indexed_by_foreign_worker[0], MPI_DOUBLE,
                intercomm);

  for(unsigned int i = 0; i < n_elements_to_recv; ++i) 
    recv_dof_value_vec[i] = temp_recv_data[recv_dof_index_vec[i]];

  //for(unsigned int i = 0; i < n_elements_to_recv; ++i)
  //if(solution[i] != recv_dof_value_vec[i] )
  //{
    //std::cout 
      //<< "Worker" << worker_id << " says: "
      //<< "What I have: " 
      //<< solution[i] 
      //<< " what I need: " 
      //<< recv_dof_value_vec[i] << std::endl;
  //}
  
  if(solution == recv_dof_value_vec)
  {
    std::cout << "Работник " << worker_id << " сказал очень хорошо" << std::endl;
  }
                 
}
//--------------------------------------------
//--------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::verify_domains_match()
{
  if(owns_cells_on_the_interface == false)
    return;

  std::vector<Point<spacedim> > support_points ( target_dof_handler.n_dofs() );
  DoFTools::map_dofs_to_support_points(MappingQ1<dim,spacedim>(), 
                                       target_dof_handler, 
                                       support_points);
  //std::vector<std::size_t> hash_vec (support_points.size());
  local_hash_vec.resize ( support_points.size() );

  //Hash points
  std::transform(support_points.begin(), support_points.end(), local_hash_vec.begin(),
      [&](const Point<spacedim> &p)->std::size_t {return this->point_hasher(p);});

  //Populate local hash map
  unsigned int counter = 0;
  for(const auto &val : local_hash_vec)
    local_hash_map[val] = counter++;

  //if(local_hash_map.size() != target_dof_handler.n_dofs() )
    //Assert(false, ExcInternalError());
  

  //Exchange size within SD group
  int size = support_points.size();
  std::vector<int> in_size_vec (n_sd_workers);
  MPI_Gather(&size, 1, MPI_INT, &in_size_vec[0], 1, MPI_INT, 0, sd_comm);

  //Build disp
  std::vector<int> disp (n_sd_workers,0);
  for(unsigned int i = 1; i < disp.size(); ++i)
    disp[i] = disp[i-1] + in_size_vec[i-1];

  //Exchange local hash map within group
  unsigned int total_size = std::accumulate(in_size_vec.begin(), in_size_vec.end(), 0);
  std::vector<std::size_t> full_hash_vec (total_size);
  MPI_Gatherv(&local_hash_vec[0], local_hash_vec.size(), MPI_UNSIGNED_LONG, 
              &full_hash_vec[0], &in_size_vec[0], &disp[0], MPI_UNSIGNED_LONG, 0, sd_comm);


  if(sd_worker_id != 0)
    return;

  //Remove repetitions
  std::sort(full_hash_vec.begin(), full_hash_vec.end());
  auto last = std::unique(full_hash_vec.begin(), full_hash_vec.end());
  full_hash_vec.erase(last, full_hash_vec.end());

  size = full_hash_vec.size();

  int rec_size = -1;
  MPI_Status status;

  //Exchange size within leaders of each group
  MPI_Sendrecv(&size, 1, MPI_INT, 0, 111, 
               &rec_size, 1, MPI_INT, 0, 111, intercomm, &status);

  std::vector<std::size_t> foreign_hash_vec (rec_size);

  //Excahnge full hash within leaders of each group
  MPI_Sendrecv(&full_hash_vec[0], full_hash_vec.size(), MPI_UNSIGNED_LONG, 0, 111, 
               &foreign_hash_vec[0], rec_size, MPI_UNSIGNED_LONG, 0, 111, intercomm, &status);

  if(full_hash_vec != foreign_hash_vec)
  {
    std::cout << "Worker " << worker_id << " says that the vectors are not equal" << std::endl;
    Assert (false, ExcNotImplemented());
  }
  //for(unsigned int i = 0; i < full_hash_vec.size(); ++i)
    //std::cout << "Worker " << worker_id << " will received " << rec_size << std::endl;
}
//--------------------------------------------
//--------------------------------------------
template<int dim, int spacedim>
double ParallelMapLinker<dim,spacedim>::zero_to_pi(const double &num)
{
  //if(num == 1)
    //return numbers::PI;
  //else if (num == 0)
    //return numbers::E;
  return num;
}
//--------------------------------------------
//--------------------------------------------
template<int dim, int spacedim>
std::size_t ParallelMapLinker<dim,spacedim>::point_hasher(const Point<spacedim> &p)
{
  std::hash<double> hasher;
  std::size_t hash = hasher(zero_to_pi(p[0]));

  for(unsigned int i = 1; i <  spacedim; ++i)
    hash ^= hasher(zero_to_pi(p[i])) + 0x9e3779b9 + (hash << 6) + (hash >> 2);

  return hash;
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
