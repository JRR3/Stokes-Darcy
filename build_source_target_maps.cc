#include "parallel_map_linker.h"
#include <deal.II/lac/block_vector.h>
#include <deal.II/dofs/dof_tools.h>
//--------------------------------------------
//--------------------------------------------
DEAL_II_NAMESPACE_OPEN
//--------------------------------------------
//--------------------------------------------
template<int dim, int spacedim>
void ParallelMapLinker<dim, spacedim>::build_source_target_map()
{
  typename DoFHandler<spacedim>::active_cell_iterator 
    source_cell  = source_dof_handler->begin_active(),
    source_ecell = source_dof_handler->end();
  typename DoFHandler<dim, spacedim>::active_cell_iterator 
    target_cell   = target_dof_handler.begin_active(),
    target_ecell  = target_dof_handler.end();
  //------------------------------------------Data
  //Compute the maps (Cell,id)
  //------------------------------------------Stokes

  //------------------------------------------Flux
  unsigned int counter = 0;
  M_target_cell_id          target_cell_num;
  std::vector<Point<spacedim> >  target_center_vec;
  
  for(;target_cell != target_ecell; ++target_cell)
  {
    Point<spacedim> center = target_cell->center();
    target_center_vec.push_back(center);
    target_cell_num[target_cell] = counter;
    ++counter;

    //target_center_to_cell[center] = target_cell;
    //Used in the parallel processing
    //for(unsigned int i = 0; i < spacedim; ++i)
      //expanded_target_center_vec.push_back(center[i]);
  }
  //------------------------------------------
  //Here we build the maps (Source <-> Target)
  //-----------------------------------------------------------------Source
  for(typename M_source_cell_id::const_iterator source_it = source_cell_num.begin();
      source_it != source_cell_num.end(); ++source_it)
  {
    std::vector<DHIt> target_cell_vec;
    unsigned int source_id = source_it->second;
    source_cell  = source_it->first;

    for(typename M_target_cell_id::const_iterator target_it = target_cell_num.begin();
        target_it != target_cell_num.end(); ++target_it)
    {
      unsigned int target_id = target_it->second;
      Point<spacedim> target_center = target_center_vec[target_id];
      bool is_match = source_cell->point_inside(target_center);

      if(is_match)
      {
        unsigned int source_face= source_face_vec[source_id];
        P_cell_face source_pair = std::make_pair(source_it->first, source_face);
        target_cell_vec.push_back( target_it->first );
        target_to_source[target_it->first] = source_pair;
      }

    }//end_for_target_cell

    source_to_target[source_cell] = target_cell_vec;

  }//end_for_source_cell

}//end_build_map

//------------------------------------------------------------
template class ParallelMapLinker<1,2>;
template class ParallelMapLinker<2,3>;
//--------------------------------------------
DEAL_II_NAMESPACE_CLOSE
//--------------------------------------------
