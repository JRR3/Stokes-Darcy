foreign_worker_owns_target_cells_vec.resize(n_foreign_workers);
Point<spacedim> temp;
unsigned int counter = 0;

for(unsigned int i = 0; i < n_foreign_workers; ++i)
{
  std::vector<DHIt> temp_cell_vec;
  for(unsigned int j = 0; j < in_size_vec[i]; ++j)
  {
    for(unsigned int k = 0; k < spacedim; ++k)
    {
      temp[k] = center_data[counter++]; 
    }
    auto it = target_center_to_cell.find(temp);
    if( it != target_center_to_cell.end())
    {
      temp_cell_vec.push_back(it->second);
      //std::cout << "Worker " << worker_id << " says center " 
        //<< temp << " is his." << std::endl;
    }
  }

  foreign_worker_owns_target_cells_vec[i] = temp_cell_vec;
}//end_for_each_foreign_worker
This originates from the test of necessity of reordering
The test is not necessary in 2d but is necessary in 3D.
//for(unsigned int p = 0; p < temp_pairs.size(); ++p)
//{
  //if(temp_pairs[p].second != clone_temp_pairs[p].second)
    //std::cout << "*** Worker " << worker_id << " says they are not equal " << std::endl;

    //std::cout << "*** Worker " << worker_id << " owns p("  << p << "): " 
              //<< "( " 
              //<< temp_pairs[p].second << " , "
              //<< temp_pairs[p].first 
              //<< ") " << std::endl;
//}

