/**
* This version of the function compute_rho_flux
* uses the local averaging method which allows 
* one to use discontinuous elements of higher order.
* By default one should use DG0.
* 
* See the alternate projection version 
* in the compute_flux_projection.w file.
* 
* This function is used to compute the rho_function
* (a piecewise constant function) that arises from 
* the derivative operator N' and the operator N.
* 
* This function also returns the L2 norm of the flux.
* 
* Possible improvements:
* Remove dependency on the FEFieldFunction
*
* Should we encapsulate the effects of this function and 
* use the MapLinker class. That way we do not have to declare
* extraneous types inside the main class.
*/


template<int dim>
double SD<dim>::compute_rho_flux(   Vector<double> &flux_target,
                         const MPI_BlockVector &stokes_source, 
                         const MPI_BlockVector &darcy_source,
                         bool  print)
{
    TimerOutput::Scope timer_section(computing_timer, "Build flux function");

    pcout << "Compute rho flux" << std::endl;  

    double local_cell_flux;
    flux_target = 0;
    Vector<unsigned int> count_vector (flux_target.size());
    count_vector = 0;
    //std::ofstream out_0 ("rho_flux_q0.txt");
    //std::ofstream out_1 ("rho_flux_q1.txt");

    //std::vector<Point<dim> >      supp_points = flux_fe.get_unit_support_points();
    //Quadrature<dim> supp_q_rule ( supp_points );
    //QIterated<dim-1> face_q_rule (QTrapez<1>(), darcy_degree+2);

    //FEFaceValues<dim> stokes_fe_face_values (stokes_fe, face_q_rule, update_values
                                                        //| update_normal_vectors);
    //FEFaceValues<dim> darcy_fe_face_values (darcy_fe, face_q_rule, update_values
                                                        //| update_normal_vectors);
    FEValues<dim-1,dim> flux_fe_values (flux_fe, face_q_rule, 
                                          update_values |
                                          update_quadrature_points |
                                          update_JxW_values);

    //const unsigned int n_face_q_points    = face_q_rule.size();
    const unsigned int dofs_per_cell      = flux_fe.dofs_per_cell;
    //const unsigned int dofs_per_face      = flux_fe.dofs_per_face;
    //pcout << "Size of q rule: " << n_face_q_points << std::endl;
    //pcout << "Dofs per cell: " << dofs_per_cell  << std::endl;

    //Vector<double> local_rhs (dofs_per_cell);
    std::vector<double> flux_values (n_face_q_points);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    //Field functions
    Functions::FEFieldFunction<dim, DoFHandler<dim>, MPI_BlockVector > 
     darcy_field (darcy_dof_handler, darcy_source);
    Functions::FEFieldFunction<dim, DoFHandler<dim>, MPI_BlockVector > 
     stokes_field (stokes_dof_handler, stokes_source);


    //std::vector<Tensor<1,dim> >          darcy_face_values (n_face_q_points);
    //std::vector<Tensor<1,dim> >          stokes_face_values (n_face_q_points);
    std::vector<Vector<double> >          darcy_face_values (n_face_q_points);
    std::vector<Vector<double> >          stokes_face_values (n_face_q_points);
    std::vector<Point<dim> >              face_q_points (n_face_q_points);

    const FEValuesExtractors::Vector velocities (0);

    typename DoFHandler<dim>::active_cell_iterator stokes_cell, darcy_cell;
    typename DoFHandler<dim-1,dim>::active_cell_iterator 
      flux_cell  = flux_dof_handler.begin_active(),
      flux_ecell = flux_dof_handler.end();

    //unsigned int darcy_face, stokes_face;
    double c1, c1_JxW;
    double total_flux = 0;
    double l2_norm_sq = 0;

    for(; flux_cell != flux_ecell; ++flux_cell)
    {
       //Reinit flux
      flux_fe_values.reinit(flux_cell);
      local_cell_flux = 0;
      //--------------------Using only native elements of the class
      //The original version uses the derived types introduced in the 
      //MapLinker class.

      //CHECK: Modification in this line from using a mixed type variable to a native type
      darcy_cell  = darcy_flux_map.get_source_pair(flux_cell).first;
      stokes_cell = stokes_flux_map.get_source_pair(flux_cell).first;

      //Reinit Stokes and Darcy
      stokes_field.set_active_cell(stokes_cell);
      darcy_field.set_active_cell(darcy_cell);

      face_q_points = flux_fe_values.get_quadrature_points();
      stokes_field.vector_value_list(face_q_points, stokes_face_values);
      darcy_field.vector_value_list(face_q_points, darcy_face_values);

      for(unsigned int q = 0; q < n_face_q_points; ++q)
      {
          //Use when the stokes and darcy domains are not squares.
          //We assume that the spatial component (x_n) is the normal to the face.

          //c1 = 
          //( stokes_face_values[q] * stokes_fe_face_values.normal_vector(q) 
            //+
          //darcy_face_values[q] * darcy_fe_face_values.normal_vector(q) );
          //local_cell_flux += c1 * flux_fe_values.JxW(q);

          c1               = darcy_face_values[q](dim-1) - stokes_face_values[q](dim-1);
          c1_JxW           = c1 * flux_fe_values.JxW(q);
          local_cell_flux += c1_JxW;
          l2_norm_sq      += c1 * c1 * flux_fe_values.JxW(q);

      }//end_n_face_q_points


      //Update 
      flux_cell->get_dof_indices(local_dof_indices);
      for(unsigned int k = 0; k < dofs_per_cell; ++k)
      {
        flux_target(local_dof_indices[k])  += local_cell_flux;
        count_vector(local_dof_indices[k]) += 1;
      }


      //Print data
      if(print)
      for(unsigned int k = 0; k < dofs_per_cell; ++k)
      pcout << "Local flux(" 
                << local_dof_indices[k]
                << "): "
                << local_cell_flux
                << std::endl;

      total_flux += local_cell_flux;
     
    }//end_for_each_cell

    //Average
    for(unsigned int k = 0; k < flux_target.size(); ++k)
    {
      flux_target(k) /= count_vector(k);
    }

    if(print)
    pcout << "Total flux (sum) : " << total_flux << std::endl;
    
    return sqrt(l2_norm_sq);
    
}//end_rho_flux_fun
//---------------------------------------------------
/*!
  This function computes the pointwise flux function in a continuous
  space. Therefore it does not describe a local average 
  but a pointwise quantity.
  This function is used by the N'* operator, which returns a (LAMBDA)
  function. 
  Technically, if the meshes coincide and we are using
  a Lagrangian basis, one can simplify this step.
  We assume the target FEM space is Lagrangian. 
  Thus, intead of projecting we interpolate.
  */
//---------------------------------------------------
template<int dim>
void SD<dim>::compute_flux_function(  Vector<double> &lambda_target,
                           const MPI_BlockVector &stokes_source, 
                           const MPI_BlockVector &darcy_source)
{
    pcout << "Interpolate jump function (Q_" << lambda_degree << " space)" << std::endl;  

    Functions::FEFieldFunction<dim, DoFHandler<dim>, MPI_BlockVector > 
     darcy_field (darcy_dof_handler, darcy_source);
    Functions::FEFieldFunction<dim, DoFHandler<dim>, MPI_BlockVector > 
     stokes_field (stokes_dof_handler, stokes_source);
                                             //update_normal_vectors);

    //Q_rule using supp points.
    std::vector<Point<dim-1> > unit_cell_supp_points = lambda_fe.get_unit_support_points();
    Quadrature<dim-1> face_q_rule (unit_cell_supp_points);
    unsigned int n_face_q_points = face_q_rule.size();

    FEValues<dim-1,dim> lambda_fe_values (lambda_fe, face_q_rule, 
                                          //update_values | 
                                          //update_JxW_values |
                                          update_quadrature_points 
                                          );

    const unsigned int dofs_per_cell = lambda_fe.dofs_per_cell;


    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    //std::vector<Tensor<1,dim> >          darcy_face_values (n_face_q_points);
    //std::vector<Tensor<1,dim> >          stokes_face_values (n_face_q_points);
    std::vector<Vector<double> >         darcy_face_values (n_face_q_points);
    std::vector<Vector<double> >         stokes_face_values (n_face_q_points);
    std::vector<Point<dim> >             face_q_points (n_face_q_points);

    const FEValuesExtractors::Vector velocities (0);

    typename DoFHandler<dim>::active_cell_iterator stokes_cell, darcy_cell;
    typename DoFHandler<dim-1,dim>::active_cell_iterator 
      lambda_cell  = lambda_dof_handler.begin_active(),
      lambda_ecell = lambda_dof_handler.end();

    //unsigned int darcy_face, stokes_face;
    double point_flux;

    for(; lambda_cell != lambda_ecell; ++lambda_cell)
    {
      //Reinit lambda
      lambda_fe_values.reinit(lambda_cell);
      //CHECK: Modification in this line from using a mixed type variable to a native type
      darcy_cell  = darcy_lambda_map.get_source_pair(lambda_cell).first;
      stokes_cell = stokes_lambda_map.get_source_pair(lambda_cell).first;

      //Reinit Stokes and Darcy
      stokes_field.set_active_cell(stokes_cell);
      darcy_field.set_active_cell(darcy_cell);

      face_q_points = lambda_fe_values.get_quadrature_points();
      stokes_field.vector_value_list(face_q_points, stokes_face_values);
      darcy_field.vector_value_list(face_q_points,  darcy_face_values);

      lambda_cell->get_dof_indices(local_dof_indices);

      for(unsigned int q = 0; q < n_face_q_points; ++q)
      {
        point_flux = ( darcy_face_values[q](dim-1) - stokes_face_values[q](dim-1) );
        lambda_target(local_dof_indices[q]) = point_flux;
      }//end_q_point


    }//end_for_each_cell

}//compute_flux_function
//---------------------------------------------------------------
