template <int dim>
void SD<dim>::darcy_operator ()
{
  TimerOutput::Scope timer_section(computing_timer, "Darcy operator");

  darcy_matrix = 0;
  darcy_rhs    = 0;


  FEValues<dim> darcy_fe_values (darcy_fe, q_rule,
                                 update_values    |
                                 //update_gradients |
                                 update_quadrature_points  |
                                 update_JxW_values);

  FEValues<dim> porosity_fe_values (porosity_fe, q_rule,
                                      update_values);

  FEFaceValues<dim> darcy_fe_face_values (darcy_fe, face_q_rule,
                                          update_values    | update_normal_vectors |
                                          update_quadrature_points  | update_JxW_values);

  const unsigned int   dofs_per_cell   = darcy_fe.dofs_per_cell;

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  //const PressureRightHandSide<dim>  pressure_right_hand_side;
  //const PressureBoundaryValues<dim> pressure_boundary_values;

  //std::vector<double>               pressure_rhs_values (n_q_points);
  //std::vector<double>               boundary_values (n_face_q_points);
  std::vector<Tensor<2,dim> >       k_inverse_values (n_q_points);

  std::vector<double>               old_porosity_values (n_q_points);

  std::vector<Tensor<1,dim> >       phi_u (dofs_per_cell);
  std::vector<double>               div_phi_u (dofs_per_cell);
  std::vector<double>               phi_p (dofs_per_cell);

  const FEValuesExtractors::Vector  velocities (0);
  const FEValuesExtractors::Scalar  pressure (dim);

  typename DoFHandler<dim>::active_cell_iterator
  cell = darcy_dof_handler.begin_active(),
  endc = darcy_dof_handler.end();
  typename DoFHandler<dim>::active_cell_iterator
  porosity_cell = porosity_dof_handler.begin_active();

  for (; cell!=endc; ++cell, ++porosity_cell)
    {
      darcy_fe_values.reinit (cell);
      porosity_fe_values.reinit (porosity_cell);

      local_matrix = 0;
      //local_rhs = 0;

      porosity_fe_values.get_function_values (old_porosity_solution, old_porosity_values);

      //pressure_right_hand_side.value_list (darcy_fe_values.get_quadrature_points(),
                                           //pressure_rhs_values);
      k_inverse.value_list (darcy_fe_values.get_quadrature_points(),
                            k_inverse_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int k=0; k<dofs_per_cell; ++k)
            {
              phi_u[k]     = darcy_fe_values[velocities].value (k,q);
              div_phi_u[k] = darcy_fe_values[velocities].divergence (k,q);
              phi_p[k]     = darcy_fe_values[pressure].value (k,q);
            }
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              for (unsigned int j=0; j<=i; ++j)
                {
                  local_matrix(i,j) += (phi_u[i] * k_inverse_values[q] * phi_u[j]
                                        + grad_div * div_phi_u[i] * div_phi_u[j]
                                        - div_phi_u[i] * phi_p[j]
                                        - phi_p[i] * div_phi_u[j])
                                       * darcy_fe_values.JxW(q);
                }

              //local_rhs(i) += (-phi_p[i] * pressure_rhs_values[q])*
                              //darcy_fe_values.JxW(q);
            }
        }

      //for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
        //if (cell->at_boundary(face_no))
          //{
            //darcy_fe_face_values.reinit (cell, face_no);
//
            //pressure_boundary_values
            //.value_list (darcy_fe_face_values.get_quadrature_points(),
                         //boundary_values);
//
            //for (unsigned int q=0; q<n_face_q_points; ++q)
              //for (unsigned int i=0; i<dofs_per_cell; ++i)
                //{
                  //const Tensor<1,dim>
                  //phi_i_u = darcy_fe_face_values[velocities].value (i, q);
//
                  //local_rhs(i) += -(phi_i_u *
                                    //darcy_fe_face_values.normal_vector(q) *
                                    //boundary_values[q] *
                                    //darcy_fe_face_values.JxW(q));
                //}
          //}

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=i+1; j<dofs_per_cell; ++j)
          local_matrix(i,j) = local_matrix(j,i);

      cell->get_dof_indices (local_dof_indices);

      darcy_constraints.distribute_local_to_global (local_matrix,
                                                    local_rhs,
                                                    local_dof_indices,
                                                    darcy_matrix,
                                                    darcy_rhs);

    }//for_each_cell
}
