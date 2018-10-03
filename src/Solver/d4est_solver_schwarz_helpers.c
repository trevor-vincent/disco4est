/** 
 * Iterate over a field and restrict it to each
 * subdomain
 * 
 * @param [in] schwarz_data 
 * @param [in] field 
 * @param [in] ghost_data 
 * @param [in] field_num 
 * @param [out] field_over_subdomains 
 */

void
d4est_solver_schwarz_restrict_nodal_field_to_subdomains
(
 d4est_schwarz_data_t* schwarz_data,
 d4est_factors_t* d4est_factors,
 double* field,
 d4est_ghost_data_t* ghost_data,
 int ghost_data_num_of_field,
 double* restricted_field_over_subdomains
)
{
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[i];
    for (int j = 0; j < sub_data->num_elements; j++){
      d4est_solver_schwarz_element_data_t ed = sub_data->element_data[j];

      double* field_ed = NULL;      
      if (ed->mpi_rank == p4est->mpirank){
        field_ed = &field[d4est_factors->element_data[ed->id]->nodal_stride]
      }
      else {
        field_ed =
          d4est_ghost_data_get_field_on_element
          (
           d4est_factors->element_data[ed->id + p4est->local_num_quadrants],
           ghost_data_num_of_field,
           ghost_data
          );
      }
      d4est_operators_apply_schwarz_restrictor
        (
         d4est_ops,
         field_ed,
         (P4EST_DIM),
         &(sub_e->faces[0]),
         ed->deg,
         schwarz_data->overlap_size,
         D4OPS_NO_TRANSPOSE,
         &restricted_field_over_subdomains[sub_data->restricted_nodal_stride
                                           + ed->restricted_nodal_stride]
        );
    }
  }
}


/** 
 * 
 * 
 * @param schwarz_data 
 * @param field_over 
 */
void
d4est_solver_schwarz_restrict_transpose_restricted_field_over_subdomains
(
 d4est_schwarz_data_t* schwarz_data,
 double* restricted_field_over_subdomains,
 double* field_over_subdomains
)
{
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[i];
    for (int j = 0; j < sub_data->num_elements; j++){
      d4est_solver_schwarz_element_data_t ed = sub_data->element_data[j];

      field_ed = &field[d4est_factors->element_data[ed->id]->nodal_stride]
                 d4est_operators_apply_schwarz_restrictor
                 (
                  d4est_ops,
                  &restricted_field_over_subdomains[sub_data->restricted_nodal_stride
                                         + ed->restricted_nodal_stride]
                  (P4EST_DIM),
                  &(sub_e->faces[0]),
                  ed->deg,
                  schwarz_data->overlap_size,
                  D4OPS_TRANSPOSE,
                  &field_over_subdomains[sub_data->nodal_stride
                                         + ed->nodal_stride]
                 );
    }
  }
}

void
d4est_solver_schwarz_restrict_field_over_subdomains
(
 d4est_schwarz_data_t* schwarz_data,
 double* field_over_subdomains,
 double* restricted_field_over_subdomains
)
{
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[i];
    for (int j = 0; j < sub_data->num_elements; j++){
      d4est_solver_schwarz_element_data_t ed = sub_data->element_data[j];

      field_ed = &field[d4est_factors->element_data[ed->id]->nodal_stride]
                 d4est_operators_apply_schwarz_restrictor
                 (
                  d4est_ops,
                  &field_over_subdomains[sub_data->nodal_stride
                                         + ed->nodal_stride],
                  (P4EST_DIM),
                  &(sub_e->faces[0]),
                  ed->deg,
                  schwarz_data->overlap_size,
                  D4OPS_NO_TRANSPOSE,
                  &restricted_field_over_subdomains[sub_data->restricted_nodal_stride
                                                    + ed->restricted_nodal_stride]

                 );
    }
  }
}

void
d4est_solver_schwarz_compute_correction
(
 d4est_solver_schwarz_data_t* schwarz_data,
 double* du_restricted_field_over_subdomains,
 double* correction_field_over_subdomains
)
{
  double* weighted_du_restricted_field_over_subdomains
    = P4EST_ALLOC(double, schwarz_data->restricted_nodal_size);
  
  d4est_solver_schwarz_apply_weights_over_all_subdomains
    (
     schwarz_data,
     du_restricted_field_over_subdomains,
     weighted_du_restricted_field_over_subdomains
    );

  d4est_solver_schwarz_restrict_transpose_restricted_field_over_subdomains
    (
     schwarz_data,
     weighted_du_restricted_field_over_subdomains,
     correction_field_over_subdomains
    );
  
  P4EST_FREE(weighted_du_restricted_field_over_subdomains);
}

/** 
 * Take a single subdomain field on a certain process
 * and turn it into a global nodal field (e.g. the parts
 * on ghost elements will be returned to their given 
 * processes). This is useful for visualization.
 * This uses communication.
 * 
 */
void
d4est_solver_schwarz_convert_subdomain_field_to_global_nodal_field
(
 p4est_t* p4est,
 d4est_schwarz_data_t* schwarz_data,
 double* field_over_subdomains, /* not restricted */
 double* field,
 int mpirank,
 int sub_id,
 int local_nodes
)
{
  /* transfer subdomains if multi-process */
  /* copy over relevant data if process matches and element matches, otherwise zero */

  for (int i = 0; i < local_nodes; i++){
    field[i] = 0.;
  }

  /* Don't need to check if p4est->mpirank == mpirank in this case */
  if (p4est->mpisize == 1){    
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[sub_id];
    for (int j = 0; j < sub_data->num_elements; j++){
      d4est_solver_schwarz_element_data_t ed = sub_data->element_data[j];
      int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
      int local_nodal_stride = d4est_factors->element_data[ed->id]->nodal_stride;
      int sub_nodal_stride = sub_data->subdomain_nodal_stride + ed->nodal_stride;

      d4est_util_copy_1st_to_second(&field_over_subdomains[sub_nodal_stride],
                                    &field[local_nodal_stride]
                                    volume_nodes);
      
    }
  }
  else {
    /* TRANSER STUFF HERE BEFORE IF */
    if (p4est->mpirank == mpirank){
      /* GO THROUGH LOCAL SUBDOMAINS AND COPY DATA */
    }
    else {
      /* GO THROUGH GHOST SUBDOMAIN DATA AND COPY DATA */
    }
  }
  
}

