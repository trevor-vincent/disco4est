#define _GNU_SOURCE
#include <stdio.h>
#include <d4est_checkpoint.h>
#include <d4est_h5.h>
#include <d4est_mesh.h>
#include <d4est_element_data.h>

void
d4est_checkpoint_save
(
 int checkpoint_number,
 const char* checkpoint_prefix,
 p4est_t* p4est,
 d4est_amr_t* d4est_amr,
 d4est_mesh_data_t* storage,
 const char ** field_names,
 hid_t* field_data_types,
 int* field_sizes,
 void** fields
)
{
  zlog_category_t* c_default = zlog_get_category("d4est_checkpoint");
  char* checkpoint_folder = d4est_util_add_cwd("Checkpoints");
  d4est_util_make_directory(checkpoint_folder,0);
  
  if (checkpoint_number >= 0){
    asprintf(&checkpoint_folder,"%s%d/", checkpoint_folder, checkpoint_number);
    d4est_util_make_directory(checkpoint_folder,0);
  }
  
  char* checkpoint_folder_and_prefix = NULL;
  asprintf(&checkpoint_folder_and_prefix,"%s%s", checkpoint_folder, checkpoint_prefix);
  
  char* p4est_file_name;
  asprintf(&p4est_file_name,"%s.p4est", checkpoint_folder_and_prefix);
  
  int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
  d4est_mesh_get_array_of_degrees(p4est, (void*)deg_array, D4EST_INT);
  
  /* p4est_save(p4est_file_name, */
             /* p4est, */
             /* 0/\* do not save element data *\/); */

  int num_dg_nodes = 0;
  for (int i = 0; i < p4est->local_num_quadrants; i++){
    num_dg_nodes += d4est_lgl_get_nodes((P4EST_DIM), deg_array[i]);
  }
  
  d4est_h5_create_file(p4est->mpirank, checkpoint_folder_and_prefix);
  d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "degree", H5T_NATIVE_INT, p4est->local_num_quadrants);
  d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "degree", H5T_NATIVE_INT, deg_array);

  if (storage != NULL){
    d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "x", H5T_NATIVE_DOUBLE, num_dg_nodes);
    d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "y", H5T_NATIVE_DOUBLE, num_dg_nodes); 
    d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "x", H5T_NATIVE_DOUBLE, &storage->xyz[0]);
    d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "y", H5T_NATIVE_DOUBLE, &storage->xyz[num_dg_nodes]);
#if (P4EST_DIM)==3
    d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "z", H5T_NATIVE_DOUBLE, num_dg_nodes);
    d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "z", H5T_NATIVE_DOUBLE, &storage->xyz[num_dg_nodes*2]);
#endif
  }
  /* just in case the saved p4est file is written out wrong, which has been the case, we 
   * will save the refinement log, which can be used to restore the forest if the 
   * refinement log at each refinement level is saved (e.g. you checkpoint at each level)  */
  if (d4est_amr != NULL){
    d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "refinement_log", H5T_NATIVE_INT, d4est_amr->initial_log_size);
    d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "refinement_log", H5T_NATIVE_INT, d4est_amr->refinement_log); 

    int refinement_log_sum = d4est_util_sum_array_int(d4est_amr->refinement_log, d4est_amr->initial_log_size);
    d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "refinement_log_sum", H5T_NATIVE_INT, 1);
    d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "refinement_log_sum", H5T_NATIVE_INT, &refinement_log_sum); 
  }
  /* useful to check that restored d4est mesh has not changed from checkpoint file */
  int checkdeg = d4est_mesh_get_local_degree_sum(p4est);
  int checkp4est = p4est_checksum(p4est);
  int checkmpi = p4est->mpisize;
  int checkquadrants = p4est->local_num_quadrants;

  d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "p4est_checksum", H5T_NATIVE_INT, 1);
  d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "p4est_checksum", H5T_NATIVE_INT, &checkp4est); 
  d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "degree_sum", H5T_NATIVE_INT, 1);
  d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "degree_sum", H5T_NATIVE_INT, &checkdeg); 
  d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "mpi_size", H5T_NATIVE_INT, 1);
  d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "mpi_size", H5T_NATIVE_INT, &checkmpi); 
  d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "quadrants", H5T_NATIVE_INT, 1);
  d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "quadrants", H5T_NATIVE_INT, &checkquadrants); 

  int num_fields = 0;
  if (field_names != NULL && fields != NULL){
    for (int i = 0; field_names[i] != NULL; i++){
      num_fields++;
    }
  }

  
  for (int i = 0; i < num_fields; i++){
    d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, field_names[i], field_data_types[i], field_sizes[i]);
    d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, field_names[i], field_data_types[i], fields[i]);

    /* save field sums to do checks of field/file integrity */
    if (field_sizes[i] > 0){
      char* field_sum_name = NULL;
      asprintf(&field_sum_name,"%s_sum",field_names[i]);
      double field_sum = d4est_util_sum_array_dbl(fields[i], field_sizes[i]);
      d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, field_sum_name, field_data_types[i], 1);
      /* d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, field_sum_name, H5T_NATIVE_DOUBLE, &field_sum); */

      if (field_data_types[i] == H5T_NATIVE_DOUBLE){
        double field_sum = d4est_util_sum_array_dbl(fields[i], field_sizes[i]);
        d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, field_sum_name, field_data_types[i], &field_sum);
      }
      else if (field_data_types[i] == H5T_NATIVE_INT) {
        int field_sum = d4est_util_sum_array_int(fields[i], field_sizes[i]);
        d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, field_sum_name, field_data_types[i], &field_sum);
      }
      else {
        zlog_info(c_default, "not a supported hdf5 type");
        D4EST_ABORT("");
      }
      free(field_sum_name);
    }
  }
  
  
  /* int num_dg_fields = 0; */
  /* if (dg_field_names != NULL && dg_fields != NULL){ */
  /*   for (int i = 0; dg_field_names[i] != NULL; i++){ */
  /*     num_dg_fields++; */
  /*   } */
  /* } */
  
  /* for (int i = 0; i < num_dg_fields; i++){ */
  /*   d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, dg_field_names[i], H5T_NATIVE_DOUBLE, num_dg_nodes); */
  /*   d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, dg_field_names[i], H5T_NATIVE_DOUBLE, dg_fields[i]); */
  /*   /\* save field sums to do checks of field/file integrity *\/ */
  /*   char* field_sum_name = NULL; */
  /*   asprintf(&field_sum_name,"%s_sum",dg_field_names[i]); */
  /*   double field_sum = d4est_util_sum_array_dbl(dg_fields[i], num_dg_nodes); */
  /*   d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, field_sum_name, H5T_NATIVE_DOUBLE, 1); */
  /*   d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, field_sum_name, H5T_NATIVE_DOUBLE, &field_sum); */
  /*   free(field_sum_name); */
  /* } */

  /* int num_element_fields = 0; */
  /* if (element_field_names != NULL && element_fields != NULL){ */
  /*   for (int i = 0; element_field_names[i] != NULL; i++){ */
  /*     num_element_fields++; */
  /*   } */
  /* } */

  /* int num_element_nodes = p4est->local_num_quadrants; */
  /* for (int i = 0; i < num_element_fields; i++){ */
  /*   d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, element_field_names[i], element_field_types[i], num_element_nodes); */
  /*   d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, element_field_names[i], element_field_types[i], element_fields[i]); */

  /*   /\* save field sums to do checks of field/file integrity *\/ */
  /*   char* field_sum_name = NULL; */
  /*   asprintf(&field_sum_name,"%s_sum",element_field_names[i]); */
  /*   d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, field_sum_name, element_field_types[i], 1); */

  /*   if (element_field_types[i] == H5T_NATIVE_DOUBLE){ */
  /*     double field_sum = d4est_util_sum_array_dbl(element_fields[i], p4est->local_num_quadrants); */
  /*     d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, field_sum_name, element_field_types[i], &field_sum); */
  /*   } */
  /*   else if (element_field_types[i] == H5T_NATIVE_INT) { */
  /*     int field_sum = d4est_util_sum_array_int(element_fields[i], p4est->local_num_quadrants); */
  /*   d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, field_sum_name, element_field_types[i], &field_sum); */
  /*   } */
  /*   else { */
  /*     zlog_info(c_default, "not a supported hdf5 type"); */
  /*     D4EST_ABORT(""); */
  /*   } */
    
  /*   free(field_sum_name); */
  /* } */

  /* int num_constant = 0; */
  /* if (constant_names != NULL && constant != NULL){ */
  /*   for (int i = 0; constant_names[i] != NULL; i++){ */
  /*     num_constant++; */
  /*   } */
  /* } */

  /* int num_constant_nodes = 1; */
  /* for (int i = 0; i < num_constant; i++){ */
  /*   d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, constant_names[i], constant_types[i], num_constant_nodes); */
  /*   d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, constant_names[i], constant_types[i], constant[i]); */
  /* } */


  free(checkpoint_folder);
  free(checkpoint_folder_and_prefix);
  free(p4est_file_name);
  P4EST_FREE(deg_array);
}


p4est_t*
d4est_checkpoint_load_p4est_from_file
(
 sc_MPI_Comm mpicomm,
 const char* checkpoint_prefix,
 p4est_connectivity_t** connectivity
)
{
  printf("[D4EST_CHECKPOINT]: Loading meshing from checkpoint file %s.p4est\n", checkpoint_prefix);

  char* p4est_file_name;
  asprintf(&p4est_file_name,"%s.p4est", checkpoint_prefix);
  
  p4est_t* p4est = p4est_load (p4est_file_name, mpicomm, 0, 0, NULL, connectivity);

  free(p4est_file_name);
  return p4est;
}


static void 
d4est_checkpoint_load_mesh_from_amr_history_mark_elements
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  d4est_amr_t* d4est_amr = (d4est_amr_t*) info->p4est->user_pointer;
  int* refinement_log = (int*) (d4est_amr->scheme->amr_scheme_data);
  d4est_element_data_t* elem_data = (d4est_element_data_t*) info->quad->p.user_data;
  d4est_amr->refinement_log[elem_data->id] = refinement_log[elem_data->id];
  /* printf("refinement_log[%d] = %d\n", elem_data->id, refinement_log[elem_data->id]); */
}

/** 
 * Requires p4est to be initialized
 * as it was at the beginning of the
 * simulation. e.g. don't change initial
 * mesh part of the options input file
 * if you want to load a mesh this way
 * this will also set the degrees.
 * Make sure the checkpoints are in a
 * Checkpoints subfolder and labeled
 * 0 ... checkpoint_history_size - 1
 * This function will iteratively refine the mesh
 * until it reaches the final configuration
 * 
 * @param mpicomm 
 * @param checkpoint_prefix 
 * @param p4est 
 */
void
d4est_checkpoint_load_mesh_from_amr_history
(
 sc_MPI_Comm mpicomm,
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_mesh_initial_extents_t* initial_grid_input
)
{
  int initial_checkpoint_number = initial_grid_input->initial_checkpoint_number;
  int final_checkpoint_number = initial_grid_input->checkpoint_number;
  const char* checkpoint_prefix = initial_grid_input->checkpoint_prefix;

  
  for (int i = initial_checkpoint_number; i < final_checkpoint_number + 1; i++){
    
    int* refinement_log = P4EST_ALLOC(int, p4est->local_num_quadrants);
    
    d4est_checkpoint_read_dataset
      (
       p4est,
       checkpoint_prefix,
       "refinement_log",
       H5T_NATIVE_INT,
       refinement_log,
       i
      );

    int refinement_log_sum = d4est_util_sum_array_int(refinement_log, p4est->local_num_quadrants);

    d4est_checkpoint_check_dataset
      (
       p4est,
       checkpoint_prefix,
       "refinement_log",
       H5T_NATIVE_INT,
       (void*)&refinement_log_sum,
       i
      );
    
    /* create amr scheme */
    d4est_amr_t* d4est_amr = P4EST_ALLOC(d4est_amr_t, 1);
    d4est_amr_scheme_t* scheme = P4EST_ALLOC(d4est_amr_scheme_t, 1);
    scheme->post_balance_callback = NULL;
    scheme->pre_refine_callback = NULL;
    scheme->refine_replace_callback_fcn_ptr = NULL;
    scheme->balance_replace_callback_fcn_ptr = NULL;
    scheme->mark_elements = d4est_checkpoint_load_mesh_from_amr_history_mark_elements;
    scheme->amr_scheme_data = refinement_log;
    scheme->destroy = NULL;
    
    d4est_amr->mpirank = p4est->mpirank;
    d4est_amr->scheme = scheme;
    d4est_amr->balance_log = NULL;
    d4est_amr->refinement_log = NULL;
    d4est_amr->initial_log = NULL;
    d4est_amr->max_degree = 1000;
    
    d4est_amr_step
      (
       p4est,
       NULL,
       d4est_amr,
       NULL,
       NULL,
       NULL
      );

    d4est_mesh_update
      (
       p4est,
       NULL,
       d4est_ops,
       d4est_geom,
       NULL,
       NULL,
       initial_grid_input,
       DO_NOT_INITIALIZE_GHOST,
       DO_NOT_INITIALIZE_QUADRATURE_DATA,
       DO_NOT_INITIALIZE_GEOMETRY_DATA,
       DO_NOT_INITIALIZE_GEOMETRY_ALIASES,
       d4est_mesh_set_quadratures_after_amr,
       initial_grid_input
      );
    

    int quadrants;
    int mpi_size;
    d4est_checkpoint_read_dataset
      (
       p4est,
       checkpoint_prefix,
       "quadrants",
       H5T_NATIVE_INT,
       &quadrants,
       i
      );

    d4est_checkpoint_read_dataset
      (
       p4est,
       checkpoint_prefix,
       "mpi_size",
       H5T_NATIVE_INT,
       &mpi_size,
       i
      );

    printf("p4est->mpisize, mpi_size from file = %d, %d\n", p4est->mpisize, mpi_size);
    printf("local num quadrants from p4est, from file = %d, %d\n", p4est->local_num_quadrants, quadrants);
    
    if(!(p4est->local_num_quadrants == quadrants && p4est->mpisize == mpi_size)){
      D4EST_ABORT("!(p4est->local_num_quadrants == quadrants && p4est->mpisize == mpi_size)");
    }

  }
}

void
d4est_checkpoint_read_dataset
(
 p4est_t* p4est,
 const char* checkpoint_prefix,
 const char* dataset_name,
 hid_t dataset_type,
 void* dataset,
 int checkpoint_number
)
{
  char* checkpoint_folder = d4est_util_add_cwd("Checkpoints");
  asprintf(&checkpoint_folder,"%s%d/", checkpoint_folder, checkpoint_number);
  
  char* checkpoint_folder_and_prefix = NULL;
  asprintf(&checkpoint_folder_and_prefix,"%s%s", checkpoint_folder, checkpoint_prefix);
  
  d4est_h5_read_dataset(p4est->mpirank,
                        checkpoint_folder_and_prefix,
                        dataset_name,
                        dataset_type,
                        dataset);

  free(checkpoint_folder_and_prefix);
  free(checkpoint_folder);
}


void
d4est_checkpoint_check_dataset
(
 p4est_t* p4est,
 const char* checkpoint_prefix,
 const char* dataset_name,
 hid_t dataset_type,
 void* dataset_sum,
 int checkpoint_number
){
  zlog_category_t *c_default = zlog_get_category("d4est_h5");
  char* dataset_sum_name = NULL;
  asprintf(&dataset_sum_name,"%s_sum",dataset_name);

  char* checkpoint_folder = d4est_util_add_cwd("Checkpoints");
  asprintf(&checkpoint_folder,"%s%d/", checkpoint_folder, checkpoint_number);
  
  char* checkpoint_folder_and_prefix = NULL;
  asprintf(&checkpoint_folder_and_prefix,"%s%s", checkpoint_folder, checkpoint_prefix);

  double dataset_sum_dbl_check = -1;
  int dataset_sum_int_check = -1;
  
  if (dataset_type == H5T_NATIVE_DOUBLE){
    d4est_h5_read_dataset(
                          p4est->mpirank,
                          checkpoint_folder_and_prefix,
                          dataset_sum_name,
                          dataset_type,
                          (void*)&dataset_sum_dbl_check
    );
  }
  else if (dataset_type == H5T_NATIVE_INT){
    d4est_h5_read_dataset(
                          p4est->mpirank,
                          checkpoint_folder_and_prefix,
                          dataset_sum_name,
                          dataset_type,
                          (void*)&dataset_sum_int_check
    );
  }
  else {
    zlog_error(c_default, "Not a supported Hdf5 datase2t type");
    zlog_error(c_default, "type = %d", (int)dataset_type);
    D4EST_ABORT("");
  }

  
  
  free(checkpoint_folder_and_prefix);
  free(checkpoint_folder);
  free(dataset_sum_name);  
  if (dataset_type == H5T_NATIVE_DOUBLE){
    double dataset_sum_dbl = *(double*)dataset_sum;
    if (fabs(dataset_sum_dbl - dataset_sum_dbl_check) > 1.e-14){
      zlog_error(c_default, "dataset = %s", dataset_name);
      zlog_error(c_default, "dataset_sum_dbl = %.20f", dataset_sum_dbl);
      zlog_error(c_default, "dataset_sum_dbl_check = %.20f", dataset_sum_dbl_check);
      zlog_error(c_default, "fabs(dataset_sum_dbl - dataset_sum_dbl_check) = %.20f", fabs(dataset_sum_dbl - dataset_sum_dbl_check));
      D4EST_ABORT("");
    }

  }
  else if (dataset_type == H5T_NATIVE_INT){
    int dataset_sum_int = *(int*)dataset_sum;
    if (dataset_sum_int != dataset_sum_int_check){
      zlog_error(c_default, "dataset = %s", dataset_name);
      zlog_error(c_default, "dataset_sum_int = %d", dataset_sum_int);
      zlog_error(c_default, "dataset_sum_int_check = %d", dataset_sum_int_check);
      zlog_error(c_default, "dataset_sum_int - dataset_sum_int_check = %d", dataset_sum_int - dataset_sum_int_check);
      D4EST_ABORT("");
    }
  }
  else {
    zlog_error(c_default, "Not a supported Hdf5 dataset type");
    zlog_error(c_default, "type = %d", (int)dataset_type);
    D4EST_ABORT("");
  }

}
