#define _GNU_SOURCE
#include <d4est_checkpoint.h>
#include <d4est_mesh.h>
#include <d4est_element_data.h>


d4est_checkpoint_options_t*
d4est_checkpoint_options_init
(
 const char* input_file
){
  d4est_checkpoint_options_t* checkpoint_options = P4EST_ALLOC(d4est_checkpoint_options_t,1);

  checkpoint_options->check_id = 0;
  checkpoint_options->check_at_walltime_minus = -1;
  checkpoint_options->check_every_n_amr_steps = -1;
  checkpoint_options->check_every_n_krylov_steps = -1;
  checkpoint_options->check_every_n_newton_steps = -1;

  /* d4est_checkpoint_options_parse_input_file(checkpoint_options,input_file);  */
}

void
d4est_checkpoint_save
(){

}

void
d4est_checkpoint_options_destroy
(
 d4est_checkpoint_options_t* checkpoint_options
){
  P4EST_FREE(checkpoint_options);
}

void
d4est_checkpoint_save_aux
(
 const char* checkpoint_prefix,
 p4est_t* p4est,
 d4est_elliptic_data_t* elliptic_data,
 d4est_mesh_geometry_storage_t* storage,
 int number_of_user_datasets,
 void(*user_fcn)(p4est_t*, d4est_elliptic_data_t*, const char*, void*),
 void* user_ctx
){

  char* p4est_file_name;
  asprintf(&p4est_file_name,"%s.p4est", checkpoint_prefix);
  
  int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
  d4est_mesh_get_array_of_degrees(p4est, (void*)deg_array, D4EST_INT);
  
  p4est_save(p4est_file_name,
             p4est,
             0/* do not save element data */);

  d4est_h5_create_file(p4est->mpirank, checkpoint_prefix);
  d4est_h5_create_dataset(p4est->mpirank, checkpoint_prefix, "degree", H5T_NATIVE_INT, p4est->local_num_quadrants);
  d4est_h5_create_dataset(p4est->mpirank, checkpoint_prefix, "x", H5T_NATIVE_DOUBLE, elliptic_data->local_nodes);
  d4est_h5_create_dataset(p4est->mpirank, checkpoint_prefix, "y", H5T_NATIVE_DOUBLE, elliptic_data->local_nodes); 
  d4est_h5_write_dataset(p4est->mpirank, checkpoint_prefix, "degree", H5T_NATIVE_INT, deg_array);
  d4est_h5_write_dataset(p4est->mpirank, checkpoint_prefix, "x", H5T_NATIVE_DOUBLE, &storage->xyz[0]);
  d4est_h5_write_dataset(p4est->mpirank, checkpoint_prefix, "y", H5T_NATIVE_DOUBLE, &storage->xyz[elliptic_data->local_nodes]);
#if (P4EST_DIM)==3
  d4est_h5_create_dataset(p4est->mpirank, checkpoint_prefix, "z", H5T_NATIVE_DOUBLE, elliptic_data->local_nodes);
  d4est_h5_write_dataset(p4est->mpirank, checkpoint_prefix, "z", H5T_NATIVE_DOUBLE, &storage->xyz[elliptic_data->local_nodes*2]);
#endif

  user_fcn(p4est, elliptic_data, checkpoint_prefix, user_ctx);

  P4EST_FREE(deg_array);
  free(p4est_file_name);
}


p4est_t*
d4est_checkpoint_load_mesh
(
 sc_MPI_Comm mpicomm,
 const char* checkpoint_prefix,
 p4est_connectivity_t** connectivity
){
  printf("[D4EST_CHECKPOINT]: Loading meshing from checkpoint file %s.p4est\n", checkpoint_prefix);

  char* p4est_file_name;
  asprintf(&p4est_file_name,"%s.p4est", checkpoint_prefix);
  
  p4est_t* p4est = p4est_load (p4est_file_name, mpicomm, 0, 0, NULL, connectivity);

  free(p4est_file_name);
  return p4est;
}
