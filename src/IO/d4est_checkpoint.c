#define _GNU_SOURCE
#include <d4est_checkpoint.h>
#include <d4est_mesh.h>
#include <d4est_element_data.h>


void
d4est_checkpoint_save
(
 int checkpoint_number,
 const char* checkpoint_prefix,
 p4est_t* p4est,
 d4est_mesh_data_t* storage,
 const char ** dg_field_names,
 double ** dg_fields
)
{
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
  
  p4est_save(p4est_file_name,
             p4est,
             0/* do not save element data */);


  int num_dg_nodes = 0;
  for (int i = 0; i < p4est->local_num_quadrants; i++){
    num_dg_nodes += d4est_lgl_get_nodes((P4EST_DIM), deg_array[i]);
  }
  
  d4est_h5_create_file(p4est->mpirank, checkpoint_folder_and_prefix);
  d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "degree", H5T_NATIVE_INT, p4est->local_num_quadrants);
  d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "x", H5T_NATIVE_DOUBLE, num_dg_nodes);
  d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "y", H5T_NATIVE_DOUBLE, num_dg_nodes); 
  d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "degree", H5T_NATIVE_INT, deg_array);
  d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "x", H5T_NATIVE_DOUBLE, &storage->xyz[0]);
  d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "y", H5T_NATIVE_DOUBLE, &storage->xyz[num_dg_nodes]);
#if (P4EST_DIM)==3
  d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "z", H5T_NATIVE_DOUBLE, num_dg_nodes);
  d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, "z", H5T_NATIVE_DOUBLE, &storage->xyz[num_dg_nodes*2]);
#endif

  int num_dg_fields = 0;
  if (dg_field_names != NULL && dg_fields != NULL){
    for (int i = 0; dg_field_names[i] != NULL; i++){
      num_dg_fields++;
    }
  }

  
  for (int i = 0; i < num_dg_fields; i++){
    d4est_h5_create_dataset(p4est->mpirank, checkpoint_folder_and_prefix, dg_field_names[i], H5T_NATIVE_DOUBLE, num_dg_nodes);
    d4est_h5_write_dataset(p4est->mpirank, checkpoint_folder_and_prefix, dg_field_names[i], H5T_NATIVE_DOUBLE, dg_fields[i]);  
  }
  
  free(checkpoint_folder);
  free(checkpoint_folder_and_prefix);
  free(p4est_file_name);
  P4EST_FREE(deg_array);
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
