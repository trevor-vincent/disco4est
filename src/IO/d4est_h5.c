#define _GNU_SOURCE
#include <d4est_h5.h>
#include <hdf5.h>
#include <stdlib.h>

void
d4est_h5_create_file
(
 int mpirank,
 const char* file_name_prefix
){

  char* file_name;
  asprintf(&file_name,"%s_%d.h5", file_name_prefix, mpirank);
  
  /* Create a new file using default properties. */
  hid_t file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /* Close the file. */
  H5Fclose(file_id);
  
  free(file_name);
   
}

void
d4est_h5_create_dataset
(
 int mpirank,
 const char* file_name_prefix,
 const char* dataset_name,
 hid_t dataset_type,
 int size
){

  char* file_name;
  asprintf(&file_name,"%s_%d.h5", file_name_prefix, mpirank);

  char* dataset_path;
  asprintf(&dataset_path,"/%s", dataset_name);
  
  /* Create a new file using default properties. */
  hid_t file_id = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);

  /* Create the data space for the dataset. */
  hsize_t dims [] = {size}; 
  hid_t dataspace_id = H5Screate_simple(1, dims, NULL);

  /* Create the dataset. */
  hid_t dataset_id = H5Dcreate2(file_id, dataset_path, dataset_type, dataspace_id, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* End access to the dataset and release resources used by it. */
  H5Dclose(dataset_id);

  /* Terminate access to the data space. */ 
  H5Sclose(dataspace_id);

  /* Close the file. */
  H5Fclose(file_id);

  free(file_name);
  free(dataset_path);
}

void
d4est_h5_write_dataset
(
 int mpirank,
 const char* file_name_prefix,
 const char* dataset_name,
 hid_t dataset_type,
 void* dataset
)
{
  char* file_name;
  asprintf(&file_name,"%s_%d.h5", file_name_prefix, mpirank);

  char* dataset_path;
  asprintf(&dataset_path,"/%s", dataset_name);
  
  /* Open an existing file. */
  hid_t file_id = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);

  /* Open an existing dataset. */
  hid_t dataset_id = H5Dopen2(file_id, dataset_path, H5P_DEFAULT);

  /* Write the dataset. */
  H5Dwrite(dataset_id, dataset_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset);
  /* H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset); */

  /* Close the dataset. */
  H5Dclose(dataset_id);

  /* Close the file. */
  H5Fclose(file_id);

  free(file_name);
  free(dataset_path);
}

void
d4est_h5_read_dataset
(
 int mpirank,
 const char* file_name_prefix,
 const char* dataset_name,
 hid_t dataset_type,
 void* dataset
)
{
  char* file_name;
  asprintf(&file_name,"%s_%d.h5", file_name_prefix, mpirank);

  char* dataset_path;
  asprintf(&dataset_path,"/%s", dataset_name);
  
  /* Open an existing file. */
  hid_t file_id = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);

  /* Open an existing dataset. */
  hid_t dataset_id = H5Dopen2(file_id, dataset_path, H5P_DEFAULT);

  /* Write the dataset. */
  /* H5Dwrite(dataset_id, dataset_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset); */
  H5Dread(dataset_id, dataset_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset);

  /* Close the dataset. */
  H5Dclose(dataset_id);

  /* Close the file. */
  H5Fclose(file_id);

  free(file_name);
  free(dataset_path);
}
