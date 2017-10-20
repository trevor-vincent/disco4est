#ifndef D4EST_H5_H
#define D4EST_H5_H 

#include <hdf5.h>

/* This file was automatically generated.  Do not edit! */
void d4est_h5_read_dataset(int mpirank,const char *file_name_prefix,const char *dataset_name,hid_t dataset_type,void *dataset);
void d4est_h5_write_dataset(int mpirank,const char *file_name_prefix,const char *dataset_name,hid_t dataset_type,void *dataset);
void d4est_h5_create_dataset(int mpirank,const char *file_name_prefix,const char *dataset_name,hid_t dataset_type,int size);
void d4est_h5_create_file(int mpirank,const char *file_name_prefix);

#endif
