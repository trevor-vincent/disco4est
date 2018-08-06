#ifndef D4EST_CHECKPOINT_H
#define D4EST_CHECKPOINT_H 

#include <sc.h>
#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_h5.h>
#include <d4est_amr.h>
#include <d4est_mesh.h>


/* This file was automatically generated.  Do not edit! */
void d4est_checkpoint_check_dataset(p4est_t *p4est,const char *checkpoint_prefix,const char *dataset_name,hid_t dataset_type,void *dataset_sum,int checkpoint_number);
void d4est_checkpoint_read_dataset(p4est_t *p4est,const char *checkpoint_prefix,const char *dataset_name,hid_t dataset_type,void *dataset,int checkpoint_number);
void d4est_checkpoint_load_mesh_from_amr_history(sc_MPI_Comm mpicomm,p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_mesh_initial_extents_t *initial_grid_input);
p4est_t *d4est_checkpoint_load_p4est_from_file(sc_MPI_Comm mpicomm,const char *checkpoint_prefix,p4est_connectivity_t **connectivity);
void d4est_checkpoint_save(int checkpoint_number,const char *checkpoint_prefix,p4est_t *p4est,d4est_amr_t *d4est_amr,d4est_mesh_data_t *storage,const char **dg_field_names,double **dg_fields,const char **element_field_names,double **element_fields);

#endif
