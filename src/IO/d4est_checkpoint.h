#ifndef D4EST_CHECKPOINT_H
#define D4EST_CHECKPOINT_H 

#include <sc.h>
#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_h5.h>
#include <d4est_mesh.h>

typedef struct {

  int check_at_walltime_minus;
  int check_every_n_amr_steps;
  int check_every_n_krylov_steps;
  int check_every_n_newton_steps;
  int check_id;
  
} d4est_checkpoint_options_t;

/* This file was automatically generated.  Do not edit! */
p4est_t *d4est_checkpoint_load_mesh(sc_MPI_Comm mpicomm,const char *checkpoint_prefix,p4est_connectivity_t **connectivity);
void d4est_checkpoint_save(int checkpoint_number,const char *checkpoint_prefix,p4est_t *p4est,d4est_mesh_data_t *storage,const char **dg_field_names,double **dg_fields);

#endif
