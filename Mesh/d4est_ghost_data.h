#ifndef D4EST_GHOST_DATA_H
#define D4EST_GHOST_DATA_H 

#include <d4est_ghost.h>
#include <d4est_field.h>

typedef struct {

  int num_ghosts;
  int num_vecs;
  d4est_field_type_t* transfer_types;
  
  int* ghost_data_sizes;  
  int receive_size;
  int** receive_strides;
  double* receive_data;  

} d4est_ghost_data_t;

/* This file was automatically generated.  Do not edit! */
void d4est_ghost_data_exchange(p4est_t *p4est,d4est_ghost_t *d4est_ghost,d4est_ghost_data_t *d4est_ghost_data,double **transfer_vecs);
void d4est_ghost_data_destroy(d4est_ghost_data_t *d4est_ghost_data);
d4est_ghost_data_t *d4est_ghost_data_init(p4est_t *p4est,d4est_ghost_t *d4est_ghost,d4est_field_type_t *field_types,int num_vecs);


#endif
