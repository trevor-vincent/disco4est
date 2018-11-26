#ifndef D4EST_GHOST_DATA_EXT_H
#define D4EST_GHOST_DATA_EXT_H 


#include <d4est_ghost.h>
#include <d4est_field.h>


typedef
int
(*d4est_ghost_data_ext_fcn_t)
(
 d4est_ghost_t* d4est_ghost,
 int ghost_id,
 int vec_id,
 void*
);

typedef struct {

  int num_ghosts;
  int num_vecs;  
  int* ghost_data_sizes;  
  int receive_size; /* total ghost nodes over all num_vecs */
  int** receive_strides;
  char* receive_data;  

  d4est_ghost_data_ext_fcn_t field_size_of_ghost_fcn;
  d4est_ghost_data_ext_fcn_t field_size_of_mirror_fcn;
  d4est_ghost_data_ext_fcn_t field_stride_of_mirror_fcn;
  void* user_ctx;
  
} d4est_ghost_data_ext_t;

/* This file was automatically generated.  Do not edit! */
void *d4est_ghost_data_ext_get_field_on_element(d4est_element_data_t *ed,int ghost_name_id,d4est_ghost_data_ext_t *dgd);
void d4est_ghost_data_ext_exchange(p4est_t *p4est,d4est_ghost_t *d4est_ghost,d4est_ghost_data_ext_t *d4est_ghost_data_ext,char **transfer_vecs);
void d4est_ghost_data_ext_destroy(d4est_ghost_data_ext_t *d4est_ghost_data_ext);
d4est_ghost_data_ext_t *d4est_ghost_data_ext_init(p4est_t *p4est,d4est_ghost_t *d4est_ghost,int num_vecs,d4est_ghost_data_ext_fcn_t field_size_of_ghost_fcn,d4est_ghost_data_ext_fcn_t field_size_of_mirror_fcn,d4est_ghost_data_ext_fcn_t field_stride_of_mirror_fcn,void *user_ctx);

#endif
