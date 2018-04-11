#ifndef D4EST_GHOST_H
#define D4EST_GHOST_H 

#include <pXest.h>
#include <d4est_element_data.h>

typedef struct {

  p4est_ghost_t* ghost;
  d4est_element_data_t* ghost_elements;
  d4est_element_data_t** mirror_elements;
  
} d4est_ghost_t;

/* This file was automatically generated.  Do not edit! */
void d4est_ghost_update(p4est_t *p4est, d4est_ghost_t *d4est_ghost);
void d4est_ghost_destroy(d4est_ghost_t *d4est_ghost);
d4est_ghost_t *d4est_ghost_init(p4est_t *p4est);

#endif
