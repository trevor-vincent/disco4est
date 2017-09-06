#ifndef D4EST_GEOMETRY_BRICK_H
#define D4EST_GEOMETRY_BRICK_H 

#include <d4est_geometry.h>

typedef struct
{
  double X0, X1, Y0, Y1, Z0, Z1;
  char input_section [50];
  
} d4est_geometry_brick_attr_t;

/* This file was automatically generated.  Do not edit! */
int d4est_geometry_brick_which_child_of_root(p4est_qcoord_t q[(P4EST_DIM)],p4est_qcoord_t dq);
void d4est_geometry_brick_new(int mpirank,const char *input_file,const char *input_section,const char *printf_prefix,d4est_geometry_t *d4est_geom);
void d4est_geometry_brick_new_aux(d4est_geometry_t *d4est_geom,d4est_geometry_brick_attr_t *brick_attrs);
 
#endif
