#ifndef D4EST_GEOMETRY_BRICK_H
#define D4EST_GEOMETRY_BRICK_H 

#include <d4est_geometry.h>

typedef struct
{
  double X0, X1, Y0, Y1, Z0, Z1;
  char input_section [50];
  
} d4est_geometry_brick_attr_t;

/* This file was automatically generated.  Do not edit! */
void d4est_geometry_brick_new(int mpirank,const char *input_file,const char *input_section,const char *printf_prefix, d4est_geometry_t *d4est_geom);

#endif
