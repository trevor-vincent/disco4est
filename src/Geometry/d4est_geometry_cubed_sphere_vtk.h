#ifndef D4EST_GEOMETRY_CUBED_SPHERE_VTK_H
#define D4EST_GEOMETRY_CUBED_SPHERE_VTK_H 


#include <pXest.h>
#include <p8est_geometry.h>
#include <zlog.h>


typedef struct
{
  double R2_a, R2_b, R1_a, R1_b, R0;
  double outer_angle_multiplier;
  double inner_angle_multiplier;
  double inner_nonz_multiplier;
  double Clength;
  char input_section [50];
} d4est_geometry_cubed_sphere_vtk_attr_t;

/* This file was automatically generated.  Do not edit! */
void d4est_geometry_cubed_sphere_vtk_new(int mpirank,const char *input_file,const char *input_section,zlog_category_t *c_default,d4est_geometry_t *d4est_geom);

  
#endif
