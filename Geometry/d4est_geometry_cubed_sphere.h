#ifndef D4EST_GEOMETRY_SPHERE_H
#define D4EST_GEOMETRY_SPHERE_H

#include <pXest.h>
#include <p8est_geometry.h>


typedef struct
{
  double R2, R1, R0;
  double Clength;
  int compactify_outer_shell;
  int compactify_inner_shell;
    
} d4est_geometry_cubed_sphere_attr_t;

void d4est_geometry_cubed_sphere_with_cube_hole_new(int mpirank,const char *input_file,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_outer_shell_block_new(int mpirank,const char *input_file,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_inner_shell_block_new(int mpirank,const char *input_file,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_innerouter_shell_new(int mpirank,const char *input_file,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_7tree_new(int mpirank,const char *input_file,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_new(int mpirank,const char *input_file,d4est_geometry_t *d4est_geom);


#endif
