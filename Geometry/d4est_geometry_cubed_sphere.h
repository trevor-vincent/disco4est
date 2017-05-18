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
  char input_section [50];
  
} d4est_geometry_cubed_sphere_attr_t;

void d4est_geometry_cubed_sphere_with_cube_hole_new(int mpirank,const char *input_file,const char *input_section,const char *printf_prefix,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_outer_shell_block_new(int mpirank,const char *input_file,const char *input_section,const char *printf_prefix,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_inner_shell_block_new(int mpirank,const char *input_file,const char *input_section,const char *printf_prefix,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_innerouter_shell_new(int mpirank,const char *input_file,const char *input_section,const char *printf_prefix,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_7tree_new(int mpirank,const char *input_file,const char *input_section,const char *printf_prefix,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_new(int mpirank,
                                     const char *input_file,
                                     const char *input_section,
                                     const char *printf_prefix,
                                     d4est_geometry_t *d4est_geom);
p4est_connectivity_t *d4est_connectivity_new_sphere_7tree(void);

void
d4est_geometry_cubed_sphere_outer_shell_block_get_custom_quadrature
(d4est_geometry_t* d4est_geom,
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [3],
 p4est_qcoord_t dq,
 int degree,
 long double* custom_abscissas,
 long double* custom_weights,
 int test_moments
);

#endif
