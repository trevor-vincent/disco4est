#ifndef D4EST_GEOMETRY_SPHERE_H
#define D4EST_GEOMETRY_SPHERE_H

#include <pXest.h>
#include <p8est_geometry.h>
#include <zlog.h>


typedef struct
{
  double R2, R1, R0;
  double Clength;
  int compactify_outer_shell;
  int compactify_inner_shell;
  char input_section [50];
  
} d4est_geometry_cubed_sphere_attr_t;

/* This file was automatically generated.  Do not edit! */
void d4est_geometry_cubed_sphere_with_sphere_hole_new(int mpirank,const char *input_file,const char *input_section,zlog_category_t *c_default,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_with_cube_hole_new(int mpirank,const char *input_file,const char *input_section,zlog_category_t *c_default,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_outer_shell_block_new(int mpirank,const char *input_file,const char *input_section,zlog_category_t *c_default,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_outer_shell_block_new_aux(d4est_geometry_t *d4est_geom,d4est_geometry_cubed_sphere_attr_t *sphere_attrs);
void d4est_geometry_cubed_sphere_inner_shell_block_new(int mpirank,const char *input_file,const char *input_section,zlog_category_t *c_default,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_innerouter_shell_new(int mpirank,const char *input_file,const char *input_section,zlog_category_t *c_default,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_7tree_new(int mpirank,const char *input_file,const char *input_section,zlog_category_t *c_default,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_new(int mpirank,const char *input_file,const char *input_section,zlog_category_t *c_default,d4est_geometry_t *d4est_geom);
void d4est_geometry_cubed_sphere_X_aux_rotate(int which_tree,double xyz_top[(P4EST_DIM)],double xyz[(P4EST_DIM)]);
void d4est_geometry_cubed_sphere_DX_aux_rotate(int which_tree,double dxyz_drst_top[(P4EST_DIM)][(P4EST_DIM)],double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]);


#endif
