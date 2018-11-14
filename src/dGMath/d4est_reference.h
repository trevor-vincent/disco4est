#ifndef D4EST_REFERENCE_H
#define D4EST_REFERENCE_H 

#include <assert.h>
#include <d4est_util.h>
#include <pXest.h>

typedef struct {

  double* r;
  double* s;
  double* t;
  
} d4est_rst_t;

/* This file was automatically generated.  Do not edit! */
int d4est_reference_get_mirrored_face(int face);
void d4est_reference_dir_and_side_of_face(int face,int *dir,int *side);
void d4est_reference_rtox_array(double *r,double xl,double h,double *x,int nodes);
double d4est_reference_rtox(double r,double xl,double h);
int d4est_reference_corner_to_node(int dim,int deg,int corner);
int d4est_reference_reorient_face_order(int face_dim,int f_m,int f_p,int o,int i);
void d4est_reference_get_normal(int face,int dim,double *n);
void d4est_reference_child_to_parent_coords(int c,double *r);
int d4est_reference_is_child_left_or_right(int c,int dir);

#endif
