#ifndef D4EST_QUADRATURE_COMPACTIFIED_H
#define D4EST_QUADRATURE_COMPACTIFIED_H 

#include <d4est_geometry.h>

/* This file was automatically generated.  Do not edit! */
void d4est_quadrature_compactified_c1tpc2_neg3_aa_and_bb(int n,long double *aa,long double *bb,void *user);
void d4est_quadrature_compactified_c1tpc2_neg2_aa_and_bb(int n,long double *aa,long double *bb,void *user);
void d4est_quadrature_compactified_c1tpc2_neg1_aa_and_bb(int n,long double *aa,long double *bb,void *user);
long double d4est_quadrature_compactified_c1tpc2_neg1_moment_fcn(int n,void *user);
long double d4est_quadrature_compactified_c1tpc2_neg4_moment_fcn(int n,void *user);
long double d4est_quadrature_compactified_c1tpc2_neg2_moment_fcn(int n,void *user);
long double d4est_quadrature_compactified_c1tpc2_neg3_moment_fcn(int n,void *user);
long double d4est_quadrature_compactified_c1tpc2_neg1_weight_fcn(long double x,void *user);
long double d4est_quadrature_compactified_c1tpc2_neg2_weight_fcn(long double x,void *user);
long double d4est_quadrature_compactified_c1tpc2_neg3_weight_fcn(long double x,void *user);
long double d4est_quadrature_compactified_c1tpc2_neg4_weight_fcn(long double x,void *user);
void d4est_quadrature_compactified_c1tpc2_neg4_aa_and_bb(int n,long double *aa,long double *bb,void *user);
void d4est_quadrature_compactified_setup_storage(p4est_t *p4est,p4est_ghost_t *ghost,void *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad);
void d4est_quadrature_compactified_new(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,const char *input_file,const char *input_section);
void d4est_quadrature_compactified_storage_destroy(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad);
void d4est_quadrature_compactified_compute_rst_face(double *abscissas,int deg_quad,double *rst_face[(P4EST_DIM)-1]);
double *d4est_quadrature_compactified_get_rst(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,void *object,d4est_quadrature_object_type_t object_type,d4est_quadrature_integrand_type_t integrand_type,int degree,int rst_direction);
int d4est_quadrature_compactified_check_type(d4est_quadrature_t *d4est_quad);
void d4est_quadrature_compactified_compute_interpolation_matrix_and_transpose(d4est_operators_t *d4est_ops,double *abscissas,double *interp,double *interp_trans,int deg,int deg_quad);
void d4est_quadrature_compactified_compute_rst_volume(double *abscissas,int deg_quad,double *rst_volume,int rst_direction);
void d4est_quadrature_compactified_compute_abscissas_and_weights(d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,double *abscissas,double *weights,p4est_topidx_t tree,p4est_qcoord_t q0[(P4EST_DIM)],p4est_qcoord_t dq,int degree);

#endif
