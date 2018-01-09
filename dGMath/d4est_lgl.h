#ifndef D4EST_LGL_H
#define D4EST_LGL_H 

/* This file was automatically generated.  Do not edit! */
double d4est_lgl_interpolate(double rst[(P4EST_DIM)],double *lgl,double *f,int deg);
double d4est_lgl_lagrange_1d(double x,double *lgl,int j,int deg);
double d4est_lgl_gradjacobi(double r,double alpha,double beta,int N);
double d4est_lgl_jacobi(double r,double alpha,double beta,int N);
int d4est_lgl_get_nodes(int dim,int deg);

#endif
