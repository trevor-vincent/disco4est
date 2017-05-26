#ifndef D4EST_QUADRATURE_GAUSS_LOBATTO_H
#define D4EST_QUADRATURE_GAUSS_LOBATTO_H 

/* This file was automatically generated.  Do not edit! */
double *d4est_quadrature_gauss_lobatto_get_interp_trans(d4est_operators_t *d4est_ops,d4est_quadrature_t *d4est_quad,d4est_geometry_t *d4est_geom,d4est_mesh_object_t mesh_object,int deg_lobatto,int deg_quad,int rst_direction);
double *d4est_quadrature_gauss_lobatto_get_interp(d4est_operators_t *d4est_ops,d4est_quadrature_t *d4est_quad,d4est_geometry_t *d4est_geom,d4est_mesh_object_t mesh_object,int deg_lobatto,int deg_quad,int rst_direction);
double *d4est_quadrature_gauss_lobatto_get_rst(d4est_operators_t *d4est_ops,d4est_quadrature_t *d4est_quad,d4est_geometry_t *d4est_geom,d4est_mesh_object_t mesh_object,int degree,int rst_direction);
double *d4est_quadrature_gauss_lobatto_get_weights(d4est_operators_t *d4est_ops,d4est_quadrature_t *d4est_quad,d4est_geometry_t *d4est_geom,d4est_mesh_object_t mesh_object,int degree,int rst_direction);
void d4est_quadrature_gauss_lobatto_new(d4est_quadrature_t *d4est_quad,d4est_geometry_t *d4est_geom,const char *input_file);

#endif
