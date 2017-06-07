#include <d4est_operators.h>
#include <d4est_quadrature.h>
#include <d4est_quadrature_legendre.h>
#include <util.h>

void
d4est_quadrature_legendre_new
(
 d4est_quadrature_t* d4est_quad,
 d4est_geometry_t* d4est_geom,
 const char* input_file
)
{
  d4est_quad->get_weights = d4est_quadrature_legendre_get_weights;
  d4est_quad->get_rst = d4est_quadrature_legendre_get_rst;
  d4est_quad->get_interp = d4est_quadrature_legendre_get_interp;
  d4est_quad->get_interp_trans = d4est_quadrature_legendre_get_interp_trans;
  d4est_quad->user_destroy = NULL;
  d4est_quad->user_reinit = NULL;
  d4est_quad->user = NULL;
}

double*
d4est_quadrature_legendre_get_weights
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 int degree,
 int rst_direction
){
  return d4est_operators_fetch_GL_weights_1d(d4est_ops, degree);
}

double*
d4est_quadrature_legendre_get_rst
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 int degree,
 int rst_direction
){
  if (object_type == D4EST_QUAD_FACE)
    return d4est_operators_fetch_Gauss_xyz_nd(d4est_ops, (P4EST_DIM)-1, degree, rst_direction);
  else if (object_type == D4EST_QUAD_VOLUME)
    return d4est_operators_fetch_Gauss_xyz_nd(d4est_ops, (P4EST_DIM), degree, rst_direction);
  else{
    mpi_abort("[D4EST_ERROR]: Object type must be QUAD_FACE or QUAD_VOLUME");
    return NULL;
  }
}


double*
d4est_quadrature_legendre_get_interp
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 int deg_lobatto,
 int deg_quad,
 int rst_direction
){
  return d4est_operators_fetch_ref_GLL_to_GL_interp_1d(d4est_ops, deg_lobatto, deg_quad);
}

double*
d4est_quadrature_legendre_get_interp_trans
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 int deg_lobatto,
 int deg_quad,
 int rst_direction
 
){
  return d4est_operators_fetch_ref_GLL_to_GL_interp_trans_1d(d4est_ops, deg_lobatto, deg_quad);
}

