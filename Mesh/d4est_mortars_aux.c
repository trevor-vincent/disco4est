#include <pXest.h>
#include <d4est_mortars_aux.h>
#include <d4est_quadrature.h>
#include <d4est_quadrature_lobatto.h>
#include <d4est_geometry.h>
#include <d4est_operators.h>

void
d4est_mortars_compute_dxyz_drst_face_isoparametric
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int which_tree,
 int face,
 int deg,
 int deg_quad,
 double* xyz_rst_face_quad [(P4EST_DIM)][(P4EST_DIM)]
)
{
  int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), deg);
  
  double* xyz [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(xyz, volume_nodes);
  
  d4est_rst_t rst_points_lobatto;
  rst_points_lobatto.r = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, NULL, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, deg, 0);
  rst_points_lobatto.s = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, NULL, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, deg, 1);
  rst_points_lobatto.t = NULL;
#if (P4EST_DIM)==3
  rst_points_lobatto.t = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, NULL, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, deg, 2);
#endif

  d4est_geometry_compute_xyz
    (
     d4est_ops,
     d4est_geom,
     rst_points_lobatto,
     which_tree,
     deg,
     q0,
     dq,
     xyz
    );

  
  double* tmp = P4EST_ALLOC(double, volume_nodes);
  double* tmp2 = P4EST_ALLOC(double, volume_nodes);
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int d1 = 0; d1 < (P4EST_DIM); d1++){
      d4est_operators_apply_dij(d4est_ops, &xyz[d][0], (P4EST_DIM), deg, d1, tmp);
      d4est_operators_apply_slicer(d4est_ops,
                                   tmp,
                                   (P4EST_DIM),
                                   face,
                                   deg,
                                   tmp2);     
      d4est_quadrature_interpolate(d4est_ops, d4est_quad, d4est_geom, NULL, QUAD_OBJECT_MORTAR, QUAD_INTEGRAND_UNKNOWN, tmp2, deg, &xyz_rst_face_quad[d][d1][0], deg_quad);
    }
  }
  P4EST_FREE(tmp);
  P4EST_FREE(tmp2);
  D4EST_FREE_DIM_VEC(xyz);
}
/* This file was automatically generated.  Do not edit! */
