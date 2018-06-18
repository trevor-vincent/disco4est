#include <d4est_elliptic_eqns.h>

void
d4est_elliptic_eqns_apply_lhs
(
 p4est_t* p4est,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_elliptic_eqns_t* eqns,
 d4est_elliptic_data_t* vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors
)
{
  eqns->apply_lhs
    (
     p4est,
     ghost,
     ghost_data,
     vecs,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     eqns->user
    );
}


void
d4est_elliptic_eqns_build_residual
(
 p4est_t* p4est,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_elliptic_eqns_t* eqns,
 d4est_elliptic_data_t* vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors
)
{
  eqns->build_residual
    (
     p4est,
     ghost,
     ghost_data,
     vecs,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     eqns->user
    );
}
