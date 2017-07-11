#include <d4est_elliptic_eqns.h>

void
d4est_elliptic_eqns_apply_lhs
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_elliptic_eqns_t* eqns,
 d4est_elliptic_data_t* vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad
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
     eqns->user
    );
}


void
d4est_elliptic_eqns_build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_elliptic_eqns_t* eqns,
 d4est_elliptic_data_t* vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad
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
     eqns->user
    );
}
