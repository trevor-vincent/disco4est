#include <d4est_estimator_residual.h>
#include <d4est_mesh.h>

void
d4est_estimator_residual_compute
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 diam_compute_option_t diam_opt
)
{
  d4est_elliptic_eqns_build_residual
    (
     p4est,
     ghost,
     ghost_data,
     fcns,
     vecs,
     d4est_ops,
     d4est_geom,
     d4est_quad
    );
  
  d4est_mesh_compute_l2_norm_sqr
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     vecs->Au,
     vecs->local_nodes,
     STORE_LOCALLY,
     NULL,
     NULL
    );

}
