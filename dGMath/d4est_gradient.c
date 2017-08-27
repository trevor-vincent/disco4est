#include <d4est_gradient.h>
#include <pXest.h>
#include <d4est_mesh.h>
#include <d4est_geometry.h>
#include <d4est_xyz_functions_ext.h>
#include <d4est_quadrature.h>
#include <d4est_quadrature_lobatto.h>
#include <d4est_util.h>
#include <sc_reduce.h>
#include <d4est_linalg.h>

void
d4est_gradient
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 double* u,
 int deg_lobatto,
 double* grad_u_quad [(P4EST_DIM)],
 int deg_quad,
 double* rst_xyz_quad [(P4EST_DIM)][(P4EST_DIM)]
)
{
  D4EST_ASSERT(object_type == QUAD_OBJECT_VOLUME);

  int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), deg_quad);
  int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM), deg_lobatto);

  double* du_di = P4EST_ALLOC(double, volume_nodes_lobatto);
  double* du_di_quad = P4EST_ALLOC(double, volume_nodes_quad);

  for (int j = 0; j < (P4EST_DIM); j++){
    for (int k = 0; k < volume_nodes_quad; k++){
      grad_u_quad[j][k] = 0.;
    }
  }
  
  for (int i = 0; i < (P4EST_DIM); i++){
    d4est_operators_apply_dij(d4est_ops, u, (P4EST_DIM), deg_lobatto, i, du_di);

    d4est_quadrature_interpolate(
                  d4est_ops,
                  d4est_quad,
                  d4est_geom,
                  object,
                  object_type,
                  integrand_type,
                  du_di,
                  deg_lobatto,
                  du_di_quad,
                  deg_quad
                 );

    
    for (int j = 0; j < (P4EST_DIM); j++){
      for (int k = 0; k < volume_nodes_quad; k++){
        grad_u_quad[j][k] += rst_xyz_quad[i][j][k]*du_di_quad[k];
      }
    }
  }

  P4EST_FREE(du_di);
  P4EST_FREE(du_di_quad);
}

double
d4est_gradient_l2_norm
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 double* u,
 int deg_lobatto,
 int deg_quad,
 double* rst_xyz_quad [(P4EST_DIM)][(P4EST_DIM)],
 double* jac_quad 
)
{
  double grad_l2_norm = 0.;
  int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), deg_quad);
  double* grad_u_quad[(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(grad_u_quad, volume_nodes_quad);

  d4est_gradient
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     object,
     object_type,
     integrand_type,
     u,
     deg_lobatto,
     grad_u_quad,
     deg_quad,
     rst_xyz_quad
    );
  
  for (int d = 0; d < (P4EST_DIM); d++){
    grad_l2_norm += d4est_quadrature_innerproduct
                    (
                     d4est_ops,
                     d4est_geom,
                     d4est_quad,
                     object,
                     object_type,
                     integrand_type, 
                     grad_u_quad[d],
                     grad_u_quad[d],
                     jac_quad,
                     deg_quad
                    );
  }
  
  D4EST_FREE_DIM_VEC(grad_u_quad);
  return grad_l2_norm;
}
