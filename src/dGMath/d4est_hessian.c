#include <d4est_hessian.h>
#include <pXest.h>
#include <d4est_mesh.h>
#include <d4est_util.h>
#include <d4est_geometry.h>
#include <d4est_xyz_functions_ext.h>
#include <d4est_quadrature.h>
#include <d4est_quadrature_lobatto.h>
#include <d4est_util.h>
#include <sc_reduce.h>
#include <d4est_linalg.h>

void
d4est_hessian
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 double* u,
 int deg_lobatto,
 double* hessian_u_quad [(P4EST_DIM)][(P4EST_DIM)],
 double* hessian_trace_u_quad,
 int deg_quad,
 double* drdx_quad [(P4EST_DIM)][(P4EST_DIM)],
 double* d2xdrdr_quad [(P4EST_DIM)][(P4EST_DIM)][P4EST_DIM]
)
{
  D4EST_ASSERT(object_type == QUAD_OBJECT_VOLUME);

  int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), deg_quad);
  int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM), deg_lobatto);

  double* du [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(du, volume_nodes_lobatto);
  double* du_quad [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(du_quad, volume_nodes_quad);
  double* d2u [P4EST_DIM][P4EST_DIM]; D4EST_ALLOC_DBYD_MAT(d2u, volume_nodes_lobatto);
  double* d2u_quad [P4EST_DIM][P4EST_DIM]; D4EST_ALLOC_DBYD_MAT(d2u_quad, volume_nodes_quad);

  double* d2rdrdx_quad [(P4EST_DIM)][(P4EST_DIM)][(P4EST_DIM)];

  for (int m = 0; m < (P4EST_DIM); m++){
    for (int n = 0; n < (P4EST_DIM); n++){
      for (int k = 0; k < (P4EST_DIM); k++){
        d2rdrdx_quad[m][n][k] = P4EST_ALLOC_ZERO(double, volume_nodes_quad);
        for (int node = 0; node < volume_nodes_quad; node++){
          for (int a = 0; a < (P4EST_DIM); a++) {
            for (int l = 0; l < (P4EST_DIM); ++l) {
              d2rdrdx_quad[m][n][k][node] -= drdx_quad[m][l][node]*drdx_quad[a][n][node]*d2xdrdr_quad[l][a][k][node];
              /* printf("d2rdrdx_quad[%d][%d][%d] = %.15f\n", m,n,k, d2rdrdx_quad[m][n][k][node]); */
              /* printf("drdx_quad[%d][%d] = %.15f\n", m,n,k, drdx_quad[m][n][node]); */
            }
          }
        }        
      }
    }
  }

  
  if (hessian_u_quad != NULL && hessian_u_quad[0][0] != NULL){
      for (int j = 0; j < (P4EST_DIM); j++){
        for (int k = 0; k < (P4EST_DIM); k++){
          for (int l = 0; l < volume_nodes_quad; l++){
          hessian_u_quad[j][k][l] = 0.;
        }
      }
    }
  }
  
  if (hessian_trace_u_quad != NULL){
    for (int l = 0; l < volume_nodes_quad; l++){
      hessian_trace_u_quad[l] = 0.;
    }
  }
  
  for (int b = 0; b < (P4EST_DIM); b++){
    d4est_operators_apply_dij(d4est_ops, u, (P4EST_DIM), deg_lobatto, b, du[b]);

    
    /* DEBUG_PRINT_ARR_DBL(du[b], volume_nodes_lobatto); */
    d4est_quadrature_interpolate(
                                 d4est_ops,
                                 d4est_quad,
                                 d4est_geom,
                                 object,
                                 object_type,
                                 integrand_type,
                                 du[b],
                                 deg_lobatto,
                                 du_quad[b],
                                 deg_quad
    );

    for (int a = 0; a < (P4EST_DIM); a++){
      d4est_operators_apply_dij(d4est_ops, du[b], (P4EST_DIM), deg_lobatto, a, d2u[b][a]);
      d4est_quadrature_interpolate(
                                   d4est_ops,
                                   d4est_quad,
                                   d4est_geom,
                                   object,
                                   object_type,
                                   integrand_type,
                                   d2u[b][a],
                                   deg_lobatto,
                                   d2u_quad[b][a],
                                   deg_quad
      );
    }    
  }

  if (hessian_u_quad != NULL && hessian_u_quad[0][0] != NULL){
    for (int i = 0; i < (P4EST_DIM); i++){
      for (int j = 0; j < (P4EST_DIM); j++){
        for (int node = 0; node < volume_nodes_quad; node++){
          for (int a = 0; a < (P4EST_DIM); a++){
            for (int b = 0; b < (P4EST_DIM); b++){
              hessian_u_quad[j][i][node] += drdx_quad[a][i][node]*(d2rdrdx_quad[b][j][a][node]*du_quad[b][node] + drdx_quad[b][j][node]*d2u_quad[b][a][node]);
            }
          }       
        }
      }
    }

    
  }
  if (hessian_trace_u_quad != NULL){
    for (int node = 0; node < volume_nodes_quad; node++){
      for (int i = 0; i < (P4EST_DIM); i++){
        for (int a = 0; a < (P4EST_DIM); a++){
          /* printf("d2u_quad[%d][%d][node] = %.15f\n", i,a, d2u_quad[i][a][node]); */
          for (int b = 0; b < (P4EST_DIM); b++){
            hessian_trace_u_quad[node] += drdx_quad[a][i][node]*(d2rdrdx_quad[b][i][a][node]*du_quad[b][node] + drdx_quad[b][i][node]*d2u_quad[b][a][node]);
          }
        }       
      }
    }      
  }

  for (int m = 0; m < (P4EST_DIM); m++){
    for (int n = 0; n < (P4EST_DIM); n++){
      for (int k = 0; k < (P4EST_DIM); k++){
        P4EST_FREE(d2rdrdx_quad[m][n][k]);
      }
    }
  }
  
  D4EST_FREE_DIM_VEC(du);
  D4EST_FREE_DIM_VEC(du_quad);
  D4EST_FREE_DBYD_MAT(d2u);
  D4EST_FREE_DBYD_MAT(d2u_quad);
}

void
d4est_hessian_compute_d2xdrdr_on_quadrature_points_of_element
(
 d4est_element_data_t* ed,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_hessian_compute_method_t compute_method,
 double* d2xdrdr_quad [(P4EST_DIM)][(P4EST_DIM)][(P4EST_DIM)]
)
{
  int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), ed->deg_quad);
  int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);

  d4est_quadrature_volume_t mesh_object;
  mesh_object.dq =  ed->dq;
  mesh_object.tree = ed->tree;
  mesh_object.element_id = ed->id;
        
  mesh_object.q[0] = ed->q[0];
  mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
  mesh_object.q[2] = ed->q[2];
#endif

  if (compute_method == HESSIAN_ANALYTICAL){

    d4est_rst_t rst_points_quad;
    rst_points_quad = d4est_quadrature_get_rst_points
                      (
                       d4est_ops,
                       d4est_quad,
                       d4est_geom,
                       &mesh_object,
                       QUAD_OBJECT_VOLUME,
                       QUAD_INTEGRAND_UNKNOWN,
                       ed->deg_quad
                      );

    double d2xdrdr_quad_temp [(P4EST_DIM)][(P4EST_DIM)][(P4EST_DIM)];
    for (int node = 0; node < volume_nodes_quad; node++){

      if (d4est_geom->D2X != NULL){
        double rst [(P4EST_DIM)];
        rst[0] = rst_points_quad.r[node];
        rst[1] = rst_points_quad.s[node];
#if (P4EST_DIM)==3
        rst[2] = rst_points_quad.t[node];
#endif
        d4est_geom->D2X(d4est_geom,
                        ed->tree,
                        ed->q,
                        ed->dq,
                        rst,
                        d2xdrdr_quad_temp
                       );

        for (int d1 = 0; d1 < (P4EST_DIM); d1++){
          for (int d2 = 0; d2 < (P4EST_DIM); d2++){
            for (int d3 = 0; d3 < (P4EST_DIM); d3++){
              d2xdrdr_quad[d1][d2][d3][node]
                = d2xdrdr_quad_temp[d1][d2][d3];
            }
          }
        }    
      }
      else {
        D4EST_ABORT("Compute method should be numerical if D2X=NULL");
      }

    }
  }
  else if (compute_method == HESSIAN_NUMERICAL){

    d4est_mesh_data_on_element_t md_on_e = d4est_mesh_data_on_element
                                           (
                                            d4est_factors,
                                            ed
                                           );


    double* dr_tmp = P4EST_ALLOC(double, volume_nodes);
    double* drdr_tmp = P4EST_ALLOC(double, volume_nodes);
    
    for (int d1 = 0; d1 < (P4EST_DIM); d1++){
      for (int d2 = 0; d2 < (P4EST_DIM); d2++){
        for (int d3 = 0; d3 < (P4EST_DIM); d3++){
          d4est_operators_apply_dij(d4est_ops, md_on_e.xyz[d1], (P4EST_DIM), ed->deg, d2, dr_tmp);  
          d4est_operators_apply_dij(d4est_ops, dr_tmp, (P4EST_DIM), ed->deg, d3, drdr_tmp);                
          d4est_quadrature_interpolate(
                                       d4est_ops,
                                       d4est_quad,
                                       d4est_geom,
                                       &mesh_object,
                                       QUAD_OBJECT_VOLUME,
                                       QUAD_INTEGRAND_UNKNOWN,
                                       drdr_tmp,
                                       ed->deg,
                                       d2xdrdr_quad[d1][d2][d3],
                                       ed->deg_quad
          );
        }
      }
    }

    P4EST_FREE(dr_tmp);
    P4EST_FREE(drdr_tmp);
  }
  else {
    D4EST_ABORT("compute method must be NUMERICAL or ANALYTICAL");
  }
}


void
d4est_hessian_compute_hessian_trace_of_field_on_quadrature_points
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_hessian_compute_method_t compute_method,
 double* field_lobatto, /* field of size local_nodes AKA local LGL nodes */
 double* del2field /* field of size local_quad_nodes AKA local quadrature nodes */
)
{
  double* d2xdrdr [(P4EST_DIM)][(P4EST_DIM)][(P4EST_DIM)];    
  for (int d1 = 0; d1 < (P4EST_DIM); d1++){
    for (int d2 = 0; d2 < (P4EST_DIM); d2++){
      for (int d3 = 0; d3 < (P4EST_DIM); d3++){
        d2xdrdr[d1][d2][d3] = P4EST_ALLOC(double, d4est_factors->local_sizes.local_nodes_quad);
      }
    }
  }
 
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;

      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);
        
        double* d2xyz_drstdrst_quad_on_e [(P4EST_DIM)][(P4EST_DIM)][(P4EST_DIM)];
        for (int d1 = 0; d1 < (P4EST_DIM); d1++){
          for (int d2 = 0; d2 < (P4EST_DIM); d2++){
            for (int d3 = 0; d3 < (P4EST_DIM); d3++){
              d2xyz_drstdrst_quad_on_e[d1][d2][d3] = &d2xdrdr[d1][d2][d3][ed->quad_stride];
            }
          }
        }

        d4est_hessian_compute_d2xdrdr_on_quadrature_points_of_element
          (
           ed,
           d4est_ops,
           d4est_geom,
           d4est_quad,
           d4est_factors,
           compute_method,
           d2xyz_drstdrst_quad_on_e
          );

        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq =  ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = ed->q[2];
#endif

    
        d4est_mesh_data_on_element_t md_on_e = d4est_mesh_data_on_element
                                               (
                                                d4est_factors,
                                                ed
                                               );

    
        d4est_hessian
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           &mesh_object,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &field_lobatto[ed->nodal_stride],
           ed->deg,
           NULL,
           &del2field[ed->quad_stride],
           ed->deg_quad,
           md_on_e.rst_xyz_quad,
           d2xyz_drstdrst_quad_on_e
          );
      }
    }

  for (int d1 = 0; d1 < (P4EST_DIM); d1++){
    for (int d2 = 0; d2 < (P4EST_DIM); d2++){
      for (int d3 = 0; d3 < (P4EST_DIM); d3++){
        P4EST_FREE(d2xdrdr[d1][d2][d3]);
      }
    }
  } 
}
