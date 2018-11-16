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
 double* drst_dxyz_quad [(P4EST_DIM)][(P4EST_DIM)],
 double* d2xyz_drstdrst_quad [(P4EST_DIM)][(P4EST_DIM)][P4EST_DIM]
)
{
  D4EST_ASSERT(object_type == QUAD_OBJECT_VOLUME);

  int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), deg_quad);
  int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM), deg_lobatto);

  double* du_db [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(du_db, volume_nodes_lobatto);
  double* du_db_quad [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(du_db_quad, volume_nodes_quad);
  double* d2u_dadb [P4EST_DIM][P4EST_DIM]; D4EST_ALLOC_DBYD_MAT(d2u_dadb, volume_nodes_lobatto);
  double* d2u_dadb_quad [P4EST_DIM][P4EST_DIM]; D4EST_ALLOC_DBYD_MAT(d2u_dadb_quad, volume_nodes_quad);

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
    d4est_operators_apply_dij(d4est_ops, u, (P4EST_DIM), deg_lobatto, b, du_db[b]);  
    d4est_quadrature_interpolate(
                                 d4est_ops,
                                 d4est_quad,
                                 d4est_geom,
                                 object,
                                 object_type,
                                 integrand_type,
                                 du_db[b],
                                 deg_lobatto,
                                 du_db_quad[b],
                                 deg_quad
    );

    for (int a = 0; a < (P4EST_DIM); a++){
      d4est_operators_apply_dij(d4est_ops, du_db[b], (P4EST_DIM), deg_lobatto, a, d2u_dadb[a][b]);
      d4est_quadrature_interpolate(
                                   d4est_ops,
                                   d4est_quad,
                                   d4est_geom,
                                   object,
                                   object_type,
                                   integrand_type,
                                   d2u_dadb[a][b],
                                   deg_lobatto,
                                   d2u_dadb_quad[a][b],
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
              for (int g = 0; g < (P4EST_DIM); g++){
                for (int l = 0; l < (P4EST_DIM); l++){        
                  hessian_u_quad[i][j][node] += -1.*drst_dxyz_quad[a][i][node]*drst_dxyz_quad[b][l][node]
                                                *drst_dxyz_quad[g][j][node]*d2xyz_drstdrst_quad[l][g][a][node]*
                                                du_db_quad[b][node];
                }
              }
            }
          }       
        }
      }
    }

    for (int i = 0; i < (P4EST_DIM); i++){
      for (int j = 0; j < (P4EST_DIM); j++){
        for (int node = 0; node < volume_nodes_quad; node++){
          for (int a = 0; a < (P4EST_DIM); a++){
            for (int b = 0; b < (P4EST_DIM); b++){
              hessian_u_quad[i][j][node] += drst_dxyz_quad[a][i][node]*drst_dxyz_quad[b][j][node]
                                            *d2u_dadb_quad[b][a][node];
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
            for (int b = 0; b < (P4EST_DIM); b++){
              for (int g = 0; g < (P4EST_DIM); g++){
                for (int l = 0; l < (P4EST_DIM); l++){        
                  hessian_trace_u_quad[node] += -1.*drst_dxyz_quad[a][i][node]*drst_dxyz_quad[b][l][node]
                                                *drst_dxyz_quad[g][i][node]*d2xyz_drstdrst_quad[l][g][a][node]*
                                                du_db_quad[b][node];

                }
              }
            }
          }       
        }
    }


        for (int node = 0; node < volume_nodes_quad; node++){
          for (int i = 0; i < (P4EST_DIM); i++){
          for (int a = 0; a < (P4EST_DIM); a++){
            for (int b = 0; b < (P4EST_DIM); b++){
              hessian_trace_u_quad[node] += drst_dxyz_quad[a][i][node]*drst_dxyz_quad[b][i][node]
                                            *d2u_dadb_quad[b][a][node];
            }
          }
        }
    }
      
  }
  
  D4EST_FREE_DIM_VEC(du_db);
  D4EST_FREE_DIM_VEC(du_db_quad);
  D4EST_FREE_DBYD_MAT(d2u_dadb);
  D4EST_FREE_DBYD_MAT(d2u_dadb_quad);
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
 
  
  for (int i = 0; i < p4est->local_num_quadrants; i++){

    d4est_element_data_t* ed = d4est_factors->element_data[i];
    
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

  for (int d1 = 0; d1 < (P4EST_DIM); d1++){
    for (int d2 = 0; d2 < (P4EST_DIM); d2++){
      for (int d3 = 0; d3 < (P4EST_DIM); d3++){
        P4EST_FREE(d2xdrdr[d1][d2][d3]);
      }
    }
  } 
}
