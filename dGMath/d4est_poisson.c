#include <pXest.h>
#include <d4est_element_data.h>
#include <d4est_elliptic_data.h>
#include <d4est_operators.h>
#include <d4est_linalg.h>
#include <d4est_poisson.h>
#include <d4est_poisson_flux.h>
#include <d4est_xyz_functions.h>
#include <d4est_quadrature.h>
#include <d4est_mortars.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_mesh.h>
#include <d4est_util.h>

void
 d4est_poisson_build_rhs_with_strong_bc
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_elliptic_data_t* prob_vecs,
 d4est_poisson_flux_data_t* flux_fcn_data_for_build_rhs,
 double* rhs,
 d4est_xyz_fcn_t problem_rhs_fcn,
 d4est_mesh_init_field_option_t init_option,
 void* ctx
)
{
  int local_nodes = prob_vecs->local_nodes;

  double* u_eq_0 = P4EST_ALLOC_ZERO(double, local_nodes);
  double* Au_eq_0 = P4EST_ALLOC_ZERO(double, local_nodes);

  d4est_elliptic_data_t elliptic_data_for_rhs;
  d4est_elliptic_data_copy_ptrs(prob_vecs, &elliptic_data_for_rhs);
  elliptic_data_for_rhs.u = u_eq_0; 
  elliptic_data_for_rhs.Au = Au_eq_0; 

  d4est_poisson_apply_aij(p4est, ghost, ghost_data, &elliptic_data_for_rhs, flux_fcn_data_for_build_rhs, d4est_ops, d4est_geom, d4est_quad, d4est_factors); 


  double* f = NULL;
  if (init_option == INIT_FIELD_ON_LOBATTO){
    f = P4EST_ALLOC(double, local_nodes);
  }
  else if (init_option == INIT_FIELD_ON_QUAD){
    int local_quad_nodes = d4est_mesh_get_local_quad_nodes(p4est);
    f = P4EST_ALLOC(double, local_quad_nodes);
  }
  else {
    D4EST_ABORT("Not a support init option");
  }
  
  d4est_mesh_init_field
    (
     p4est,
     f,
     problem_rhs_fcn,
     d4est_ops,
     d4est_geom,
     init_option,
     ctx
    );
  
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data;

        d4est_quadrature_volume_t mesh_object = {.dq = ed->dq,
                                                 .tree = ed->tree,
                                                 .q[0] = ed->q[0],
                                                 .q[1] = ed->q[1],
#if (P4EST_DIM)==3                                              
                                                 .q[2] = ed->q[2],
#endif
                                                 .element_id = ed->id
                                                };

        if (init_option == INIT_FIELD_ON_LOBATTO){
          d4est_quadrature_apply_mass_matrix
            (
             d4est_ops,
             d4est_geom,
             d4est_quad,
             &mesh_object,
             QUAD_OBJECT_VOLUME,
             QUAD_INTEGRAND_UNKNOWN,
             &f[ed->nodal_stride],
             ed->deg,
             ed->J_quad,
             ed->deg_quad,
             &rhs[ed->nodal_stride]
            );
        }
        else if (init_option == INIT_FIELD_ON_QUAD){
          d4est_quadrature_apply_galerkin_integral
            (
             d4est_ops,
             d4est_geom,
             d4est_quad,
             &mesh_object,
             QUAD_OBJECT_VOLUME,
             QUAD_INTEGRAND_UNKNOWN,
             &f[ed->quad_stride],
             ed->deg,
             ed->J_quad,
             ed->deg_quad,
             &rhs[ed->nodal_stride]
            );
        }
        else {
          D4EST_ABORT("Not a supported init option");
        }
      }
    }

  d4est_linalg_vec_axpy(-1., Au_eq_0, rhs, local_nodes);  
  P4EST_FREE(u_eq_0);
  P4EST_FREE(Au_eq_0);
  P4EST_FREE(f);
}

void
d4est_poisson_apply_stiffness_matrix
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad
)
{
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data; 

        d4est_quadrature_volume_t mesh_vol = {.dq = ed->dq,
                                              .tree = ed->tree,
                                              .q[0] = ed->q[0],
                                              .q[1] = ed->q[1],
#if (P4EST_DIM)==3                                              
                                              .q[2] = ed->q[2],
#endif
                                              .element_id = ed->id
                                             };

        d4est_quadrature_apply_stiffness_matrix
          (
           d4est_ops,
           d4est_quad,
           d4est_geom,
           &mesh_vol,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &ed->u_elem[0],
           ed->deg,
           ed->J_quad,
           ed->rst_xyz_quad,
           ed->deg_quad,
           ed->Au_elem
          );
      }

    }
}

void
d4est_poisson_compute_dudr
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 double* dudr [(P4EST_DIM)]
){

 int stride = 0;
 for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data;         
        int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM),ed->deg);
        for (int i = 0; i < (P4EST_DIM); i++){
          ed->dudr_elem[i] = &dudr[i][stride];
          d4est_operators_apply_dij(d4est_ops, &ed->u_elem[0], (P4EST_DIM), ed->deg, i, &ed->dudr_elem[i][0]);
        }
        stride += volume_nodes_lobatto;
      }
    }


  for (int gid = 0; gid < ghost->ghosts.elem_count; gid++){
    d4est_element_data_t* ged = &ghost_data[gid];
    int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM),ged->deg);
    for (int i = 0; i < (P4EST_DIM); i++){
      ged->dudr_elem[i] = &dudr[i][stride];
      d4est_operators_apply_dij(d4est_ops, &ged->u_elem[0], (P4EST_DIM), ged->deg, i, &ged->dudr_elem[i][0]);
    }
    stride += volume_nodes_lobatto;
  }

 
}


void
d4est_poisson_apply_mortar_matrices
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_poisson_flux_data_t* flux_fcn_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors
)
{
  d4est_mortars_fcn_ptrs_t flux_fcns = d4est_poisson_flux_fetch_fcns(flux_fcn_data);

  
  d4est_mortars_compute_flux_on_local_elements
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &flux_fcns,
     DO_NOT_EXCHANGE_GHOST_DATA /* already done above */
    );
}




void
d4est_poisson_apply_aij
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_poisson_flux_data_t* flux_fcn_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors
)
{
  d4est_poisson_flux_init_element_data
    (
     p4est,
     d4est_ops,
     prob_vecs->u,
     prob_vecs->Au
    );
  
  d4est_poisson_apply_stiffness_matrix
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad
    );


  p4est_ghost_exchange_data(p4est,ghost,ghost_data);

  int ghost_nodes = d4est_mesh_get_ghost_nodes(ghost, ghost_data);
  double* dudr [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(dudr,
                                                  prob_vecs->local_nodes
                                                  + ghost_nodes
                                                 );
  
  d4est_poisson_compute_dudr
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     dudr
    );

  d4est_poisson_apply_mortar_matrices
    (
     p4est,
     ghost,
     ghost_data,
     flux_fcn_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors
    );

  D4EST_FREE_DIM_VEC(dudr);
  
  
  /* printf("after apply mortar matrices "); */
  /* DEBUG_PRINT_ARR_DBL(prob_vecs->Au, prob_vecs->local_nodes); */
  
}
