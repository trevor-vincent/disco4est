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
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_elliptic_data_t* prob_vecs,
 d4est_poisson_flux_data_t* flux_fcn_data_for_build_rhs,
 double * restrict  rhs,
 d4est_xyz_fcn_t problem_rhs_fcn,
 d4est_mesh_init_field_option_t init_option,
 void* ctx,
 int which_field
)
{
  int local_nodes = prob_vecs->local_nodes;

  double* u_eq_0 = P4EST_ALLOC_ZERO(double, local_nodes);
  double* Au_eq_0 = P4EST_ALLOC_ZERO(double, local_nodes);

  d4est_elliptic_data_t elliptic_data_for_rhs;
  d4est_elliptic_data_copy_ptrs(prob_vecs, &elliptic_data_for_rhs);
  elliptic_data_for_rhs.u = u_eq_0; 
  elliptic_data_for_rhs.Au = Au_eq_0; 

  
  d4est_poisson_apply_aij(p4est, ghost, ghost_data, &elliptic_data_for_rhs, flux_fcn_data_for_build_rhs, d4est_ops, d4est_geom, d4est_quad, d4est_factors, which_field); 

  



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
     d4est_factors,
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



        double* J_quad = d4est_mesh_get_jacobian_on_quadrature_points(d4est_factors,
                                                                      ed);


        
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
             J_quad,
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
             J_quad,
             ed->deg_quad,
             &rhs[ed->nodal_stride]
            );
        }
        else {
          D4EST_ABORT("Not a supported init option");
        }
      }
    }

  /* DEBUG_PRINT_ARR_DBL_SUM(Au_eq_0, local_nodes); */
  /* DEBUG_PRINT_ARR_DBL_SUM(rhs, local_nodes); */

  /* for (int i = 0; i < local_nodes; i++){ */
    /* Au_eq_0[i] = 0.; */
  /* } */
  
  d4est_linalg_vec_axpy(-1., Au_eq_0, rhs, local_nodes);
  /* DEBUG_PRINT_ARR_DBL_SUM(rhs, local_nodes); */
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
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 double * restrict  u,
 double * restrict  Au,
 int local_nodes,
 int which_field
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


        double* J_quad = d4est_mesh_get_jacobian_on_quadrature_points(d4est_factors,
                                                                      ed);



        d4est_mesh_data_on_element_t md_on_e = d4est_mesh_data_on_element
                                               (
                                                d4est_factors,
                                                ed
                                               );
        
        d4est_quadrature_apply_stiffness_matrix
          (
           d4est_ops,
           d4est_quad,
           d4est_geom,
           &mesh_vol,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &u[which_field*local_nodes + ed->nodal_stride],           
           ed->deg,
           J_quad,
           md_on_e.rst_xyz_quad,
           ed->deg_quad,
           &Au[which_field*local_nodes + ed->nodal_stride]
          );
      }

    }
}

void
d4est_poisson_compute_dudr
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_ghost_data_t* d4est_ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 double* dudr_local [(P4EST_DIM)],
 double* dudr_ghost [(P4EST_DIM)],
 double * restrict  u,
 int local_nodes,
 int which_field
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
          d4est_operators_apply_dij(d4est_ops, &u[which_field*local_nodes + ed->nodal_stride], (P4EST_DIM), ed->deg, i, &dudr_local[i][stride]);
        }
        stride += volume_nodes_lobatto;
      }
    }

 stride = 0;
 for (int gid = 0; gid < d4est_ghost->ghost->ghosts.elem_count; gid++){
   d4est_element_data_t* ged = &d4est_ghost->ghost_elements[gid];
   int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM),ged->deg);
   for (int i = 0; i < (P4EST_DIM); i++){
     double* u_ghost_elem = d4est_ghost_data_get_field_on_element(ged,0,d4est_ghost_data);
     d4est_operators_apply_dij(d4est_ops, u_ghost_elem, (P4EST_DIM), ged->deg, i, &dudr_ghost[i][stride]);
   }
   stride += volume_nodes_lobatto;
 }
}


void
d4est_poisson_apply_mortar_matrices
(
 p4est_t* p4est,
 d4est_ghost_t* ghost,
 /* d4est_ghost_data_t* ghost_data, */
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
     /* ghost_data, */
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &flux_fcns
     /* DO_NOT_EXCHANGE_GHOST_DATA /\* already done above *\/ */
    );
}




void
d4est_poisson_apply_aij
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_ghost_data_t* d4est_ghost_data,
 d4est_elliptic_data_t* d4est_elliptic_data,
 d4est_poisson_flux_data_t* flux_fcn_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 int which_field
)
{
  /* d4est_poisson_flux_init_element_data */
  /*   ( */
  /*    p4est, */
  /*    d4est_ops, */
  /*    d4est_elliptic_data->u, */
  /*    d4est_elliptic_data->Au */
  /*   ); */
  
  d4est_poisson_apply_stiffness_matrix
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     d4est_elliptic_data->u,
     d4est_elliptic_data->Au,
     d4est_elliptic_data->local_nodes,
     which_field
    );


  
  d4est_ghost_data_exchange(p4est,d4est_ghost,d4est_ghost_data,d4est_elliptic_data->u);
  
  /* p4est_ghost_exchange_data(p4est,ghost,ghost_data); */

  double* dudr_local [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudr_local,d4est_elliptic_data->local_nodes);

  int ghost_nodes = d4est_mesh_get_ghost_nodes(d4est_ghost);
  double* dudr_ghost [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudr_ghost, ghost_nodes);
  
  d4est_poisson_compute_dudr
    (
     p4est,
     d4est_ghost,
     d4est_ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     dudr_local,
     dudr_ghost,
     d4est_elliptic_data->u,
     d4est_elliptic_data->local_nodes,
     which_field
    );

  flux_fcn_data->d4est_ghost_data = d4est_ghost_data;
  flux_fcn_data->d4est_ghost = d4est_ghost;
  flux_fcn_data->p4est = p4est;
  flux_fcn_data->which_field = which_field;
  flux_fcn_data->local_nodes = d4est_elliptic_data->local_nodes;
  flux_fcn_data->u = d4est_elliptic_data->u;
  flux_fcn_data->Au = d4est_elliptic_data->Au;
  for (int i = 0; i < (P4EST_DIM); i++){
    flux_fcn_data->dudr_local[i] = dudr_local[i];
    flux_fcn_data->dudr_ghost[i] = dudr_ghost[i];
  }


  /* DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_elliptic_data->Au, d4est_elliptic_data->local_nodes); */
  
  d4est_poisson_apply_mortar_matrices
    (
     p4est,
     d4est_ghost,
     /* d4est_ghost_data, */
     flux_fcn_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors
    );


  /* DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_elliptic_data->Au, d4est_elliptic_data->local_nodes); */
  
  D4EST_FREE_DIM_VEC(dudr_local);
  D4EST_FREE_DIM_VEC(dudr_ghost);
}
