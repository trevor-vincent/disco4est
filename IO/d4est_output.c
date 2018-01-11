#define _GNU_SOURCE
#include <d4est_output.h>
#include <d4est_vtk.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <d4est_mesh.h>
#include <d4est_estimator_stats.h>
#include <d4est_ip_energy_norm.h>
#include <sc_reduce.h>

static void
d4est_output_calculate_analytic_error
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 double* u,
 int local_nodes,
 d4est_xyz_fcn_t analytic_solution,
 void* analytic_solution_ctx,
 double* u_analytic,
 double* error
)
{
  d4est_mesh_init_field
    (
     p4est,
     u_analytic,
     analytic_solution,
     d4est_ops,
     d4est_geom,
     INIT_FIELD_ON_LOBATTO,
     analytic_solution_ctx
    );

  /* for (int i = 0; i < prob_vecs->local_nodes; i++) */
  /* printf("u_analytic[%d] = %.25f\n",i, u_analytic[i]); */
    
  d4est_linalg_vec_axpyeqz(-1., u, u_analytic, error, local_nodes);
  d4est_linalg_vec_fabs(error, local_nodes);
}

void
d4est_output_energy_norm_fit
(
 p4est_t* p4est,
 d4est_output_energy_norm_fit_t* fit
)
{
  double slope;
  double intercept;
  d4est_util_linear_regression
    (
     fit->log_energy_norm_data,
     fit->dof_data,
     &slope,
     &intercept,
     fit->stride
    );

  if (p4est->mpirank == 0){
    printf("[D4EST_OUTPUT_FIT](2): C1 = %f, C2 = %.15f\n", intercept, slope);
  }
}

void
d4est_output_energy_norm_add_entry_and_fit
(
 p4est_t* p4est,
 d4est_output_energy_norm_fit_t* fit,
 double global_energy_norm_sqr,
 double global_dof
){
  if (global_energy_norm_sqr > 0.){
    fit->log_energy_norm_data[fit->stride] = log(sqrt(global_energy_norm_sqr));
    fit->dof_data[fit->stride] = pow(global_dof, 1./(2.*(P4EST_DIM)-1.));
    fit->stride += 1;

    if (fit->stride >= 2){
      if (p4est->mpirank == 0){
        printf("[D4EST_OUTPUT_FIT](1): ||err|| = C1*exp(C2*DOF^(1/%d))\n",2*(P4EST_DIM)-1);
        printf("[D4EST_OUTPUT_FIT](1): ||err|| = %.15f\n",sqrt(global_energy_norm_sqr));
      }
      d4est_output_energy_norm_fit(p4est,fit);
    }
  }
}

void
d4est_output_destroy_energy_norm_fit
(
 d4est_output_energy_norm_fit_t* fit
)
{
  P4EST_FREE(fit->log_energy_norm_data);
  P4EST_FREE(fit->dof_data);
  P4EST_FREE(fit);
}

d4est_output_norms_t
d4est_output_norms
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_geometry_storage_t* d4est_factors,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_ip_energy_norm_data_t* energy_norm_data,
 double estimator,
 double* error,
 d4est_output_energy_norm_fit_t* fit,
 int (*skip_element_fcn)(d4est_element_data_t*)
)
{
  int local_nodes = 0;
  int local_quad_nodes = 0;

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
        local_nodes += d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        local_quad_nodes += d4est_lgl_get_nodes((P4EST_DIM), ed->deg_vol_quad);
      }
    }
  
  double local_l2_norm_sqr = d4est_mesh_compute_l2_norm_sqr
                             (
                              p4est,
                              d4est_ops,
                              d4est_geom,
                              d4est_quad,
                              error,
                              local_nodes,
                              DO_NOT_STORE_LOCALLY,
                              skip_element_fcn,
                              NULL
                             );

  double local_Linf = d4est_mesh_compute_linf
                      (
                       p4est,
                       error,
                       skip_element_fcn
                      );
      
  double local_energy_norm_sqr = -1.;
  if (energy_norm_data != NULL)
    local_energy_norm_sqr = d4est_ip_energy_norm_compute
                            (
                             p4est,
                             error,
                             energy_norm_data,
                             ghost,
                             ghost_data,
                             d4est_ops,
                             d4est_geom,
                             d4est_quad,
                             d4est_factors
                            );

  
  double local_nodes_dbl = (double)local_nodes;
  double local_quad_nodes_dbl = (double)local_quad_nodes;
  double local_estimator = estimator;
  double local_reduce [5];
  double global_reduce [5];
  double global_Linf;
    
  local_reduce[0] = local_l2_norm_sqr;
  local_reduce[1] = local_nodes_dbl;
  local_reduce[2] = local_quad_nodes_dbl;
  local_reduce[3] = local_estimator;
  local_reduce[4] = local_energy_norm_sqr;
    
  sc_reduce
    (
     &local_reduce[0],
     &global_reduce[0],
     5,
     sc_MPI_DOUBLE,
     sc_MPI_SUM,
     0,
     sc_MPI_COMM_WORLD
    );

    sc_reduce
    (
     &local_Linf,
     &global_Linf,
     1,
     sc_MPI_DOUBLE,
     sc_MPI_MAX,
     0,
     sc_MPI_COMM_WORLD
    );

  double global_l2_norm_sqr = global_reduce[0];
  double global_nodes_dbl = global_reduce[1];
  double global_quad_nodes_dbl = global_reduce[2];
  double global_estimator = global_reduce[3];
  double global_energy_norm_sqr = global_reduce[4];

  int avg_deg = pow((int)(global_nodes_dbl/p4est->global_num_quadrants), 1./(P4EST_DIM)) - 1.f + .5f;
  int avg_deg_quad = pow((int)(global_quad_nodes_dbl/p4est->global_num_quadrants), 1./(P4EST_DIM)) - 1.f + .5f;

  if(energy_norm_data != NULL && fit != NULL){
    d4est_output_energy_norm_add_entry_and_fit(p4est,fit, global_energy_norm_sqr, global_nodes_dbl);
  }
    
  if (p4est->mpirank == 0){
    printf
      (
       "[D4EST_OUTPUT]: global_elements %d global_nodes %d avg_deg %d global_quad_nodes %d avg_deg_quad %d global_estimator %.25f global_l2 %.25f global_enorm %.25f global_Linf %.25f\n",
       (int)p4est->global_num_quadrants,
       (int)global_nodes_dbl,
       avg_deg,
       (int)global_quad_nodes_dbl,
       avg_deg_quad,
       (global_estimator < 0) ? -1. : sqrt(global_estimator),
       (global_l2_norm_sqr < 0) ? -1 : sqrt(global_l2_norm_sqr),
       (global_energy_norm_sqr < 0) ? -1. :sqrt(global_energy_norm_sqr),
       global_Linf
      );
  }


  return (d4est_output_norms_t)
    {.global_nodes = (int)global_nodes_dbl,
        .global_num_quadrants = p4est->global_num_quadrants,
        .avg_deg = avg_deg,
        .global_quad_nodes = (int)global_quad_nodes_dbl,
        .avg_deg_quad = avg_deg_quad,
        .global_estimator = (global_estimator < 0) ? -1. : sqrt(global_estimator),
        .global_l2_norm_sqr = (global_l2_norm_sqr < 0) ? -1 : sqrt(global_l2_norm_sqr),
        .global_linf = global_Linf
        };
    
}




void
d4est_output_norms_using_analytic_solution
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_geometry_storage_t* d4est_factors,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 double local_estimator,
 d4est_elliptic_data_t* prob_vecs,
 d4est_ip_energy_norm_data_t* energy_norm_data,
 d4est_xyz_fcn_t analytic_solution,
 void* ctx,
 d4est_output_energy_norm_fit_t* fit,
 int (*skip_element_fcn)(d4est_element_data_t*)
)
{
  double* error = P4EST_ALLOC(double, prob_vecs->local_nodes);
  double* u_analytic = P4EST_ALLOC(double, prob_vecs->local_nodes);
  d4est_output_calculate_analytic_error(p4est, d4est_ops, d4est_geom, d4est_quad, prob_vecs->u, prob_vecs->local_nodes, analytic_solution, ctx, u_analytic, error);
  d4est_output_norms(p4est, d4est_ops, d4est_geom, d4est_quad, d4est_factors, ghost, ghost_data, energy_norm_data, local_estimator, error, fit, skip_element_fcn);
  P4EST_FREE(error);
  P4EST_FREE(u_analytic);
}
