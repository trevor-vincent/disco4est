#define _GNU_SOURCE
#include <d4est_output.h>
#include <d4est_vtk.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <d4est_estimator_stats.h>
#include <sc_reduce.h>

static void
d4est_output_calculate_analytic_error
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_elliptic_data_t* prob_vecs,
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
       analytic_solution_ctx
      );
    
    d4est_linalg_vec_axpyeqz(-1., prob_vecs->u, u_analytic, error, prob_vecs->local_nodes);
    d4est_linalg_vec_fabs(error, prob_vecs->local_nodes);
}

static void
d4est_output_norms
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_estimator_stats_t* stats,
 double* error
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
        local_quad_nodes += d4est_lgl_get_nodes((P4EST_DIM), ed->deg_quad);
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
                                 DO_NOT_STORE_LOCALLY
                                );

    double local_nodes_dbl = (double)local_nodes;
    double local_quad_nodes_dbl = (double)local_quad_nodes;
    double local_estimator = stats->total;
    double local_reduce [4];
    double global_reduce [4];
    
    local_reduce[0] = local_l2_norm_sqr;
    local_reduce[1] = local_nodes_dbl;
    local_reduce[2] = local_quad_nodes_dbl;
    local_reduce[3] = local_estimator;
    
    sc_reduce
      (
       &local_reduce[0],
       &global_reduce[0],
       4,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       0,
       sc_MPI_COMM_WORLD
      );

    double global_l2_norm_sqr = global_reduce[0];
    double global_nodes_dbl = global_reduce[1];
    double global_quad_nodes_dbl = global_reduce[2];
    double global_estimator = global_reduce[3];

    double avg_deg = pow(global_nodes_dbl/p4est->global_num_quadrants, 1./(P4EST_DIM)) - 1;
    double avg_deg_quad = pow(global_quad_nodes_dbl/p4est->global_num_quadrants, 1./(P4EST_DIM)) - 1;
    
    if (p4est->mpirank == 0){
      printf
        (
         "[D4EST_OUTPUT]: global_elements %d global_nodes %d avg_deg %d global_quad_nodes %d avg_deg_quad %d global_estimator %.25f global_l2 %.25f\n",
         (int)p4est->global_num_quadrants,
         (int)global_nodes_dbl,
         (int)avg_deg,
         (int)global_quad_nodes_dbl,
         (int)avg_deg_quad,
         sqrt(global_estimator),
         sqrt(global_l2_norm_sqr)
        );
    }
}


void
d4est_output_norms_using_analytic_solution
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_estimator_stats_t* stats,
 d4est_elliptic_data_t* prob_vecs,
 d4est_xyz_fcn_t analytic_solution,
 void* ctx
)
{
  double* error = P4EST_ALLOC(double, prob_vecs->local_nodes);
  double* u_analytic = P4EST_ALLOC(double, prob_vecs->local_nodes);
  d4est_output_calculate_analytic_error(p4est, d4est_ops, d4est_geom, d4est_quad, prob_vecs, analytic_solution, ctx, u_analytic, error);
  d4est_output_norms(p4est, d4est_ops, d4est_geom, d4est_quad, stats, error);
  P4EST_FREE(error);
  P4EST_FREE(u_analytic);
}


static void
vtk_field_plotter
(
 d4est_vtk_context_t* vtk_ctx,
 void* user
)
{
  d4est_output_vtk_fields_t* vecs = user;
  vtk_ctx = d4est_vtk_write_dg_point_dataf
            (
             vtk_ctx,
             4,
             0,
             "u",
             vecs->u,
             "error",
             vecs->error,
             "u_compare",
             vecs->u_compare,
             "jacobian",
             vecs->jacobian,
             vtk_ctx
            );


   vtk_ctx = d4est_vtk_write_dg_cell_dataf
                (
                 vtk_ctx,
                 1,
                 1,
                 1,
                 0,
                 1,
                 0,
                 0,
                 vtk_ctx
                );

   if (vecs->eta2 != NULL){
   vtk_ctx = d4est_vtk_write_dg_cell_dataf
             (
              vtk_ctx,
              1,
              1,
              1,
              0,
              1,
              1,
              0,
              "eta",
              vecs->eta2,
              vtk_ctx
             );

   }
   else {
     vtk_ctx = d4est_vtk_write_dg_cell_dataf
               (
                vtk_ctx,
                1,
                1,
                1,
                0,
                1,
                0,
                0,
                vtk_ctx
               );
   }
}

static void
d4est_output_vtk
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 double* u,
 double* u_compare,
 double* error,
 const char* input_file,
 const char* save_as_prefix,
 int local_nodes,
 int level
)
{
    double* jacobian_lgl = P4EST_ALLOC(double, local_nodes);
    int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
    double* eta2_array = P4EST_ALLOC(double, p4est->local_num_quadrants);
    d4est_mesh_compute_jacobian_on_lgl_grid(p4est, d4est_ops, d4est_geom, jacobian_lgl);
    d4est_mesh_get_array_of_degrees(p4est, deg_array);
    d4est_mesh_get_array_of_estimators(p4est, eta2_array);

    d4est_output_vtk_fields_t vtk_nodal_vecs;
    vtk_nodal_vecs.u = u;
    vtk_nodal_vecs.u_compare = u_compare;
    vtk_nodal_vecs.error = error;
    vtk_nodal_vecs.jacobian = jacobian_lgl;
    vtk_nodal_vecs.eta2 = eta2_array;


    char* save_as;
    asprintf(&save_as,"%s_level_%d", save_as_prefix, level);
    
    /* vtk output */
    d4est_vtk_save_geometry_and_dg_fields
      (
       save_as,
       p4est,
       d4est_ops,
       deg_array,
       input_file,
       "d4est_vtk_geometry",
       vtk_field_plotter,
       (void*)&vtk_nodal_vecs
      );

    free(save_as);
    P4EST_FREE(jacobian_lgl);
    P4EST_FREE(deg_array);
    P4EST_FREE(eta2_array);

}


void
d4est_output_vtk_with_analytic_error
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_elliptic_data_t* prob_vecs,
 const char* input_file,
 const char* save_as_prefix,
 d4est_xyz_fcn_t analytic_solution,
 void* ctx,
 int level
)
{
  double* error = P4EST_ALLOC(double, prob_vecs->local_nodes);
  double* u_analytic = P4EST_ALLOC(double, prob_vecs->local_nodes);
  d4est_output_calculate_analytic_error(p4est, d4est_ops, d4est_geom, d4est_quad, prob_vecs, analytic_solution, ctx, u_analytic, error);

  d4est_output_vtk
    (
     p4est,
     d4est_ops,
     d4est_geom,
     prob_vecs->u,
     u_analytic,
     error,
     input_file,
     save_as_prefix,
     prob_vecs->local_nodes,
     level
    );
  
  P4EST_FREE(error);
  P4EST_FREE(u_analytic);
}
