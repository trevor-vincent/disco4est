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

d4est_output_energy_norm_fit_t*
d4est_output_new_energy_norm_fit
(
 int num_of_entries
)
{
  d4est_output_energy_norm_fit_t* fit = P4EST_ALLOC(d4est_output_energy_norm_fit_t, 1);

  fit->log_energy_norm_data = P4EST_ALLOC(double, num_of_entries);
  fit->dof_data = P4EST_ALLOC(double, num_of_entries);
  fit->stride = 0;
  
  return fit;
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


void
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


static void
vtk_dg_field_plotter
(
 d4est_vtk_context_t* vtk_ctx,
 void* user
)
{
  d4est_output_vtk_dg_fields_t* vecs = user;

  if (vecs->residual == NULL){
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
  }
  else {
  vtk_ctx = d4est_vtk_write_dg_point_dataf
            (
             vtk_ctx,
             5,
             0,
             "u",
             vecs->u,
             "error",
             vecs->error,
             "u_compare",
             vecs->u_compare,
             "jacobian",
             vecs->jacobian,
             "residual",
             vecs->residual,
             vtk_ctx
            );
  }

  /* vtk_ctx = d4est_vtk_write_dg_cell_dataf */
  /*              ( */
  /*               vtk_ctx, */
  /*               1, */
  /*               1, */
  /*               1, */
  /*               0, */
  /*               1, */
  /*               0, */
  /*               0, */
  /*               vtk_ctx */
  /*              ); */

  if (vecs->eta2 != NULL){
    vtk_ctx = d4est_vtk_write_dg_cell_dataf
              (
               vtk_ctx,
               1,
               1,
               1,
               0,
               1,
               2,
               0,
               "eta",
               vecs->eta2,
               "l2",
               vecs->l2,
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
               1,
               0,
               "l2",
               vecs->l2,
               vtk_ctx
              );
  }
}

void
d4est_output_vtk_with_no_fields
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 const char* input_file,
 const char* save_as_prefix,
 int level
){
  int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
  d4est_mesh_get_array_of_degrees(p4est, (void*)deg_array, D4EST_INT);

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
     NULL,
     NULL
    );

  free(save_as);
  P4EST_FREE(deg_array);
}

void
d4est_output_vtk
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 double* u,
 double* u_compare,
 double* error,
 double* residual,
 const char* input_file,
 const char* save_as_prefix,
 int local_nodes,
 int level,
 int save_estimator
)
{
  double* jacobian_lgl = P4EST_ALLOC(double, local_nodes);
  int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);

  d4est_mesh_compute_jacobian_on_lgl_grid(p4est, d4est_ops, d4est_geom, jacobian_lgl);
  d4est_mesh_get_array_of_degrees(p4est, (void*)deg_array, D4EST_INT);

    
  double* eta2_array = NULL;
  if(save_estimator){
    eta2_array = P4EST_ALLOC(double, p4est->local_num_quadrants);
    d4est_mesh_get_array_of_estimators(p4est, eta2_array);
  }

  double* l2_array = NULL;
  l2_array = P4EST_ALLOC(double, p4est->local_num_quadrants);
  d4est_mesh_compute_l2_norm_sqr(p4est,
                                 d4est_ops,
                                 d4est_geom,
                                 d4est_quad,
                                 error,
                                 local_nodes,
                                 DO_NOT_STORE_LOCALLY,
                                 NULL,
                                 l2_array);
    
  d4est_output_vtk_dg_fields_t vtk_nodal_vecs;
  vtk_nodal_vecs.u = u;
  vtk_nodal_vecs.u_compare = u_compare;
  vtk_nodal_vecs.error = error;
  vtk_nodal_vecs.jacobian = jacobian_lgl;
  vtk_nodal_vecs.eta2 = eta2_array;
  vtk_nodal_vecs.l2 = l2_array;
  vtk_nodal_vecs.residual = residual;


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
     vtk_dg_field_plotter,
     (void*)&vtk_nodal_vecs
    );

  free(save_as);
  P4EST_FREE(jacobian_lgl);
  P4EST_FREE(deg_array);
  P4EST_FREE(eta2_array);
  P4EST_FREE(l2_array);

}


static void
vtk_cell_field_plotter
(
 d4est_vtk_context_t* vtk_ctx,
 void* user
)
{
  d4est_output_vtk_cell_fields_t* vecs = user;

  if (vecs->eta2 != NULL){
    vtk_ctx = d4est_vtk_write_dg_cell_dataf
              (
               vtk_ctx,
               1,
               1,
               1,
               0,
               1,
               4,
               0,
               "degree_dbl",
               vecs->deg,
               "degree_quad_dbl",
               vecs->deg_quad,
               "l2",
               vecs->l2,
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
               3,
               0,
               "degree_dbl",
               vecs->deg,
               "degree_quad_dbl",
               vecs->deg_quad,
               "l2",
               vecs->l2,
               vtk_ctx
              );
  }
}

void
d4est_output_vtk_degree_mesh_with_analytic_error
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_elliptic_data_t* prob_vecs,
 d4est_xyz_fcn_t analytic_solution,
 void* ctx,
 int local_nodes,
 const char* input_file,
 const char* save_as_prefix,
 int save_estimator,
 int level
){

  double* error = P4EST_ALLOC(double, prob_vecs->local_nodes);
  double* u_analytic = P4EST_ALLOC(double, prob_vecs->local_nodes);
  d4est_output_calculate_analytic_error(p4est, d4est_ops, d4est_geom, d4est_quad, prob_vecs->u, prob_vecs->local_nodes, analytic_solution, ctx, u_analytic, error);

  d4est_output_vtk_degree_mesh
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     error,
     local_nodes,
     input_file,
     save_as_prefix,
     save_estimator,
     level
    );

  P4EST_FREE(error);
  P4EST_FREE(u_analytic);
}

void
d4est_output_vtk_degree_mesh
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 double* error,
 int local_nodes,
 const char* input_file,
 const char* save_as_prefix,
 int save_estimator,
 int level
){
  double* deg_array = P4EST_ALLOC(double, p4est->local_num_quadrants);
  double* deg_quad_array = P4EST_ALLOC(double, p4est->local_num_quadrants);
  d4est_mesh_get_array_of_degrees(p4est, (void*)deg_array, D4EST_DOUBLE);
  d4est_mesh_get_array_of_quadrature_degrees(p4est, (void*)deg_quad_array, D4EST_DOUBLE);
    
  double* eta2_array = NULL;
  if(save_estimator){
    eta2_array = P4EST_ALLOC(double, p4est->local_num_quadrants);
    d4est_mesh_get_array_of_estimators(p4est, eta2_array);
  }

  double* l2_array = NULL;
  l2_array = P4EST_ALLOC(double, p4est->local_num_quadrants);
  d4est_mesh_compute_l2_norm_sqr(p4est,
                                 d4est_ops,
                                 d4est_geom,
                                 d4est_quad,
                                 error,
                                 local_nodes,
                                 DO_NOT_STORE_LOCALLY,
                                 NULL,
                                 l2_array);

  
  d4est_output_vtk_cell_fields_t vtk_cell_vecs;
  vtk_cell_vecs.eta2 = eta2_array;
  vtk_cell_vecs.l2 = l2_array;
  vtk_cell_vecs.deg = deg_array;
  vtk_cell_vecs.deg_quad = deg_quad_array;
    
  char* save_as;
  asprintf(&save_as,"%s_level_%d", save_as_prefix, level);
    
  /* vtk output */
  d4est_vtk_save_geometry_and_cell_fields
    (
     save_as,
     p4est,
     d4est_ops,
     input_file,
     "d4est_vtk_geometry",
     vtk_cell_field_plotter,
     (void*)&vtk_cell_vecs
    );

  free(save_as);
  P4EST_FREE(deg_array);
  P4EST_FREE(deg_quad_array);
  P4EST_FREE(eta2_array); 
  P4EST_FREE(l2_array); 
  
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
 int save_estimator,
 int level
)
{
  double* error = P4EST_ALLOC(double, prob_vecs->local_nodes);
  double* u_analytic = P4EST_ALLOC(double, prob_vecs->local_nodes);
  d4est_output_calculate_analytic_error(p4est, d4est_ops, d4est_geom, d4est_quad, prob_vecs->u, prob_vecs->local_nodes, analytic_solution, ctx, u_analytic, error);

  /* double error_sum = 0.; */
  /* for (int i = 0; i < prob_vecs->local_nodes; i++) */
  /*   error_sum += error[i]; */
  /* printf("error sum = %.25f\n", error_sum); */

  
  d4est_output_vtk
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     prob_vecs->u,
     u_analytic,
     error,
     NULL,
     input_file,
     save_as_prefix,
     prob_vecs->local_nodes,
     level,
     save_estimator
    );

  P4EST_FREE(error);
  P4EST_FREE(u_analytic);
}


/* } */
