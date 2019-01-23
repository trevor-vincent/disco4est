#define _GNU_SOURCE
#include <d4est_norms.h>
#include <d4est_vtk.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <d4est_mesh.h>
#include <d4est_estimator_stats.h>
#include <d4est_ip_energy_norm.h>
#include <sc_reduce.h>
#include <zlog.h>


/*
  Computes the L_2 norm.
  
  - TODO: Add definition
  - TODO: Merge with `d4est_mesh_compute_l2_norm_sqr`
*/
double
d4est_norms_fcn_L2
(
  p4est_t* p4est,
  double *field_value_errors,
  int num_nodes_local,
  void *norm_fcn_ctx,
  int (*skip_element_fcn)(d4est_element_data_t*)
)
{
  d4est_norms_fcn_L2_ctx_t *ctx = norm_fcn_ctx;
  
  double norm_sq_local = d4est_mesh_compute_l2_norm_sqr(
    ctx->p4est,
    ctx->d4est_ops,
    ctx->d4est_geom,
    ctx->d4est_quad,
    ctx->d4est_factors,
    field_value_errors,
    num_nodes_local,
    skip_element_fcn,
    NULL
  );
    
  // Reduce over all parallel processes
  double norm_sq;
  sc_reduce(
    &norm_sq_local,
    &norm_sq,
    1,
    sc_MPI_DOUBLE,
    sc_MPI_SUM,
    0,
    sc_MPI_COMM_WORLD
  );

  return sqrt(norm_sq);
}


/*
  Computes the L_infty norm defined by:
  
  L_infty(\vec{x}) = sup_i(x_i)
*/
double
d4est_norms_fcn_Linfty
(
 p4est_t* p4est,
 double *field_value_errors,
 int num_nodes_local,
 void *norm_fcn_ctx,
 int (*skip_element_fcn)(d4est_element_data_t*)
)
{
  double norm_local = 0;

  if (skip_element_fcn == NULL){
    for (int i = 0; i < num_nodes_local; i++) {
      if (field_value_errors[i] > norm_local)
        norm_local = field_value_errors[i];
    }
  }
  else {
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
          int skip_element = (skip_element_fcn != NULL) ? skip_element_fcn(ed) : 0;
          if (!skip_element){
            int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
            for (int i = 0; i < volume_nodes; i++) {
              if (field_value_errors[ed->nodal_stride + i] > norm_local)
                norm_local = field_value_errors[ed->nodal_stride + i];
            }
          }        
        }
      }
  }
  // Reduce over all parallel processes
  double norm;
  sc_reduce(
    &norm_local,
    &norm,
    1,
    sc_MPI_DOUBLE,
    sc_MPI_MAX,
    0,
    sc_MPI_COMM_WORLD
  );

  return norm;
}

d4est_norms_linear_fit_t*
d4est_norms_linear_fit_init
(
)
{
  d4est_norms_linear_fit_t* fit = P4EST_ALLOC(d4est_norms_linear_fit_t, 1);
  fit->log_norm_data = NULL;
  fit->dof_data = NULL;
  fit->num_of_data_entries = 0;
  fit->stride = 0;
  return fit;
}

void
d4est_norms_linear_fit
(
 p4est_t* p4est,
 d4est_norms_linear_fit_t* fit
)
{
  double slope;
  double intercept;
  d4est_util_linear_regression
    (
     fit->log_norm_data,
     fit->dof_data,
     &slope,
     &intercept,
     fit->stride
    );

  if (p4est->mpirank == 0) {
    zlog_category_t *c_default = zlog_get_category("d4est_norms_linear_fit");
    zlog_info(c_default, "C1 = %f, C2 = %.15f", intercept, slope);
  }
}

void
d4est_norms_linear_fit_add_entry_and_fit
(
 p4est_t* p4est,
 d4est_norms_linear_fit_t* fit,
 double global_norm_sqr,
 double global_dof
){
  if (global_norm_sqr > 0.){
    fit->num_of_data_entries += 1;
    fit->log_norm_data = P4EST_REALLOC(fit->log_norm_data, double, fit->num_of_data_entries);
    fit->dof_data = P4EST_REALLOC(fit->dof_data, double, fit->num_of_data_entries);
    
    fit->log_norm_data[fit->stride] = log(sqrt(global_norm_sqr));
    fit->dof_data[fit->stride] = global_dof;
    fit->stride += 1;

    if (fit->stride >= 2){
      if (p4est->mpirank == 0){
        zlog_category_t *c_default = zlog_get_category("d4est_norms_linear_fit");
        zlog_info(c_default, "||err|| = C1*exp(C2*DOF)");
        zlog_info(c_default, "||err|| = %.15f", sqrt(global_norm_sqr));
      }
      d4est_norms_linear_fit(p4est,fit);
    }
  }
}

void
d4est_norms_linear_fit_destroy
(
 d4est_norms_linear_fit_t* fit
)
{
  P4EST_FREE(fit->log_norm_data);
  P4EST_FREE(fit->dof_data);
  P4EST_FREE(fit);
}


/*
  Computes the energy norm for elliptic PDEs.
  
  - TODO: Add definition
*/
double
d4est_norms_fcn_energy
(
  p4est_t* p4est,
  double *field_value_errors,
  int num_nodes_local,
  void *norm_fcn_ctx,
  int (*skip_element_fcn)(d4est_element_data_t*)
)
{
  if (skip_element_fcn != NULL){
    D4EST_ABORT("Do not use d4est_norms_fcn_energy when skip_element_fcn != NULL, (not supported)");
  }
  
  d4est_norms_fcn_energy_ctx_t *ctx = norm_fcn_ctx;
  
  double norm_sq_local = d4est_ip_energy_norm_compute
                         (
                          ctx->p4est,
                          field_value_errors,
                          ctx->energy_norm_data,
                          ctx->ghost,
                          ctx->ghost_data,
                          ctx->d4est_ops,
                          ctx->d4est_geom,
                          ctx->d4est_quad,
                          ctx->d4est_factors,
                          ctx->which_field
                         );

  // Reduce over all parallel processes
  double norm_sq;
  sc_reduce(
    &norm_sq_local,
    &norm_sq,
    1,
    sc_MPI_DOUBLE,
    sc_MPI_SUM,
    0,
    sc_MPI_COMM_WORLD
  );

  return sqrt(norm_sq);
}

/*
  Extracts the energy estimator from the context to save alongside norms.
*/
double
d4est_norms_fcn_energy_estimator
(
  p4est_t* p4est,
  double *field_value_errors,
  int num_nodes_local,
  void *norm_fcn_ctx,
  int (*skip_element_fcn)(d4est_element_data_t*)
)
{
  double norm_sq_local = 0.;
  d4est_norms_fcn_energy_ctx_t *ctx = norm_fcn_ctx;
  double* estimator = ctx->energy_estimator;
  
  if (skip_element_fcn == NULL){
    norm_sq_local = ctx->energy_estimator_sq_local;
  }
  else {
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
          int skip_element = (skip_element_fcn != NULL) ? skip_element_fcn(ed) : 0;
          if (!skip_element){
            norm_sq_local += estimator[ed->id];
          }        
        }
      }
  }
  // Reduce over all parallel processes
  double norm_sq;
  sc_reduce(
    &norm_sq_local,
    &norm_sq,
    1,
    sc_MPI_DOUBLE,
    sc_MPI_SUM,
    0,
    sc_MPI_COMM_WORLD
  );

  return sqrt(norm_sq);
}


/*
  Writes column headers to norm data output files.
  
  > Note: Call from one process only to avoid duplicate lines.
  
  - field_names: Null terminated list of names for the fields on the computational grid. The norms are written in a separate file for each field.
  - norm_names: Null terminated list of names for the norms to be computed for each field. These are the column headers association to each norm.
*/
void
d4est_norms_write_headers
(
  const char** field_names,
  const char** norm_names,
  const char* zlog_output_category_prefix
){
  int i_fields = -1;
  while (field_names[++i_fields] != NULL) {

    char *norms_output_category;
    if (zlog_output_category_prefix == NULL)
      asprintf(&norms_output_category, "d4est_region_all_norms_%s", field_names[i_fields]);
    else
      asprintf(&norms_output_category, "d4est_%s_norms_%s", zlog_output_category_prefix, field_names[i_fields]);
      
    zlog_category_t *c_norms_output = zlog_get_category(norms_output_category);
    free(norms_output_category);
    
    char *header;
    asprintf(&header, "num_quadrants num_nodes num_quad_nodes");
    int i_norms = -1;
    while (norm_names[++i_norms] != NULL) {
      asprintf(&header, "%s %s",header, norm_names[i_norms]);
    }
    zlog_info(c_norms_output, "%s", header);
    free(header);
  }
}

/** 
 * Writes norms out for either the entire grid, or a region
 * based on a non-null skip_element_fcn passed to d4est_norms_save
 * 
 * @param p4est 
 * @param d4est_factors 
 * @param field_names Null terminated list of names for the fields on the computational grid. The norms are written in a separate file for each field. Either provide `field_values_compare` or an `analytical_solution` for each field.
 * @param field_values 
 * @param field_values_compare 
 * @param analytical_solutions 
 * @param analytical_solution_ctxs 
 * @param norm_names Null terminated list of names for the norms to be computed for each field. These should correspond to the `norm_names` passed to an initial call to `d4est_norms_write_headers`
 * @param norm_fcns norm_fcns: A function pointer for each norm declared in `norm_names`. These functions are responsible for actually computing the norms. You can provide custom functions, or use the implementations supplied in the `d4est_norms_fcn` namespace.
 * @param norm_fcn_ctxs 
 * @param linear_fits 
 * @param zlog_output_category_prefix 
 * @param skip_element_fcn returns 1 if an element's norm should not be included, 0 otherwise, this is
 * useful when you want to only compute norms on subspaces of the mesh.
 */
void
d4est_norms_save
(
  p4est_t *p4est,
  d4est_mesh_data_t* d4est_factors,
  const char** field_names,
  double **field_values,
  double **field_values_compare,
  d4est_xyz_fcn_t *analytical_solutions,
  void **analytical_solution_ctxs,
  const char **norm_names,
  d4est_norm_fcn_t *norm_fcns,
  void **norm_fcn_ctxs,
  d4est_norms_linear_fit_t** linear_fits,
  const char* zlog_output_category_prefix,
  int (*skip_element_fcn)(d4est_element_data_t*) /* for skipping elements not in region */
)
{
  zlog_category_t *c_default = zlog_get_category("d4est_norms");
  if (p4est->mpirank == 0)
    zlog_debug(c_default, "Computing norms...");

  // Get the number of fields and norms
  int num_fields = -1;
  while (field_names[++num_fields] != NULL)
    continue;
  int num_norms = -1;
  while (norm_names[++num_norms] != NULL)
    continue;


  int local_nodes_in_region = 0;
  int local_quad_nodes_in_region = 0;
  int local_elements_in_region = 0;

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
        int skip_element = (skip_element_fcn != NULL) ? skip_element_fcn(ed) : 0;

        if (!skip_element){
          local_elements_in_region += 1;
          local_nodes_in_region += d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
          local_quad_nodes_in_region += d4est_lgl_get_nodes((P4EST_DIM), ed->deg_quad);
        }
      }
    }

  int num_nodes_local_combined[3] = {
    local_elements_in_region,
    local_nodes_in_region,
    local_quad_nodes_in_region
  };
  
  int num_nodes_combined[3];
  sc_reduce(
    &num_nodes_local_combined,
    &num_nodes_combined,
    3,
    sc_MPI_INT,
    sc_MPI_SUM,
    0,
    sc_MPI_COMM_WORLD
  );
  int num_elements = num_nodes_combined[0];
  int num_nodes = num_nodes_combined[1];
  int num_quad_nodes = num_nodes_combined[2];

  d4est_xyz_fcn_t analytical_solution_i = NULL; // Points to the provided function when iterating through the fields, if one is available.
  double *field_values_analytical = NULL; // Holds computed analytical values when iterating through the fields, if no compare values were provided.
  double *field_values_compare_i = NULL; // Points to the compare values, either computed and stored in field_values_analytical, or provided by the function argument.
  int num_nodes_local = d4est_factors->local_sizes.local_nodes;
  double *field_value_errors = P4EST_ALLOC(double,
                                           num_nodes_local); // Holds the computed errors between values and compare values when iterating through the fields.

  // Collect local measure of norms
  for (int i_fields = 0; i_fields < num_fields; i_fields++){
    
    // Use compare values if available
    field_values_compare_i = field_values_compare[i_fields];
    analytical_solution_i = analytical_solutions[i_fields];
    // Else, compute compare values from analytic function
    if (field_values_compare_i == NULL && analytical_solution_i != NULL) {
      if (field_values_analytical == NULL) {
        field_values_analytical = P4EST_ALLOC(double, num_nodes_local);
      }
      d4est_mesh_init_field(
        p4est,
        field_values_analytical,
        analytical_solution_i,
        NULL,
        NULL,
        d4est_factors,
        INIT_FIELD_ON_LOBATTO,
        analytical_solution_ctxs[i_fields]
      );
      field_values_compare_i = field_values_analytical;
    } else if (field_values_compare_i == NULL) {
      zlog_error(c_default, "Could not compute norms for %s since neither compare values nor an analytic solution is available.", field_names[i_fields]);
      continue;
    }
    
    // Compute errors
    d4est_linalg_vec_axpyeqz(-1., field_values[i_fields], field_values_compare_i, field_value_errors, num_nodes_local);
    d4est_linalg_vec_fabs(field_value_errors, num_nodes_local);
    
    // Compute norms
    double norms[num_norms];
    for (int i_norms=0; i_norms < num_norms; i_norms++) {

      norms[i_norms] = (*(norm_fcns[i_norms]))
                       (
                        p4est,
                        field_value_errors,
                        num_nodes_local,
                        norm_fcn_ctxs[i_norms],
                        skip_element_fcn
                       );

    }

    if (linear_fits != NULL){
      for (int i_norms=0; i_norms < num_norms; i_norms++) {
        if (linear_fits[i_norms] != NULL){
          d4est_norms_linear_fit_add_entry_and_fit
            (
             p4est,
             linear_fits[i_norms],
             norms[i_norms],
             num_nodes
            );
        }
      }        
    }
    

    // Write to output file
    if (p4est->mpirank == 0) {

      // Construct filename
 
    char *norms_output_category;
    if (zlog_output_category_prefix == NULL)
      asprintf(&norms_output_category, "d4est_region_all_norms_%s", field_names[i_fields]);
    else
      asprintf(&norms_output_category, "d4est_%s_norms_%s", zlog_output_category_prefix, field_names[i_fields]);
      
      zlog_category_t *c_norms_output = zlog_get_category(norms_output_category);
      free(norms_output_category);

      char *output;
      asprintf(
        &output,
        "%d %d %d",
        (int)num_elements,
        num_nodes,
        num_quad_nodes
      );

      for (int i_norms=0; i_norms < num_norms; i_norms++) {
        asprintf(
          &output,
          "%s %.25f",
          output,
          norms[i_norms]
        );
      }

      zlog_info(c_norms_output, "%s", output);
    }
    
  }
  
  P4EST_FREE(field_values_analytical);
  P4EST_FREE(field_value_errors);
  
  if (p4est->mpirank == 0)
    zlog_debug(c_default, "Completed computing norms.");
}
