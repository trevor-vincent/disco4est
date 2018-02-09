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
  double *field_value_errors,
  int num_nodes_local,
  void *norm_fcn_ctx
)
{
  d4est_norms_fcn_L2_ctx_t *ctx = norm_fcn_ctx;
  
  double norm_sq_local = d4est_mesh_compute_l2_norm_sqr(
    ctx->p4est,
    ctx->d4est_ops,
    ctx->d4est_geom,
    ctx->d4est_quad,
    field_value_errors,
    num_nodes_local,
    DO_NOT_STORE_LOCALLY,
    NULL,
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
  double *field_value_errors,
  int num_nodes_local,
  void *norm_fcn_ctx
)
{
  double norm_local = 0;
  for (int i = 0; i < num_nodes_local; i++) {
    if (field_value_errors[i] > norm_local)
      norm_local = field_value_errors[i];
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


void
d4est_norms_fcn_energy_fit
(
 p4est_t* p4est,
 d4est_norms_fcn_energy_fit_t* fit
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

  if (p4est->mpirank == 0) {
    zlog_category_t *c_default = zlog_get_category("d4est_norms_fcn_energy_fit");
    zlog_info(c_default, "(2): C1 = %f, C2 = %.15f", intercept, slope);
  }
}

void
d4est_norms_fcn_energy_add_entry_and_fit
(
 p4est_t* p4est,
 d4est_norms_fcn_energy_fit_t* fit,
 double global_energy_norm_sqr,
 double global_dof
){
  if (global_energy_norm_sqr > 0.){
    fit->log_energy_norm_data[fit->stride] = log(sqrt(global_energy_norm_sqr));
    fit->dof_data[fit->stride] = pow(global_dof, 1./(2.*(P4EST_DIM)-1.));
    fit->stride += 1;

    if (fit->stride >= 2){
      if (p4est->mpirank == 0){
        zlog_category_t *c_default = zlog_get_category("d4est_norms_fcn_energy_fit");
        zlog_info(c_default, "(1): ||err|| = C1*exp(C2*DOF^(1/%d))", 2*(P4EST_DIM)-1);
        zlog_info(c_default, "(1): ||err|| = %.15f", sqrt(global_energy_norm_sqr));
      }
      d4est_norms_fcn_energy_fit(p4est,fit);
    }
  }
}

void
d4est_norms_fcn_energy_destroy_fit
(
 d4est_norms_fcn_energy_fit_t* fit
)
{
  P4EST_FREE(fit->log_energy_norm_data);
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
  double *field_value_errors,
  int num_nodes_local,
  void *norm_fcn_ctx
)
{
  d4est_norms_fcn_energy_ctx_t *ctx = norm_fcn_ctx;
  
  double norm_sq_local = d4est_ip_energy_norm_compute(
    ctx->p4est,
    field_value_errors,
    ctx->energy_norm_data,
    ctx->ghost,
    ctx->ghost_data,
    ctx->d4est_ops,
    ctx->d4est_geom,
    ctx->d4est_quad,
    ctx->d4est_factors
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

  // Need num_nodes for fit
  int num_nodes;
  if (ctx->fit != NULL) {
    sc_reduce(
      &num_nodes_local,
      &num_nodes,
      1,
      sc_MPI_INT,
      sc_MPI_SUM,
      0,
      sc_MPI_COMM_WORLD
    );
  }
  
  // Perform energy norm fit
  if (ctx->fit != NULL && ctx->p4est->mpirank == 0) {
    d4est_norms_fcn_energy_add_entry_and_fit(
      ctx->p4est,
      ctx->fit,
      norm_sq,
      num_nodes
    );
  }

  return sqrt(norm_sq);
}

/*
  Extracts the energy estimator from the context to save alongside norms.
*/
double
d4est_norms_fcn_energy_estimator
(
  double *field_value_errors,
  int num_nodes_local,
  void *norm_fcn_ctx
)
{
  d4est_norms_fcn_energy_ctx_t *ctx = norm_fcn_ctx;
  double norm_sq_local = ctx->energy_estimator_sq_local;

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
  const char** norm_names
){
  int i_fields = -1;
  while (field_names[++i_fields] != NULL) {

    char *norms_output_category;
    asprintf(&norms_output_category, "norms_%s", field_names[i_fields]);
    zlog_category_t *c_norms_output = zlog_get_category(norms_output_category);
    free(norms_output_category);
    
    char *header;
    asprintf(&header, "num_quadrants num_nodes num_quad_nodes");
    int i_norms = -1;
    while (norm_names[++i_norms] != NULL) {
      asprintf(&header, "%s %s",header, norm_names[i_norms]);
    }
    zlog_info(c_norms_output, header);
    free(header);
  }
}

/*
  Computes norms over the computational grid and outputs them to files.
  
  - field_names: Null terminated list of names for the fields on the computational grid. The norms are written in a separate file for each field. Either provide `field_values_compare` or an `analytical_solution` for each field.
  - norm_names: Null terminated list of names for the norms to be computed for each field. These should correspond to the `norm_names` passed to an initial call to `d4est_norms_write_headers`.
  - norm_fcns: A function pointer for each norm declared in `norm_names`. These functions are responsible for actually computing the norms. You can provide custom functions, or use the implementations supplied in the `d4est_norms_fcn` namespace.
*/
void
d4est_norms_save
(
  p4est_t *p4est,
  const char** field_names,
  double **field_values,
  double **field_values_compare,
  d4est_xyz_fcn_t *analytical_solutions,
  void **analytical_solution_ctxs,
  const char **norm_names,
  d4est_norm_fcn_t *norm_fcns,
  void **norm_fcn_ctxs
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

  // Get number of nodes on the local process
  int num_nodes_local = d4est_mesh_get_local_nodes(p4est);
  int num_nodes_local_combined[2] = {
    num_nodes_local,
    d4est_mesh_get_local_quad_nodes(p4est)
  };
  int num_nodes_combined[2];
  sc_reduce(
    &num_nodes_local_combined,
    &num_nodes_combined,
    2,
    sc_MPI_INT,
    sc_MPI_SUM,
    0,
    sc_MPI_COMM_WORLD
  );
  int num_nodes = num_nodes_combined[0];
  int num_quad_nodes = num_nodes_combined[1];


  d4est_xyz_fcn_t analytical_solution_i = NULL; // Points to the provided function when iterating through the fields, if one is available.
  double *field_values_analytical = NULL; // Holds computed analytical values when iterating through the fields, if no compare values were provided.
  double *field_values_compare_i = NULL; // Points to the compare values, either computed and stored in field_values_analytical, or provided by the function argument.
  double *field_value_errors = P4EST_ALLOC(double, num_nodes_local); // Holds the computed errors between values and compare values when iterating through the fields.

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

      norms[i_norms] = (*(norm_fcns[i_norms]))(
        field_value_errors,
        num_nodes_local,
        norm_fcn_ctxs[i_norms]
      );

    }

    // Write to output file
    if (p4est->mpirank == 0) {

      // Construct filename
      char *norms_output_category;
      asprintf(&norms_output_category, "norms_%s", field_names[i_fields]);
      zlog_category_t *c_norms_output = zlog_get_category(norms_output_category);
      free(norms_output_category);

      char *output;
      asprintf(
        &output,
        "%d %d %d",
        (int)p4est->global_num_quadrants,
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

      zlog_info(c_norms_output, output);
    }
    
  }
  
  P4EST_FREE(field_values_analytical);
  P4EST_FREE(field_value_errors);
  
  if (p4est->mpirank == 0)
    zlog_debug(c_default, "Completed computing norms.");
}
