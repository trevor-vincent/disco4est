#define _GNU_SOURCE
#include <d4est_output.h>
#include <d4est_vtk.h>
#include <d4est_util.h>
#include <sc_reduce.h>

void
d4est_output_error
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 int level,
 double* error,
 double local_estimator,
 int local_nodes
)
{

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
    double local_reduce [3];
    double global_reduce [3];

    local_reduce[0] = local_l2_norm_sqr;
    local_reduce[1] = local_nodes_dbl;
    local_reduce[2] = local_estimator;
    
    sc_reduce
      (
       &local_reduce[0],
       &global_reduce[0],
       3,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       0,
       sc_MPI_COMM_WORLD
      );

    double global_l2_norm_sqr = global_reduce[0];
    double global_nodes_dbl = global_reduce[1];
    double global_estimator = global_reduce[2];

    if (p4est->mpirank == 0){
      printf
        (
         "[D4EST_OUTPUT]: level %d global_elements %d global_nodes %d global_estimator %.25f global_l2 %.25f\n",
         level,
         (int)p4est->global_num_quadrants,
         (int)global_nodes_dbl,
         sqrt(global_estimator),
         sqrt(global_l2_norm_sqr)
        );
    }
    
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
             3,
             0,
             "u",
             vecs->u,
             "error",
             vecs->error,
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

void
d4est_output_vtk
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 int level,
 double* u,
 double* error,
 int local_nodes,
 const char* input_file,
 const char* save_as_prefix,
 d4est_output_estimator_option_t eta2_option
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
    vtk_nodal_vecs.error = error;
    vtk_nodal_vecs.jacobian = jacobian_lgl;
    
    if (eta2_option == SAVE_ESTIMATOR){
      vtk_nodal_vecs.eta2 = eta2_array;
    }
    else{
      vtk_nodal_vecs.eta2 = NULL;
    }

    char* save_as;
    D4EST_ASPRINTF(save_as,"%s_level_%d", save_as_prefix, level);
    
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

    P4EST_FREE(save_as);
    P4EST_FREE(jacobian_lgl);
    P4EST_FREE(deg_array);
    P4EST_FREE(eta2_array);

}
