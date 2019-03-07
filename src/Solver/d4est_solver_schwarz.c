#include <d4est_solver_schwarz.h>
#include <d4est_solver_schwarz_helpers.h>
#include <d4est_solver_schwarz_transfer_ghost_data.h>
#include <d4est_vtk.h>
/** 
 * Init all the data needed for solving
 * with the schwarz method (additive)
 * 
 * @param p4est 
 * @param d4est_ops 
 * @param d4est_ghost 
 * @param d4est_factors 
 * @param input_file 
 * @param schwarz_ops if user does not provide operators, 
  set this pointer to NULL so that it will be computed
 * 
 * @return 
 */

d4est_solver_schwarz_t*
d4est_solver_schwarz_init
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_ghost_t* d4est_ghost,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_solver_schwarz_apply_lhs_t* apply_lhs,
 const char* input_file,
 const char* input_section
)
{
  d4est_solver_schwarz_t* schwarz
    = P4EST_ALLOC(d4est_solver_schwarz_t, 1);


  schwarz->metadata =
    d4est_solver_schwarz_metadata_init
    (
     p4est,
     d4est_ghost,
     input_file,
     input_section
    );

  schwarz->used_as_smoother = 0;

  if (schwarz_ops == NULL){
    schwarz->operators = d4est_solver_schwarz_operators_init
                              (
                               d4est_ops
                              );
    schwarz->operators_external_alloc = 0;
  }
  else {
    schwarz->operators_external_alloc = 1;
    schwarz->operators = schwarz_ops;
  }
  
  schwarz->geometric_data =
    d4est_solver_schwarz_geometric_data_init
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_ghost,
     d4est_factors,
     schwarz->metadata,
     input_file,
     input_section
    );

  schwarz->apply_lhs = apply_lhs;
  schwarz->subdomain_solver = d4est_solver_schwarz_subdomain_solver_init
                              (
                               p4est,
                               input_file,
                               input_section
                              );
  
  d4est_field_type_t field_type = NODAL;
  schwarz->residual_ghost_data = d4est_ghost_data_init(
                                                       p4est,
                                                       d4est_ghost,
                                                       &field_type,
                                                       1
                                                      );
;
  schwarz->correction_ghost_data = NULL;

  /* schwarz->collect_subdomain_solve_info = 1; */
  /* if (schwarz->collect_subdomain_solve_info){ */
    schwarz->subdomain_solve_residuals = P4EST_ALLOC_ZERO(double,
                                                     p4est->local_num_quadrants);

    schwarz->subdomain_solve_iterations = P4EST_ALLOC_ZERO(int,
                                                      p4est->local_num_quadrants);
  /* } */


  
  return schwarz;
}

void
d4est_solver_schwarz_destroy
(
 d4est_solver_schwarz_t* schwarz
)
{
  d4est_solver_schwarz_metadata_destroy(schwarz->metadata);

  if(schwarz->operators_external_alloc == 0){
    d4est_solver_schwarz_operators_destroy(schwarz->operators);
  }
  
  d4est_ghost_data_ext_destroy
    (
     schwarz->correction_ghost_data
    );

  d4est_ghost_data_destroy
    (
     schwarz->residual_ghost_data
    );
  
  d4est_solver_schwarz_geometric_data_destroy
    (
     schwarz->geometric_data
    );
  
  d4est_solver_schwarz_subdomain_solver_destroy
    (
     schwarz->subdomain_solver
    );


  /* if (schwarz->collect_subdomain_solve_info){ */
    P4EST_FREE(schwarz->subdomain_solve_residuals);
    P4EST_FREE(schwarz->subdomain_solve_iterations);
  /* } */

  
  P4EST_FREE(schwarz);
}


/** 
 * Solves Au = b
 * by solving Ae = r on subdomains
 * and then adding the correction e to u
 * User must provide both u and r, b is not
 * needed because r is given and used to solve
 * Ae = r. Both u and r are nodal mesh arrays
 * i.e. they live on the nodal space of the mesh
 * not on the schwarz subdomains. Currently
 * cannot solve coupled equations, that is to be added later.
 * 
 * @param p4est 
 * @param d4est_geom 
 * @param d4est_quad 
 * @param d4est_factors 
 * @param ghost 
 * @param schwarz_data 
 * @param subdomain 
 * @param u  nodal initial guess for Au = b 
 * @param r  nodal r residual = b - Au 
 * @param r_ghost_data ghost data for nodal r 
 * 
 */
void
d4est_solver_schwarz_iterate
(
 p4est_t* p4est,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* ghost,
 d4est_solver_schwarz_t* schwarz,
 d4est_elliptic_data_t* vecs,
 /* double* u, */
 double* r
)
{
  double* u = vecs->u;
  d4est_solver_schwarz_operators_t* schwarz_ops
    = schwarz->operators;

  d4est_solver_schwarz_metadata_t* schwarz_metadata
    = schwarz->metadata;
    
  int local_nodes = d4est_factors->local_sizes.local_nodes;
  double* du = P4EST_ALLOC_ZERO(double, schwarz_metadata->restricted_nodal_size); 
  double* restricted_r = P4EST_ALLOC_ZERO(double, schwarz_metadata->restricted_nodal_size);
  double* correction_field_over_subdomains = P4EST_ALLOC_ZERO(double, schwarz_metadata->nodal_size);

  d4est_ghost_data_exchange
    (
     p4est,
     ghost,
     schwarz->residual_ghost_data,
     r
    );
  d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomains
    (
     p4est,
     schwarz_ops->d4est_ops,
     d4est_factors,
     ghost,
     schwarz->residual_ghost_data,
     schwarz_metadata,
     schwarz_ops,
     r,
     0, /* only one field for uncoupled equations */
     restricted_r
    );

  if(schwarz->apply_lhs->get_ctx_fcn != NULL){
    schwarz->apply_lhs->get_ctx_fcn
      (
       p4est,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       ghost,
       schwarz_ops,
       schwarz_metadata,
       schwarz->geometric_data,
       vecs,
       schwarz->apply_lhs->apply_lhs_ctx
      );
  }
  
  for (int i = 0; i < schwarz_metadata->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_metadata->subdomain_metadata[i];
    /* if (i == 1){ */

      /* double* du_temp = &du[sub_data.restricted_nodal_stride]; */
      /* DEBUG_PRINT_ARR_DBL(du_temp, sub_data.restricted_nodal_size); */

    d4est_solver_schwarz_subdomain_solver_info_t info = schwarz->subdomain_solver->solver_fcn
      (
       p4est,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       ghost,
       schwarz_ops,
       schwarz_metadata,
       schwarz->geometric_data,
       schwarz->apply_lhs,
       &du[sub_data.restricted_nodal_stride],
       &restricted_r[sub_data.restricted_nodal_stride],
       i,
       schwarz->subdomain_solver->solver_params
      );

    /* printf("sub_data.coreId, info.final_res, info.final_iter = %d, %f, %d\n", sub_data.core_id, info.final_res, info.final_iter); */
    schwarz->subdomain_solve_residuals[sub_data.core_id] = info.final_res;
    schwarz->subdomain_solve_iterations[sub_data.core_id] = info.final_iter;
    /* } */
  }

  d4est_solver_schwarz_compute_correction
    (
     schwarz_metadata,
     schwarz_ops,
     du,
     correction_field_over_subdomains
    );

  d4est_solver_schwarz_transfer_ghost_data_and_add_corrections
    (
     p4est,
     ghost,
     schwarz_metadata,
     &schwarz->correction_ghost_data,
     u,
     correction_field_over_subdomains
    );
  
  P4EST_FREE(restricted_r);
  P4EST_FREE(correction_field_over_subdomains);
  P4EST_FREE(du);  
}

void
d4est_solver_schwarz_debug_vtk
(
 p4est_t* p4est,
 d4est_solver_schwarz_t* schwarz,
 char* input_file,
 char* input_section,
 char* save_as,
 char* folder,
 int sub_folder_number,
 const char ** dg_field_names,
 double ** dg_fields
)
{
  d4est_vtk_save_aux
    (
     p4est,
     schwarz->operators->d4est_ops,
     input_file,
     input_section,
     dg_field_names,
     dg_fields,
     (const char * []){"subdomain_solve_residuals",NULL},
     (double* []){schwarz->subdomain_solve_residuals},
     (const char * []){"subdomain_solve_iterations",NULL},
     (int* []){schwarz->subdomain_solve_iterations},
     folder,
     sub_folder_number
    );    
}
