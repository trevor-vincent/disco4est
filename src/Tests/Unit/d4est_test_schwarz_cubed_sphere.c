#include <pXest.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_operators.h>
#include <d4est_amr.h>
#include <d4est_amr_random.h>
#include <d4est_vtk.h>
#include <d4est_h5.h>
#include <d4est_checkpoint.h>
#include <d4est_element_data.h>
#include <d4est_laplacian.h>
#include <d4est_laplacian_flux.h>
#include <d4est_solver_schwarz.h>
#include <d4est_solver_schwarz_helpers.h>
#include <petscsnes.h>
#include <zlog.h>

static void
d4est_test_build_residual
(
 p4est_t* p4est,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* user
)
{
  d4est_laplacian_flux_data_t* flux_data_for_apply_lhs = user;  
  d4est_laplacian_apply_aij(p4est, ghost, ghost_data, prob_vecs, flux_data_for_apply_lhs, d4est_ops, d4est_geom, d4est_quad, d4est_factors, 0);
  d4est_linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);
}

static void
d4est_test_apply_lhs
(
 p4est_t* p4est,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* user
)
{
  d4est_laplacian_flux_data_t* flux_data_for_apply_lhs = user;  
  d4est_laplacian_apply_aij(p4est, ghost, ghost_data, prob_vecs, flux_data_for_apply_lhs, d4est_ops, d4est_geom, d4est_quad, d4est_factors, 0);
}

double
poly_vec_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
){
  return exp(x + y + z)*(x*x + y*y + z*z - 1.);
}

double
boundary_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void *user
)
{
  return poly_vec_fcn(x,
                      y,
#if(P4EST_DIM)==3
                      z,
#endif
                      user);
}

double
neg_laplacian_poly_vec_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
){
  /* return -4.; */
  /* double pi = 3.1415926535897932384626433832795; */
  /* return 4*pi*(-cos(pi*(x*x + y*y)) + pi*(x*x + y*y)*sin(pi*(x*x + y*y))); */
  return -1.*exp(x + y + z)* (3. + x *(4. + 3.* x) + y* (4. + 3.* y) + z* (4. + 3. * z));
  /* return -2*exp(x + y)*(x*(2. + x) + 1. + y*(2.+y)); */
}

static void 
d4est_test_iterate_refine_only_one_element
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  d4est_amr_t* d4est_amr = (d4est_amr_t*) info->p4est->user_pointer;
  int* element_to_refine = (int*) (d4est_amr->scheme->amr_scheme_data);
  d4est_element_data_t* elem_data = (d4est_element_data_t*) info->quad->p.user_data;
  if (elem_data->id == *element_to_refine){
    d4est_amr->refinement_log[elem_data->id] = -elem_data->deg;
  }
  else {
    d4est_amr->refinement_log[elem_data->id] = elem_data->deg;
  }
}

int main(int argc, char *argv[])
{
  sc_MPI_Comm mpicomm;
  PetscInitialize(&argc,&argv,(char*)0,NULL);
  mpicomm = PETSC_COMM_WORLD;
  
  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  
#ifndef NDEBUG
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DEBUG MODE ON\n");
  p4est_init(NULL, SC_LP_ERROR);
  /* p4est_init(NULL, SC_LP_ALWAYS); */
#else
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DEBUG MODE OFF\n");
  p4est_init(NULL, SC_LP_ERROR);
#endif
  
#if (P4EST_DIM)==3
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DIM = 3\n");
#else
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DIM = 2\n");
#endif
  const char* default_input_file = "d4est_test_schwarz_cubed_sphere.input";
  
  if (proc_rank == 0)
    printf("[D4EST_INFO]: options file = %s\n", (argc == 2) ? argv[1] : default_input_file);
 
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    (argc == 2) ? argv[1] : default_input_file,
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse((argc == 2) ? argv[1] :      (argc == 2) ? argv[1] : default_input_file, d4est_geom);

  p4est_t* p4est;
  p4est = p4est_new_ext
          (
           mpicomm,
           d4est_geom->p4est_conn,
           initial_grid_input->min_quadrants,
           initial_grid_input->min_level,
           initial_grid_input->fill_uniform,
           sizeof(d4est_element_data_t),
           NULL,
           NULL
          );
  

  printf("num_quadrants = %d\n", p4est->local_num_quadrants);
  p4est_partition(p4est, 1, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);  
  d4est_ghost_t* d4est_ghost = NULL;   

  if (proc_rank == 0){
    printf("[D4EST_INFO]: mpisize = %d\n", proc_size);
  }
  if (proc_rank == 0 && initial_grid_input->load_from_checkpoint == 0){
    printf("[D4EST_INFO]: min_quadrants = %d\n", initial_grid_input->min_quadrants);
    printf("[D4EST_INFO]: min_level = %d\n", initial_grid_input->min_level);
    printf("[D4EST_INFO]: fill_uniform = %d\n", initial_grid_input->fill_uniform);
  }
  
  sc_MPI_Barrier(mpicomm);
  printf("[D4EST_INFO]: elements on proc %d = %d\n", proc_rank, p4est->local_num_quadrants);
  sc_MPI_Barrier(mpicomm);
  
  /* start just-in-time dg-math */
  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_data_t* d4est_factors = d4est_mesh_data_init(p4est);
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, (argc == 2) ? argv[1] : default_input_file, "quadrature");

  /* BEGIN AMR AND MESH STUFF */
  /* BEGIN AMR AND MESH STUFF */
  /* BEGIN AMR AND MESH STUFF */
  /* BEGIN AMR AND MESH STUFF */
  /* BEGIN AMR AND MESH STUFF */
  
  d4est_mesh_local_sizes_t local_sizes= d4est_mesh_update
                                        (
                                         p4est,
                                         &d4est_ghost,
                                         d4est_ops,
                                         d4est_geom,
                                         d4est_quad,
                                         d4est_factors,
                                         initial_grid_input,
                                         INITIALIZE_GHOST,
                                         INITIALIZE_QUADRATURE_DATA,
                                         INITIALIZE_GEOMETRY_DATA,
                                         INITIALIZE_GEOMETRY_ALIASES,
                                         d4est_mesh_set_initial_extents,
                                         (void*)initial_grid_input
                                        );



  /* create amr scheme */
  int add_amr_in_one_element = 0;
  if (add_amr_in_one_element){
    int* refinement_log = P4EST_ALLOC(int, p4est->local_num_quadrants);
    int element_to_refine = 0;
    d4est_amr_t* d4est_amr = P4EST_ALLOC(d4est_amr_t, 1);
    d4est_amr_scheme_t* scheme = P4EST_ALLOC(d4est_amr_scheme_t, 1);
    scheme->post_balance_callback = NULL;
    scheme->pre_refine_callback = NULL;
    scheme->refine_replace_callback_fcn_ptr = NULL;
    scheme->balance_replace_callback_fcn_ptr = NULL;
    scheme->mark_elements = d4est_test_iterate_refine_only_one_element;
    scheme->amr_scheme_data = &element_to_refine;
    scheme->destroy = NULL;
    
    d4est_amr->mpirank = p4est->mpirank;
    d4est_amr->scheme = scheme;
    d4est_amr->balance_log = NULL;
    d4est_amr->refinement_log = NULL;
    d4est_amr->initial_log = NULL;
    d4est_amr->max_degree = 1000;
    
    d4est_amr_step
      (
       p4est,
       NULL,
       d4est_amr,
       NULL,
       NULL,
       NULL
      );

    P4EST_FREE(refinement_log);
    P4EST_FREE(d4est_amr);
    P4EST_FREE(scheme);

    local_sizes = d4est_mesh_update
                  (
                   p4est,
                   &d4est_ghost,
                   d4est_ops,
                   d4est_geom,
                   d4est_quad,
                   d4est_factors,
                   initial_grid_input,
                   INITIALIZE_GHOST,
                   INITIALIZE_QUADRATURE_DATA,
                   INITIALIZE_GEOMETRY_DATA,
                   INITIALIZE_GEOMETRY_ALIASES,
                   d4est_mesh_set_quadratures_after_amr,
                   (void*)initial_grid_input
                  );
  
    initial_grid_input->initial_nodes = local_sizes.local_nodes;
  }

  d4est_field_type_t field_type = NODAL;
  d4est_ghost_data_t* d4est_ghost_data = d4est_ghost_data_init(p4est,
                                                               d4est_ghost,
                                                               &field_type,
                                                               1);

  
  /* END AMR AND MESH STUFF */
  /* END AMR AND MESH STUFF */
  /* END AMR AND MESH STUFF */
  /* END AMR AND MESH STUFF */
  /* END AMR AND MESH STUFF */
  
  d4est_solver_schwarz_data_t* schwarz_data
    = d4est_solver_schwarz_init(
                                p4est,
                                d4est_ghost,
                                (argc == 2) ? argv[1] : default_input_file
    );


  d4est_solver_schwarz_operators_t* schwarz_ops
    = d4est_solver_schwarz_operators_init
    (d4est_ops);

  d4est_solver_schwarz_print
    (
     p4est,
     schwarz_data
    );

  /* BEGIN SETUP EQUATIONS AND DATA */
  /* BEGIN SETUP EQUATIONS AND DATA */
  /* BEGIN SETUP EQUATIONS AND DATA */
  /* BEGIN SETUP EQUATIONS AND DATA */
  /* BEGIN SETUP EQUATIONS AND DATA */
  
  dirichlet_bndry_eval_method_t eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;

  d4est_laplacian_dirichlet_bc_t bc_data_for_residual;
  bc_data_for_residual.dirichlet_fcn = poly_vec_fcn;
  bc_data_for_residual.eval_method = eval_method;  

  d4est_laplacian_dirichlet_bc_t bc_data_for_lhs;
  bc_data_for_lhs.dirichlet_fcn = zero_fcn;
  bc_data_for_lhs.eval_method = eval_method;  
    
  d4est_laplacian_flux_data_t* flux_data_for_apply_lhs = d4est_laplacian_flux_new(p4est, (argc == 2) ? argv[1] : default_input_file, BC_DIRICHLET, &bc_data_for_lhs);

  d4est_laplacian_flux_data_t* flux_data_for_residual = d4est_laplacian_flux_new(p4est, (argc == 2) ? argv[1] : default_input_file, BC_DIRICHLET, &bc_data_for_residual);

  d4est_elliptic_eqns_t prob_fcns_for_lhs;
  prob_fcns_for_lhs.build_residual = d4est_test_build_residual;
  prob_fcns_for_lhs.apply_lhs = d4est_test_apply_lhs;
  prob_fcns_for_lhs.user = flux_data_for_apply_lhs;

  d4est_elliptic_eqns_t prob_fcns_for_residual;
  prob_fcns_for_residual.build_residual = d4est_test_build_residual;
  prob_fcns_for_residual.apply_lhs = d4est_test_apply_lhs;
  prob_fcns_for_residual.user = flux_data_for_residual;

    
  double* u = P4EST_ALLOC_ZERO(double, local_sizes.local_nodes);
  double* Au = P4EST_ALLOC_ZERO(double, local_sizes.local_nodes);
  double* rhs = P4EST_ALLOC(double, local_sizes.local_nodes);
  double* f = P4EST_ALLOC(double, local_sizes.local_nodes);
  double* r = P4EST_ALLOC(double, local_sizes.local_nodes);

  d4est_elliptic_data_t elliptic_data;
  elliptic_data.u = u;
  elliptic_data.Au = Au;
  elliptic_data.rhs = rhs;
  elliptic_data.local_nodes = local_sizes.local_nodes;
  elliptic_data.field_types = &field_type;
  elliptic_data.num_of_fields = 1;

  /* END SETUP EQUATIONS AND DATA */
  /* END SETUP EQUATIONS AND DATA */
  /* END SETUP EQUATIONS AND DATA */
  /* END SETUP EQUATIONS AND DATA */
  /* END SETUP EQUATIONS AND DATA */  

    
  d4est_mesh_init_field
    (
     p4est,
     u,
     poly_vec_fcn,
     d4est_ops, // unnecessary?
     d4est_geom, // unnecessary?
     d4est_factors,
     INIT_FIELD_ON_LOBATTO,
     NULL
    );
    

  double* u_over_subdomains = P4EST_ALLOC(double,
                                          schwarz_data->restricted_nodal_size);

  double* rhs_over_subdomains = P4EST_ALLOC(double,
                                            schwarz_data->restricted_nodal_size);



  double* zeroed_u = P4EST_ALLOC(double, local_sizes.local_nodes);
  double* zeroed_Au = P4EST_ALLOC(double, local_sizes.local_nodes);
  double* error = P4EST_ALLOC(double, local_sizes.local_nodes);

  double* sol_over_subdomains = P4EST_ALLOC_ZERO(double,
                                            schwarz_data->restricted_nodal_size);

  
  d4est_laplacian_build_rhs_with_strong_bc
    (
     p4est,
     d4est_ghost,
     d4est_ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &elliptic_data,
     flux_data_for_residual,
     rhs,
     neg_laplacian_poly_vec_fcn,
     INIT_FIELD_ON_LOBATTO,
     NULL,
     0
    );

  
  d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomains
    (
     p4est,
     d4est_ops,
     d4est_factors,
     d4est_ghost_data,
     schwarz_data,
     schwarz_ops,
     rhs,
     -1,
     rhs_over_subdomains
    );
 
  d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomains
    (
     p4est,
     d4est_ops,
     d4est_factors,
     d4est_ghost_data,
     schwarz_data,
     schwarz_ops,
     u,
     -1,
     u_over_subdomains
    );



  /* for (int i = 0; i < schwarz_data->num_subdomains; i++){ */

  /*   printf("stride = %d\n", schwarz_data->subdomain_data[i].restricted_nodal_stride); */
    
  /*   d4est_solver_schwarz_cg_solve_subdomain_single_core */
  /*     ( */
  /*      p4est, */
  /*      d4est_geom, */
  /*      d4est_quad, */
  /*      d4est_factors, */
  /*      d4est_ghost, */
  /*      d4est_ghost_data, */
  /*      schwarz_data, */
  /*      schwarz_ops, */
  /*      &elliptic_data, */
  /*      &prob_fcns, */
  /*      &sol_over_subdomains[schwarz_data->subdomain_data[i].restricted_nodal_stride], */
  /*      &rhs_over_subdomains[schwarz_data->subdomain_data[i].restricted_nodal_stride], */
  /*      zeroed_u, */
  /*      zeroed_Au, */
  /*      1000, */
  /*      1e-15, */
  /*      1e-15, */
  /*      i */
  /*     ); */
  /* } */

  /* DEBUG_PRINT_2ARR_DBL(rhs_over_subdomains, sol_over_subdomains, schwarz_data->restricted_nodal_size); */

  d4est_vtk_helper_array_t* helper_array = d4est_vtk_helper_array_init
    (
     p4est,
     d4est_ops,
     local_sizes.local_nodes,
     25,
     500
    );
  
  int iterations = 20;
  double* sol = P4EST_ALLOC_ZERO(double, local_sizes.local_nodes);
  
  for (int i = 0; i < iterations; i++){

    double* pre_solve_Au = d4est_vtk_helper_array_alloc_and_add_nodal_dbl_field
      (
       helper_array,
       "pre_solve_Au",
       i
      );    
    
    elliptic_data.u = sol;
    elliptic_data.Au = pre_solve_Au;
    elliptic_data.rhs = rhs;


    d4est_elliptic_eqns_apply_lhs
      (
       p4est,
       d4est_ghost,
       d4est_ghost_data,
       &prob_fcns_for_lhs,
       &elliptic_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors
    );

    elliptic_data.Au = r;
    
    /* if (i == 0){ */
    d4est_elliptic_eqns_build_residual
      (
       p4est,
       d4est_ghost,
       d4est_ghost_data,
       &prob_fcns_for_lhs,
       &elliptic_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors
    );
    /* } */
    
    double r2 = d4est_linalg_vec_dot(r, r, local_sizes.local_nodes);  
    printf("r2 = %.15f\n", r2);


    double* pre_solve_r = d4est_vtk_helper_array_alloc_and_add_nodal_dbl_field
      (
       helper_array,
       "pre_solve_r",
       i
      );
    
    double* pre_solve_u = d4est_vtk_helper_array_alloc_and_add_nodal_dbl_field
      (
       helper_array,
       "pre_solve_u",
       i
      );

    d4est_util_copy_1st_to_2nd(sol, pre_solve_u
                               ,local_sizes.local_nodes);
    d4est_util_copy_1st_to_2nd(r, pre_solve_r
                               ,local_sizes.local_nodes);


    printf("Solving\n");
    d4est_solver_schwarz_compute_and_add_correction
      (
       p4est,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       d4est_ghost,
       d4est_ghost_data,
       schwarz_data,
       schwarz_ops,
       &prob_fcns_for_lhs,
       &elliptic_data,
       sol,
       r,
       1000,
       1e-15,
       1e-15,
       NULL,//helper_array,
       -1//i
      );


    double* post_solve_u = d4est_vtk_helper_array_alloc_and_add_nodal_dbl_field
      (
       helper_array,
       "post_solve_u",
       i
      );

    double* post_solve_r = d4est_vtk_helper_array_alloc_and_add_nodal_dbl_field
      (
       helper_array,
       "post_solve_r",
       i
      );

    d4est_util_copy_1st_to_2nd(sol, post_solve_u
                               ,local_sizes.local_nodes);
    d4est_util_copy_1st_to_2nd(r, post_solve_r
                               ,local_sizes.local_nodes);
    
    /* for (int n = 0; n < local_sizes.local_nodes; n++){ */
      /* sol[n] = zeroed_u[n] */
    /* } */
    
    d4est_util_compute_error_array(sol, u, error, local_sizes.local_nodes);
    double l2 = d4est_mesh_compute_l2_norm_sqr
                (
                 p4est,
                 d4est_ops,
                 d4est_geom,
                 d4est_quad,
                 d4est_factors,
                 error,
                 local_sizes.local_nodes,
                 NULL,
                 NULL);
    printf("l2 norm = %.15f\n", l2);
    
  }

  P4EST_FREE(sol);
    
  /* BEGIN COLLECT DATA TO VISUALIZE SUBDOMAINS AND WEIGHTS */
  /* BEGIN COLLECT DATA TO VISUALIZE SUBDOMAINS AND WEIGHTS */
  /* BEGIN COLLECT DATA TO VISUALIZE SUBDOMAINS AND WEIGHTS */
  /* BEGIN COLLECT DATA TO VISUALIZE SUBDOMAINS AND WEIGHTS */
  /* BEGIN COLLECT DATA TO VISUALIZE SUBDOMAINS AND WEIGHTS */


  for (int i = 0; i < p4est->local_num_quadrants; i++){

    double* d_sub = d4est_vtk_helper_array_alloc_and_add_nodal_dbl_field
      (
       helper_array,
       "d_sub",
       i
      );
       
    double* w_d_sub = d4est_vtk_helper_array_alloc_and_add_nodal_dbl_field
      (
       helper_array,
       "w_d_sub",
       i
      );

      int* i_sub = d4est_vtk_helper_array_alloc_and_add_cell_int_field
      (
       helper_array,
       "i_sub",
       i
      );
    
    d4est_solver_schwarz_visualize_subdomain_helper_single_core
      (
       p4est,
       d4est_factors,
       schwarz_data,
       schwarz_ops,
       i,
       i_sub,
       d_sub,
       w_d_sub
      );    
  }

  for (int i = 0; i < p4est->local_num_quadrants; i++){

    d4est_solver_schwarz_subdomain_data_t* sub_data = &schwarz_data->subdomain_data[i];

    double* u_sub = d4est_vtk_helper_array_alloc_and_add_nodal_dbl_field
                    (
                     helper_array,
                     "u_sub",
                     i
                    );
      
    double* rhs_sub = d4est_vtk_helper_array_alloc_and_add_nodal_dbl_field
                      (
                       helper_array,
                       "rhs_sub",
                       i
                      );
       
    d4est_solver_schwarz_convert_restricted_subdomain_field_to_global_nodal_field
      (
       p4est,
       d4est_factors,
       schwarz_data,
       schwarz_ops,
       &u_over_subdomains[sub_data->restricted_nodal_stride],
       u_sub,
       p4est->mpirank,
       i,
       local_sizes.local_nodes,
       FIELD_NOT_ZEROED
      );

    d4est_solver_schwarz_convert_restricted_subdomain_field_to_global_nodal_field
      (
       p4est,
       d4est_factors,
       schwarz_data,
       schwarz_ops,
       &rhs_over_subdomains[sub_data->restricted_nodal_stride],
       rhs_sub,
       p4est->mpirank,
       i,
       local_sizes.local_nodes,
       FIELD_NOT_ZEROED
      );
    
    /* d4est_solver_schwarz_convert_restricted_subdomain_field_to_global_nodal_field */
    /*   ( */
    /*    p4est, */
    /*    d4est_factors, */
    /*    schwarz_data, */
    /*    schwarz_ops, */
    /*    &sol_over_subdomains[sub_data->restricted_nodal_stride], */
    /*    sol_sub, */
    /*    p4est->mpirank, */
    /*    i, */
    /*    local_sizes.local_nodes, */
    /*    FIELD_NOT_ZEROED */
    /*   ); */
    
  }     
  
  d4est_vtk_helper_array_add_nodal_dbl_field
    (
     helper_array,
     "u",
     -1,
     u
    );


  d4est_vtk_helper_array_add_nodal_dbl_field
    (
     helper_array,
     "rhs",
     -1,
     rhs
    );
  
  d4est_vtk_save_helper_array
    (
     helper_array,
     (argc == 2) ? argv[1] : (char*)default_input_file
    );


  d4est_vtk_helper_array_destroy
    (
     helper_array
    );
  /* BEGIN FREE ELLIPTIC DATA */
  /* BEGIN FREE ELLIPTIC DATA */
  /* BEGIN FREE ELLIPTIC DATA */
  /* BEGIN FREE ELLIPTIC DATA */
    
  P4EST_FREE(u);
  P4EST_FREE(Au);
  P4EST_FREE(error);
  P4EST_FREE(rhs);
  P4EST_FREE(f);
  P4EST_FREE(zeroed_u);
  P4EST_FREE(zeroed_Au);
  P4EST_FREE(r);
    

  /* BEGIN END ELLIPTIC DATA */
  /* BEGIN END ELLIPTIC DATA */
  /* BEGIN END ELLIPTIC DATA */
  /* BEGIN END ELLIPTIC DATA */

  /* BEGIN FREE SUBDOMAIN DATA */
  /* BEGIN FREE SUBDOMAIN DATA */
  /* BEGIN FREE SUBDOMAIN DATA */
  /* BEGIN FREE SUBDOMAIN DATA */
  
  P4EST_FREE(sol_over_subdomains);
  P4EST_FREE(u_over_subdomains);
  P4EST_FREE(rhs_over_subdomains);


  /* BEGIN END SUBDOMAIN DATA */
  /* BEGIN END SUBDOMAIN DATA */
  /* BEGIN END SUBDOMAIN DATA */
  /* BEGIN END SUBDOMAIN DATA */


  if (d4est_ghost_data != NULL){
    d4est_ghost_data_destroy(d4est_ghost_data);
    d4est_ghost_data = NULL;
  } 
  
  d4est_solver_schwarz_destroy
    (
     schwarz_data
    );

  d4est_solver_schwarz_operators_destroy
    (
     schwarz_ops
    );
  
  d4est_mesh_initial_extents_destroy(initial_grid_input);
  d4est_mesh_data_destroy(d4est_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  
  if (d4est_ghost) {
    d4est_ghost_destroy(d4est_ghost);
  }

  d4est_laplacian_flux_destroy(flux_data_for_apply_lhs);
  d4est_laplacian_flux_destroy(flux_data_for_residual);
  
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);
  PetscFinalize();
  return 0;
}

