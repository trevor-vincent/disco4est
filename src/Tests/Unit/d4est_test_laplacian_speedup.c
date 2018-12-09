#include <pXest.h>
#include <d4est_quadrature.h>
#include <d4est_element_data.h>
#include <d4est_quadrature_legendre.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_geometry_brick.h>
#include <petscsnes.h>
#include <d4est_linalg.h>
#include <d4est_mortars.h>
#include <d4est_amr.h>
#include <d4est_laplacian_with_opt.h>
#include <d4est_laplacian.h>
#include <d4est_solver_matrix_symmetry.h>
#include <d4est_util.h>
#include <d4est_norms.h>
#include <d4est_vtk.h>
#include <limits.h>
#include <zlog.h>
#include <time.h>

#define D4EST_REAL_EPS 100*1e-15
#if (P4EST_DIM)==2
#define TEST_DEG_INIT 2
#else
#define TEST_DEG_INIT 2
#endif

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


static void
d4est_test_apply_lhs_with_opt
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
  d4est_laplacian_with_opt_flux_data_t* flux_data_for_apply_lhs = user;  
  d4est_laplacian_with_opt_apply_aij(p4est, ghost, ghost_data, prob_vecs, flux_data_for_apply_lhs, d4est_ops, d4est_geom, d4est_quad, d4est_factors, 0);
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

#if (P4EST_DIM)==3
  return x*x + y*y + z*z;
#else
  return x*x + y*y;
#endif
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
#if (P4EST_DIM)==3
  return -6.;
#else
  return -4.;
#endif
}


static void
problem_set_degrees_amr
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  elem_data->deg_quad = elem_data->deg;

}



int main(int argc, char *argv[])
{

#ifndef D4EST_TEST
  D4EST_ABORT("D4EST_TEST not defined");
#endif
  
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

  if (proc_rank == 0)
    printf("[D4EST_INFO]: options file = %s\n", (argc == 2) ? argv[1] :      "d4est_test_laplacian_speedup.input");
 
 
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    (argc == 2) ? argv[1] :      "d4est_test_laplacian_speedup.input",
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse((argc == 2) ? argv[1] :      "d4est_test_laplacian_speedup.input", d4est_geom);

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


  p4est_partition(p4est, 1, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  
  d4est_ghost_t* d4est_ghost = d4est_ghost_init(p4est);
     
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
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, (argc == 2) ? argv[1] :      "d4est_test_laplacian_speedup.input", "quadrature");
  

  d4est_mesh_local_sizes_t local_sizes = d4est_mesh_update
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

  initial_grid_input->initial_nodes = local_sizes.local_nodes;

  p4est_partition(p4est, 1, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  int same2 = 0;
  int same = 0;
  
  /* d4est_amr_t* d4est_amr_random = d4est_amr_init_uniform_h(p4est, 7, 1); */
  int num_of_amr_steps = 2;
  d4est_amr_t* d4est_amr_random = d4est_amr_init_random_hp(p4est, num_of_amr_steps);
  int bad_element = 0;
  for (int i = 0; i < num_of_amr_steps; i++){

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
    else {
      d4est_amr_step
        (
         p4est,
         d4est_ops,
         d4est_amr_random,
         NULL,
         NULL,
         NULL
        );
    }    
    d4est_mesh_local_sizes_t local_sizes = d4est_mesh_update
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

    d4est_field_type_t field_type = NODAL;
    d4est_ghost_data_t* d4est_ghost_data = d4est_ghost_data_init(p4est,
                                                                 d4est_ghost,
                                                                 &field_type,
                                                                 1);

    
    int local_nodes = local_sizes.local_nodes;
    
    printf("level = %d, elements = %d, nodes = %d\n", i, p4est->local_num_quadrants, local_nodes);


    /* d4est_elliptic_data_t prob_vecs; */
    /* prob_vecs.Au = P4EST_ALLOC(double, local_nodes); */
    /* prob_vecs.u = P4EST_ALLOC(double, local_nodes); */
    /* /\* prob_vecs.rhs = P4EST_ALLOC(double, local_nodes); *\/ */
    /* prob_vecs.local_nodes = local_nodes; */
  

        
    dirichlet_bndry_eval_method_t eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;

    /* / Setup boundary conditions */
    d4est_laplacian_dirichlet_bc_t bc_data_for_lhs;
    bc_data_for_lhs.dirichlet_fcn = poly_vec_fcn;
    bc_data_for_lhs.eval_method = eval_method;  

    
    d4est_laplacian_flux_data_t* flux_data_for_apply_lhs = d4est_laplacian_flux_new(p4est, "d4est_test_laplacian_speedup.input", BC_DIRICHLET, &bc_data_for_lhs);
    d4est_laplacian_with_opt_flux_data_t* flux_data_for_apply_lhs_with_opt = d4est_laplacian_with_opt_flux_new(p4est, "d4est_test_laplacian_speedup.input", BC_DIRICHLET, &bc_data_for_lhs);


    d4est_elliptic_eqns_t prob_fcns;
    prob_fcns.build_residual = NULL;
    prob_fcns.apply_lhs = d4est_test_apply_lhs;
    prob_fcns.user = flux_data_for_apply_lhs;


    d4est_elliptic_eqns_t prob_fcns_with_opt;
    prob_fcns_with_opt.build_residual = NULL;
    prob_fcns_with_opt.apply_lhs = d4est_test_apply_lhs_with_opt;
    prob_fcns_with_opt.user = flux_data_for_apply_lhs_with_opt;


    
    double* poly_vec = P4EST_ALLOC(double, local_nodes);


    d4est_mesh_init_field(
                          p4est,
                          poly_vec,
                          poly_vec_fcn,
                          d4est_ops, // unnecessary?
                          d4est_geom, // unnecessary?
                          d4est_factors,
                          INIT_FIELD_ON_LOBATTO,
                          NULL
    );
    
    double* Apoly_vec = P4EST_ALLOC_ZERO(double, local_nodes);
    double* Apoly_vec_with_opt = P4EST_ALLOC_ZERO(double, local_nodes);

    double* Abc_poly_vec = P4EST_ALLOC(double, local_nodes);
    double* Apoly_vec_compare = P4EST_ALLOC(double, local_nodes);
    double* tmp = P4EST_ALLOC(double, local_nodes);
    d4est_elliptic_data_t elliptic_data;
    elliptic_data.u = poly_vec;
    elliptic_data.Au = Apoly_vec;
    elliptic_data.local_nodes = local_nodes;

    elliptic_data.field_types = &field_type;
    elliptic_data.num_of_fields = 1;

    srand(102321);
    for (int i = 0; i < local_sizes.local_nodes; i++){
      poly_vec[i] += (double)rand() / (double)RAND_MAX;
    }

    clock_t start = clock();
    d4est_elliptic_eqns_apply_lhs
      (
       p4est,
       d4est_ghost,
       d4est_ghost_data,
       &prob_fcns,
       &elliptic_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors
      );
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    elliptic_data.Au = Apoly_vec_with_opt;
    start = clock();
    d4est_elliptic_eqns_apply_lhs
      (
       p4est,
       d4est_ghost,
       d4est_ghost_data,
       &prob_fcns_with_opt,
       &elliptic_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors
      );
    end = clock();    
    double time_spent_2 = (double)(end - start) / CLOCKS_PER_SEC;    

    printf("**********DONE**********\n");
    printf("before and after opt times = %f, %f\n", time_spent, time_spent_2);

  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;

      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);
        
      int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);

      for (int j = 0; j < volume_nodes; j++){
        if((fabs(Apoly_vec[ed->nodal_stride + j] - Apoly_vec_with_opt[ed->nodal_stride + j]) > 1e-13)){
          int on_bndry = d4est_factors->node_touches_boundary[ed->nodal_stride + j];
          
          printf("ed id %d, on_bndy %d, %d, Apoly_vec[j], Apoly_vec_with_opt[j], err = %.15f, %.15f, %.15f\n",ed->id, on_bndry, i,Apoly_vec[ed->nodal_stride + j],Apoly_vec_with_opt[ed->nodal_stride + j]);
          bad_element = 1;
        }
      }
      if(bad_element == 1){
        printf("ed->deg = %d\n", ed->deg);
        double*vec_opt = &Apoly_vec_with_opt[ed->nodal_stride];
        DEBUG_PRINT_ARR_DBL(vec_opt, volume_nodes);
        D4EST_ABORT("ABORTING BAD ELEMENT");
      }
    }
  }
    double u_sum = d4est_util_sum_array_dbl(poly_vec, local_nodes);
    double Au_sum = d4est_util_sum_array_dbl(Apoly_vec, local_nodes);
    double Au_sum_with_opt = d4est_util_sum_array_dbl(Apoly_vec_with_opt, local_nodes);

    printf("local_nodes = %d\n", local_nodes);
    printf("quadrants = %d\n", p4est->local_num_quadrants);
    printf("u sum = %.25f\n", u_sum);
    printf("Au sum = %.25f\n", Au_sum);
    printf("Au sum with opt = %.25f\n", Au_sum_with_opt);

    P4EST_FREE(Apoly_vec);
    P4EST_FREE(Apoly_vec_with_opt);
    P4EST_FREE(poly_vec);
    P4EST_FREE(Abc_poly_vec);
    P4EST_FREE(Apoly_vec_compare);
    P4EST_FREE(tmp);

    if (d4est_ghost_data != NULL){
      d4est_ghost_data_destroy(d4est_ghost_data);
      d4est_ghost_data = NULL;
    } 

    d4est_laplacian_flux_destroy(flux_data_for_apply_lhs);
    d4est_laplacian_with_opt_flux_destroy(flux_data_for_apply_lhs_with_opt);
  }

  d4est_ghost_destroy(d4est_ghost);
  d4est_mesh_initial_extents_destroy(initial_grid_input);

  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_amr_destroy(d4est_amr_random);
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);
  
  PetscFinalize();
  /* if (same && same2) */
  return bad_element;
  /* else */
    /* return 1; */
}
