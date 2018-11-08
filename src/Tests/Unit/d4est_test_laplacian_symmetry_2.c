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
#include <d4est_laplacian.h>
#include <d4est_laplacian_flux.h>
#include <d4est_solver_matrix_symmetry.h>
#include <d4est_util.h>
#include <d4est_norms.h>
#include <limits.h>
#include <zlog.h>


#define D4EST_REAL_EPS 100*1e-5
#if (P4EST_DIM)==2
#define TEST_DEG_INIT 2
#else
#define TEST_DEG_INIT 2
#endif


static void
d4est_test_poisson_symmetry_apply_lhs
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
problem_set_degrees_amr
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  elem_data->deg_quad = elem_data->deg;

}


static p4est_t*
problem_build_p4est
(
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t* conn,
 p4est_locidx_t min_quadrants,
 int min_level,
 int fill_uniform
)
{
  return p4est_new_ext
    (
     mpicomm,
     conn,
     min_quadrants,
     min_level,
     fill_uniform,
     sizeof(d4est_element_data_t),
     NULL,
     NULL
    );
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
    printf("[D4EST_INFO]: options file = %s\n", (argc == 2) ? argv[1] :      "d4est_test_laplacian_symmetry_2.input");
 
 
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    (argc == 2) ? argv[1] :      "d4est_test_laplacian_symmetry_2.input",
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse((argc == 2) ? argv[1] :      "d4est_test_laplacian_symmetry_2.input", d4est_geom);

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
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, (argc == 2) ? argv[1] :      "d4est_test_laplacian_symmetry_2.input", "quadrature");
  


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

  int num_of_amr_steps = 1;
  d4est_amr_t* d4est_amr_random = d4est_amr_init_random_hp(p4est, num_of_amr_steps);

  int same = 0;
  int same2 = 0;
  int nodes = -1;
  for (int i = 0; i < num_of_amr_steps; i++){

    
    /* d4est_amr_step */
    /*   ( */
    /*    p4est, */
    /*    d4est_ops, */
    /*    d4est_amr_random, */
    /*    NULL, */
    /*    NULL, */
    /*    NULL */
    /*   ); */
    

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


    d4est_elliptic_data_t prob_vecs;
    prob_vecs.Au = P4EST_ALLOC(double, local_nodes);
    prob_vecs.u = P4EST_ALLOC(double, local_nodes);
    /* prob_vecs.rhs = P4EST_ALLOC(double, local_nodes); */
    prob_vecs.local_nodes = local_nodes;
  

    prob_vecs.field_types = &field_type;
    prob_vecs.num_of_fields = 1;
    
    dirichlet_bndry_eval_method_t eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;

    /* / Setup boundary conditions */
    d4est_laplacian_dirichlet_bc_t bc_data_for_lhs;
    bc_data_for_lhs.dirichlet_fcn = zero_fcn;
    bc_data_for_lhs.eval_method = eval_method;  

    
    d4est_laplacian_flux_data_t* flux_data_for_apply_lhs = d4est_laplacian_flux_new(p4est, "d4est_test_laplacian_symmetry_2.input", BC_DIRICHLET, &bc_data_for_lhs);


    d4est_elliptic_eqns_t prob_fcns;
    prob_fcns.build_residual = NULL;
    prob_fcns.apply_lhs = d4est_test_poisson_symmetry_apply_lhs;
    prob_fcns.user = flux_data_for_apply_lhs;

    d4est_mesh_init_field(
      p4est,
      prob_vecs.u,
      sinpix_fcn,
      d4est_ops, // unnecessary?
      d4est_geom, // unnecessary?
      d4est_factors,
      INIT_FIELD_ON_LOBATTO,
      NULL
    );

    d4est_elliptic_eqns_apply_lhs
      (
       p4est,
       d4est_ghost,
       d4est_ghost_data,
       &prob_fcns,
       &prob_vecs,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors
      );

    printf("Au[2] = %.15f\n", prob_vecs.Au[2]);    

    same2 = fabs(prob_vecs.Au[2] - 0.34630989) < (D4EST_REAL_EPS);
    
    same = d4est_solver_matrix_symmetry
               (
                p4est,
                d4est_ghost,
                d4est_ghost_data,
                &prob_fcns,
                &prob_vecs,
                d4est_ops,
                d4est_geom,
                d4est_quad,
                d4est_factors,
                SYM_PRINT_MAT,
                D4EST_REAL_EPS
               );



    
    P4EST_FREE(prob_vecs.Au);
    P4EST_FREE(prob_vecs.u);
    /* P4EST_FREE(prob_vecs.rhs); */


    if (d4est_ghost_data != NULL){
      d4est_ghost_data_destroy(d4est_ghost_data);
      d4est_ghost_data = NULL;
    } 

  d4est_laplacian_flux_destroy(flux_data_for_apply_lhs);
    
  }  


  d4est_ghost_destroy(d4est_ghost);

  
  d4est_mesh_initial_extents_destroy(initial_grid_input);

  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_amr_destroy(d4est_amr_random);
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);
  
  PetscFinalize();
  /* sc_finalize (); */
  if (same && same2)
    return 0;
  else
    return 1;
}
/*  */
