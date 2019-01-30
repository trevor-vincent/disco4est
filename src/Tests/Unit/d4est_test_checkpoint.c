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
#include <d4est_checkpoint.h>
#include <d4est_laplacian_flux.h>
#include <d4est_solver_matrix_symmetry.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_solver_krylov_petsc.h>
#include <d4est_util.h>
#include <d4est_norms.h>
#include <limits.h>
#include <zlog.h>
#include "../../Problems/Poisson/poisson_lorentzian_fcns.h"


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

  /* if (proc_rank == 0) */
    /* printf("[D4EST_INFO]: options file = %s\n", (argc == 2) ? argv[1] :      "d4est_test_laplacian_symmetry.input"); */
 


 if (zlog_init("logging.conf") != 0){
    printf("Initializing logging failed.\n");
    D4EST_ABORT("");
  }

  const char* input_file = (argc == 2) ? argv[1] : "d4est_test_checkpoint.input";
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    input_file,
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse(input_file, d4est_geom);

  p4est_t* p4est = d4est_checkpoint_load_p4est_from_file
                   (
                    mpicomm,
                    initial_grid_input,
                    &d4est_geom->p4est_conn
                   );

  p4est_reset_data (p4est,
                    sizeof(d4est_element_data_t),
                    NULL,
                    NULL);

  int* checkpoint_deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
  
  d4est_checkpoint_read_dataset
    (
     p4est,
     initial_grid_input->checkpoint_prefix,
     "degree",
     H5T_NATIVE_INT,
     checkpoint_deg_array,
     initial_grid_input->checkpoint_number
    );

  initial_grid_input->checkpoint_deg_array = checkpoint_deg_array;

  
  d4est_ghost_t* d4est_ghost = d4est_ghost_init(p4est);     
  D4EST_ASSERT(initial_grid_input->load_from_checkpoint == 1);
  D4EST_ASSERT(initial_grid_input->checkpoint_deg_array != NULL);
  
  /* start just-in-time dg-math */
  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_data_t* d4est_factors = d4est_mesh_data_init(p4est);
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, input_file, "quadrature");
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
  
  d4est_field_type_t field_type = NODAL;  
  d4est_ghost_data_t* d4est_ghost_data = d4est_ghost_data_init(p4est,
                                           d4est_ghost,
                                           &field_type,
                                           1);


  initial_grid_input->initial_nodes = local_sizes.local_nodes;
  int initial_nodes = local_sizes.local_nodes;
  
  d4est_solver_krylov_petsc_params_t d4est_solver_krylov_petsc_params;
  d4est_solver_krylov_petsc_input(p4est, input_file, "d4est_solver_krylov_petsc", &d4est_solver_krylov_petsc_params);

 dirichlet_bndry_eval_method_t eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;

  lorentzian_params_t lorentzian_params;
  lorentzian_params.R_surface = ((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user)->R2;

  d4est_laplacian_dirichlet_bc_t bc_data_dirichlet_for_lhs;
  bc_data_dirichlet_for_lhs.dirichlet_fcn = zero_fcn;
  bc_data_dirichlet_for_lhs.eval_method = eval_method;
  bc_data_dirichlet_for_lhs.user = &lorentzian_params;
  
  d4est_laplacian_dirichlet_bc_t bc_data_dirichlet_for_rhs;
  bc_data_dirichlet_for_rhs.dirichlet_fcn = poisson_lorentzian_boundary_fcn;
  bc_data_dirichlet_for_rhs.eval_method = eval_method;
  bc_data_dirichlet_for_rhs.user = &lorentzian_params;
  
  d4est_laplacian_flux_data_t* flux_data_for_lhs = NULL;
  d4est_laplacian_flux_data_t* flux_data_for_rhs = NULL; 

  flux_data_for_lhs
    = d4est_laplacian_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_dirichlet_for_lhs);
  
  flux_data_for_rhs
    = d4est_laplacian_flux_new(p4est, input_file,  BC_DIRICHLET, &bc_data_dirichlet_for_rhs);

  problem_ctx_t ctx;
  ctx.flux_data_for_apply_lhs = flux_data_for_lhs;
  ctx.flux_data_for_build_rhs = flux_data_for_rhs;

  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = poisson_lorentzian_build_residual;
  prob_fcns.apply_lhs = poisson_lorentzian_apply_lhs;
  prob_fcns.user = &ctx;

  d4est_elliptic_data_t prob_vecs;
  prob_vecs.Au = P4EST_ALLOC_ZERO(double, initial_nodes);
  prob_vecs.u = P4EST_ALLOC_ZERO(double, initial_nodes);
  prob_vecs.rhs = P4EST_ALLOC_ZERO(double, initial_nodes);
  prob_vecs.local_nodes = initial_nodes;
  prob_vecs.field_types = &field_type;
  prob_vecs.num_of_fields = 1;
  prob_vecs.local_nodes = local_sizes.local_nodes;
 
  d4est_checkpoint_read_dataset
    (
     p4est,
     initial_grid_input->checkpoint_prefix,
     "u",
     H5T_NATIVE_DOUBLE,
     prob_vecs.u,
     initial_grid_input->checkpoint_number
    );

    double sum = d4est_util_sum_array_dbl(prob_vecs.u, prob_vecs.local_nodes);
    d4est_checkpoint_check_dataset(p4est,
                           initial_grid_input->checkpoint_prefix,
                           "u",
                           H5T_NATIVE_DOUBLE,
                           (void*)&sum,
                           initial_grid_input->checkpoint_number
                          );
  
  d4est_laplacian_build_rhs_with_strong_bc
    (
     p4est,
     d4est_ghost,
     d4est_ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &prob_vecs,
     flux_data_for_rhs,
     prob_vecs.rhs,
     poisson_lorentzian_rhs_fcn,
     INIT_FIELD_ON_LOBATTO,
     &ctx,
     0
    );
  
  d4est_solver_krylov_petsc_solve
    (
     p4est,
     &prob_vecs,
     &prob_fcns,
     &d4est_ghost,
     &d4est_ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &d4est_solver_krylov_petsc_params,
     NULL,
     12
    );
    
  d4est_norms_fcn_L2_ctx_t L2_norm_ctx;
  L2_norm_ctx.p4est = p4est;
  L2_norm_ctx.d4est_ops = d4est_ops;
  L2_norm_ctx.d4est_geom = d4est_geom;
  L2_norm_ctx.d4est_quad = d4est_quad;
  L2_norm_ctx.d4est_factors = d4est_factors;
  
  d4est_norms_save(
                   p4est,
                   d4est_factors,
                   (const char * []){ "u", NULL },
                   (double * []){ prob_vecs.u },
                   (double * []){ NULL },
                   (d4est_xyz_fcn_t []){ poisson_lorentzian_analytic_solution },
                   (void * []){ &ctx },
                   (const char * []){"L_2", "L_infty", NULL},
                   (d4est_norm_fcn_t[]){ &d4est_norms_fcn_L2, &d4est_norms_fcn_Linfty },
                   (void * []){ &L2_norm_ctx, NULL},
                   NULL,
                   NULL,
                   NULL
  );

  P4EST_FREE(prob_vecs.u);
  P4EST_FREE(prob_vecs.Au);
  P4EST_FREE(prob_vecs.rhs);


  if (d4est_ghost_data != NULL){
    d4est_ghost_data_destroy(d4est_ghost_data);
    d4est_ghost_data = NULL;
  }
  
  zlog_fini();
  printf("Nodes = %d\n", local_sizes.local_nodes);
  printf("Elements = %d\n", p4est->local_num_quadrants);
  P4EST_FREE(checkpoint_deg_array);
  return 0;
}
/*  */
