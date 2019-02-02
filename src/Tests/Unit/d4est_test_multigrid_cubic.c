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
#include <d4est_hessian.h>
#include <d4est_laplacian_flux.h>
#include <d4est_solver_matrix_symmetry.h>
#include <d4est_util.h>
#include <d4est_norms.h>
#include <d4est_vtk.h>
#include <sc_reduce.h>
#include <limits.h>
#include <zlog.h>
#include <ini.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_solver_schwarz_transfer_ghost_data.h>
#include <d4est_solver_schwarz_helpers.h>
#include <d4est_solver_schwarz_geometric_data.h>
#include <d4est_solver_schwarz_laplacian_ext.h>
#include <d4est_solver_schwarz.h>
#include <d4est_solver_multigrid.h>

static void
d4est_test_schwarz_apply_lhs
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* d4est_ghost,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data,
 int subdomain,
 double* u_restricted_field_over_subdomain,
 double* Au_restricted_field_over_subdomain,
 void* ctx
){

  d4est_solver_schwarz_laplacian_ext_apply_over_subdomain
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     d4est_ghost,
     schwarz_data,
     schwarz_ops,
     schwarz_geometric_data,
     ctx,
     u_restricted_field_over_subdomain,
     Au_restricted_field_over_subdomain,
     subdomain
    );
}


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

  double pi = 3.1415926535897932384626433832795;
#if (P4EST_DIM)==3
  /* return (x-1.)*x*(y-1.)*y*(z-1)*z; */
  return sin(pi*x)*sin(pi*y)*sin(pi*z);
#else
  /* return (x-1.)*x*(y-1.)*y; */
  return sin(pi*x)*sin(pi*y);
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
  double pi = 3.1415926535897932384626433832795;
#if (P4EST_DIM)==3
  /* return - (2 * ((-1 + x) * x * (-1 + y) * y + (-1 + x) * x * (-1 + z) * z + (-1 + y) * y * (-1 + z)*z)); */
  return  3.*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z);
#else
  return  2.*pi*pi*sin(pi*x)*sin(pi*y);
  /* return - (2.*(x - 1.)*x + 2.*(y - 1.)*y); */
#endif
}

int main(int argc, char *argv[])
{

#ifndef D4EST_TEST
  D4EST_ABORT("D4EST_TEST not defined");
#endif
  
  sc_MPI_Comm mpicomm;
  PetscInitialize(&argc,&argv,(char*)0,NULL);
  mpicomm = PETSC_COMM_WORLD;

 if (zlog_init("logging.conf") != 0){
    printf("Initializing logging failed.\n");
    D4EST_ABORT("");
  }
  
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
  /*  */
  char* input_file = P4EST_ALLOC(char, 100);
  sprintf(input_file, "%s", (argc == 2) ? argv[1] : "d4est_test_multigrid_cubic.input");
  
  if (proc_rank == 0)
    printf("[D4EST_INFO]: options file = %s\n", input_file);
    
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    input_file,
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse(input_file, d4est_geom);

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
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, input_file, "quadrature");
  


  
  d4est_ghost_t* d4est_ghost = NULL;
  
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

  dirichlet_bndry_eval_method_t eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;

  d4est_laplacian_dirichlet_bc_t bc_data_for_residual;
  bc_data_for_residual.dirichlet_fcn = poly_vec_fcn;
  bc_data_for_residual.eval_method = eval_method;  

  d4est_laplacian_dirichlet_bc_t bc_data_for_lhs;
  bc_data_for_lhs.dirichlet_fcn = zero_fcn;
  bc_data_for_lhs.eval_method = eval_method;  
    
  d4est_laplacian_flux_data_t* flux_data_for_apply_lhs = d4est_laplacian_flux_new(p4est, (argc == 2) ? argv[1] : input_file, BC_DIRICHLET, &bc_data_for_lhs);

  d4est_laplacian_flux_data_t* flux_data_for_residual = d4est_laplacian_flux_new(p4est, (argc == 2) ? argv[1] : input_file, BC_DIRICHLET, &bc_data_for_residual);

  
  d4est_elliptic_eqns_t prob_fcns_for_lhs;
  prob_fcns_for_lhs.build_residual = d4est_test_build_residual;
  prob_fcns_for_lhs.apply_lhs = d4est_test_apply_lhs;
  prob_fcns_for_lhs.user = flux_data_for_apply_lhs;

  d4est_elliptic_eqns_t prob_fcns_for_residual;
  prob_fcns_for_residual.build_residual = d4est_test_build_residual;
  prob_fcns_for_residual.apply_lhs = d4est_test_apply_lhs;
  prob_fcns_for_residual.user = flux_data_for_residual;

    
  double* u = P4EST_ALLOC_ZERO(double, local_sizes.local_nodes);
  double* rhs = P4EST_ALLOC(double, local_sizes.local_nodes);
  double* r = P4EST_ALLOC(double, local_sizes.local_nodes);
  double* sol = P4EST_ALLOC(double, local_sizes.local_nodes);
  double* error = P4EST_ALLOC(double, local_sizes.local_nodes);
  d4est_field_type_t field_type = NODAL;
  d4est_elliptic_data_t elliptic_data;
  elliptic_data.u = u;
  elliptic_data.Au = r;
  elliptic_data.rhs = rhs;
  elliptic_data.local_nodes = local_sizes.local_nodes;
  elliptic_data.field_types = &field_type;
  elliptic_data.num_of_fields = 1;

  d4est_mesh_init_field
    (
     p4est,
     sol,
     poly_vec_fcn,
     d4est_ops, // unnecessary?
     d4est_geom, // unnecessary?
     d4est_factors,
     INIT_FIELD_ON_LOBATTO,
     NULL
    );  

  d4est_ghost_data_t* d4est_ghost_data = d4est_ghost_data_init(p4est,
                                                               d4est_ghost,
                                                               &field_type,
                                                               1);

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

    
    d4est_solver_multigrid_t* mg_data = d4est_solver_multigrid_data_init
                                        (
                                         p4est,
                                         d4est_ops,
                                         d4est_geom,
                                         d4est_quad,
                                         &d4est_ghost,
                                         &d4est_ghost_data,
                                         d4est_factors,
                                         initial_grid_input,
                                         input_file
                                        );

    d4est_solver_multigrid_solve
      (
       p4est,
       &elliptic_data,
       &prob_fcns_for_lhs,
       mg_data
      );

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

    if(p4est->mpirank == 0){
      printf("l2 norm = %.15f\n",l2);
    }

    int analyze_matrix;
    printf("Analyze matrix? (0/1): ");
    scanf("%d",&analyze_matrix);
    
   if(analyze_matrix){
  zlog_category_t* matrix_analysis = zlog_get_category("d4est_multigrid_matrix_analysis");
      zlog_info(matrix_analysis,"ANALYZING THE MATRIX\n");
      
      int local_nodes = d4est_factors->local_sizes.local_nodes;
      
      D4EST_ASSERT(p4est->mpisize == 1);
      double* mg_mat = P4EST_ALLOC(double, local_nodes*local_nodes);
      double* mg_mat_trans = P4EST_ALLOC(double, local_nodes*local_nodes);

      d4est_solver_multigrid_get_matrix
        (
         p4est,
         &prob_fcns_for_lhs,
         &elliptic_data,
         mg_data,
         mg_mat
        );
      
      d4est_linalg_mat_transpose(mg_mat, mg_mat_trans, local_nodes);

      double biggest_err;
      int biggest_id;
      double* eig_vals = P4EST_ALLOC(double, local_nodes);
      
      d4est_util_find_biggest_error(mg_mat, mg_mat_trans, local_nodes*local_nodes, &biggest_err, &biggest_id);
      printf("biggest_err = %.15f\n", biggest_err);
      printf("biggest_id = %d\n", biggest_id);
      printf("mg_mat[id] = %.15f\n", mg_mat[biggest_id]);
      printf("mg_mat_trans[id] = %.15f\n", mg_mat_trans[biggest_id]);
      if (biggest_err < 1e-10){
        d4est_linalg_sym_eigvecs
          (mg_mat,eig_vals, local_nodes);

        double min = d4est_util_min_dbl_array
                     (
                      eig_vals,
                      local_nodes);
        zlog_info(matrix_analysis, "Minimum eigenvalue = %f\n", min);
        zlog_info(matrix_analysis, "Biggest symmetry err = %f\n", biggest_err);
        if (min <= 0){        
          DEBUG_PRINT_ARR_DBL(eig_vals, local_nodes);
          printf("Not SPD\n");
        }
      }
      /* else { */
        /* printf("biggest_err = %.15f\n", biggest_err); */
        /* printf("biggest_id = %d\n", biggest_id); */
        /* printf("mg_mat[id] = %.15f\n", mg_mat[biggest_id]); */
        /* printf("mg_mat_trans[id] = %.15f\n", mg_mat_trans[biggest_id]); */
        /* D4EST_ABORT("Mg_mat not symmetric"); */
      /* } */

      
      P4EST_FREE(mg_mat);
      P4EST_FREE(mg_mat_trans);
    }
 

    
    d4est_solver_multigrid_data_destroy(mg_data);

  d4est_laplacian_flux_destroy(flux_data_for_apply_lhs);
  d4est_laplacian_flux_destroy(flux_data_for_residual);

  
  P4EST_FREE(u);
  P4EST_FREE(rhs);
  P4EST_FREE(r);  
  P4EST_FREE(sol);  
  P4EST_FREE(error);  

  d4est_mesh_initial_extents_destroy(initial_grid_input);
  d4est_mesh_data_destroy(d4est_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  
  if (d4est_ghost) {
    d4est_ghost_destroy(d4est_ghost);
  }


  /* if (d4est_ghost_data != NULL){ */
  /*   d4est_ghost_data_destroy(d4est_ghost_data); */
  /*   d4est_ghost_data = NULL; */
  /* } */

  
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);

  P4EST_FREE(input_file);
  zlog_fini();
  PetscFinalize();
  return 0;
}
