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
#include <d4est_vtk.h>
#include <limits.h>
#include <zlog.h>

#define D4EST_REAL_EPS 100*1e-15

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
  
  /* sc_MPI_Comm mpicomm; */
  /* PetscInitialize(&argc,&argv,(char*)0,NULL); */
  /* mpicomm = PETSC_COMM_WORLD; */

  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;
  
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
    printf("[D4EST_INFO]: options file = %s\n", (argc == 2) ? argv[1] :      "d4est_test_initial_mesh.input");
 
 
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    (argc == 2) ? argv[1] :      "d4est_test_initial_mesh.input",
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse((argc == 2) ? argv[1] :      "d4est_test_initial_mesh.input", d4est_geom);

     
  p4est_t* p4est = p4est_new_ext
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

  if (initial_grid_input->refine_per_region_at_start){
    for (int i = 0; i < initial_grid_input->number_of_regions;
         i++){

      int reg = -1;
      for (int k = 0; k < initial_grid_input->number_of_regions; k++){
        if (initial_grid_input->per_region_order[k] == -1){
          D4EST_ABORT("initial_grid_input->per_region_order[k] == -1");
        }
        if (initial_grid_input->per_region_order[k] == i){
          reg = k;
        }
      }
        
      int num_of_times = initial_grid_input->per_region_number_of_refines[reg];
      if (num_of_times == -1){
        D4EST_ABORT("num_of_times == -1");
      }
      int partition = initial_grid_input->per_region_partition[reg];
      if (partition == -1){
        D4EST_ABORT("partition == -1");
      }
                
      for (int j = 0; j < num_of_times; j++){
        void* tmp = p4est->user_pointer;
        d4est_mesh_refine_in_region_data_t rird;
        rird.d4est_geom = d4est_geom;
        rird.region = reg;        
        p4est->user_pointer = &rird;
        p4est_refine(
                     p4est,
                     0,
                     d4est_mesh_refine_in_region_callback,
                     NULL
        );

        p4est->user_pointer = tmp;
      }

      if (partition){
        if (initial_grid_input->keep_quad_fams_together){
          p4est_partition(p4est, 1, NULL);
        }
        else {
          p4est_partition(p4est, 0, NULL);
        }          
      }
        
    }
  }
  else {
    if (initial_grid_input->keep_quad_fams_together){
      p4est_partition(p4est, 1, NULL);
    }
    else {
      p4est_partition(p4est, 0, NULL);
    }
  }
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  
  d4est_ghost_t* d4est_ghost = NULL;
     
  if (proc_rank == 0){
    printf("[D4EST_INFO]: mpisize = %d\n", proc_size);
  }

  sc_MPI_Barrier(mpicomm);
  printf("[D4EST_INFO]: elements on proc %d = %d\n", proc_rank, p4est->local_num_quadrants);
  sc_MPI_Barrier(mpicomm);
  
  /* start just-in-time dg-math */
  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_data_t* d4est_factors = d4est_mesh_data_init(p4est);
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, (argc == 2) ? argv[1] :      "d4est_test_initial_mesh.input", "quadrature");
  

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


  d4est_field_type_t field_type = NODAL;
  d4est_ghost_data_t* d4est_ghost_data = d4est_ghost_data_init(p4est,
                                                               d4est_ghost,
                                                               &field_type,
                                                               1);

    
  int local_nodes = local_sizes.local_nodes;
  dirichlet_bndry_eval_method_t eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;

  /* / Setup boundary conditions */
  d4est_laplacian_dirichlet_bc_t bc_data_for_lhs;
  bc_data_for_lhs.dirichlet_fcn = poly_vec_fcn;
  bc_data_for_lhs.eval_method = eval_method;

    
  d4est_laplacian_flux_data_t* flux_data_for_apply_lhs = d4est_laplacian_flux_new(p4est, "d4est_test_initial_mesh.input", BC_DIRICHLET, &bc_data_for_lhs);


  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = NULL;
  prob_fcns.apply_lhs = d4est_test_apply_lhs;
  prob_fcns.user = flux_data_for_apply_lhs;



  double* poly_vec = P4EST_ALLOC_ZERO(double, local_nodes);
  double* Apoly_vec = P4EST_ALLOC_ZERO(double, local_nodes);

  d4est_elliptic_data_t elliptic_data;
  elliptic_data.u = poly_vec;
  elliptic_data.Au = Apoly_vec;
  elliptic_data.local_nodes = local_nodes;

  elliptic_data.field_types = &field_type;
  elliptic_data.num_of_fields = 1;

  d4est_mesh_init_field
    (
     p4est,
     poly_vec,
     poly_vec_fcn,
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
     &elliptic_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors
    );

  d4est_vtk_save
    (
     p4est,
     d4est_ops,
     "d4est_test_initial_mesh.input",
     "d4est_vtk",
     (const char*[]){"prob_vecs.u", "prob_vecs.Au", NULL},
     (double**)((const double*[]){elliptic_data.u, elliptic_data.Au, NULL}),
     NULL,
     NULL,
     NULL,
     NULL,
     -1
    );

  P4EST_FREE(Apoly_vec);
  P4EST_FREE(poly_vec);

  if (d4est_ghost_data != NULL){
    d4est_ghost_data_destroy(d4est_ghost_data);
    d4est_ghost_data = NULL;
  }

  d4est_laplacian_flux_destroy(flux_data_for_apply_lhs);
    
  d4est_ghost_destroy(d4est_ghost);
  d4est_mesh_initial_extents_destroy(initial_grid_input);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  
  //PetscFinalize();
}
