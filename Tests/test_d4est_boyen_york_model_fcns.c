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
#include <d4est_poisson.h>
#include <d4est_poisson_flux.h>
#include <d4est_poisson_flux_sipg.h>
#include <d4est_util.h>
#include <ini.h>
#include <limits.h>
#include <zlog.h>
#include "../Problems/BoyenYorkModel/boyen_york_model_fcns.h"

#define D4EST_REAL_EPS 100*1e-15
#if (P4EST_DIM)==2
#define TEST_DEG_INIT 2
#else
#define TEST_DEG_INIT 6
#endif

static void
problem_set_degrees_init
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  /* d4est_element_data_t* elem_data = elem_data_tmp; */
  elem_data->deg = TEST_DEG_INIT;
  elem_data->deg_quad = TEST_DEG_INIT;
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
  sc_MPI_Comm mpicomm;
  PetscInitialize(&argc,&argv,(char*)0,NULL);
  mpicomm = PETSC_COMM_WORLD;
  
  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  p4est_init(NULL, SC_LP_ERROR);
  
 
  const char* input_file = "test_d4est_boyen_york_model_fcns.input";


  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    input_file,
                                                    "geometry",
                                                    c_geom);
  
  p4est_t* p4est = problem_build_p4est
                   (
                    mpicomm,
                    d4est_geom->p4est_conn,
                    -1,
                    1,
                    1
                   );
  
  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_data_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t, 1);
  d4est_quad->quad_type = QUAD_TYPE_GAUSS_LEGENDRE;
  d4est_quadrature_legendre_new(d4est_quad, d4est_geom,"");
  /*  */

  
  boyen_york_model_params_t boyen_york_model_input = boyen_york_model_params(input_file);
  
  d4est_poisson_flux_data_t* flux_data_for_jac = d4est_poisson_flux_new(p4est, input_file, zero_fcn, NULL);
  d4est_poisson_flux_data_t* flux_data_for_residual = d4est_poisson_flux_new(p4est, input_file, boyen_york_model_boundary_fcn, NULL);

  boyen_york_model_input.flux_data_for_jac = flux_data_for_jac;
  boyen_york_model_input.flux_data_for_residual = flux_data_for_residual;

  d4est_poisson_flux_sipg_params_t* sipg_params_for_jac = flux_data_for_jac->user;
  d4est_poisson_flux_sipg_params_t* sipg_params_for_residual = flux_data_for_residual->user;
  sipg_params_for_jac->user = &boyen_york_model_input;
  sipg_params_for_residual->user = &boyen_york_model_input;
  
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);

  d4est_amr_t* d4est_amr = d4est_amr_init
                           (
                            p4est,
                            input_file,
                            NULL
                           );


  int local_nodes = d4est_mesh_update
                    (
                     p4est,
                     ghost,
                     ghost_data,
                     d4est_ops,
                     d4est_geom,
                     d4est_quad,
                     geometric_factors,
                     INITIALIZE_QUADRATURE_DATA,
                     INITIALIZE_GEOMETRY_DATA,
                     INITIALIZE_GEOMETRY_ALIASES,
                     problem_set_degrees_init,
                     NULL
                    );
  
  double* by_vec = P4EST_ALLOC(double, local_nodes);
  int same = 1;
  int same2 = 1;
  
  for (int level = 0; level < d4est_amr->num_of_amr_steps; ++level){

    local_nodes = d4est_mesh_update
                  (
                   p4est,
                   ghost,
                   ghost_data,
                   d4est_ops,
                   d4est_geom,
                   d4est_quad,
                   geometric_factors,
                   INITIALIZE_QUADRATURE_DATA,
                   INITIALIZE_GEOMETRY_DATA,
                   INITIALIZE_GEOMETRY_ALIASES,
                   problem_set_degrees_amr,
                   NULL
                  );

    
      d4est_mesh_init_field(p4est, by_vec, boyen_york_model_initial_guess, d4est_ops, d4est_geom, NULL);
      double* Aby_vec = P4EST_ALLOC(double, local_nodes);

      d4est_elliptic_data_t elliptic_data;
      elliptic_data.u = by_vec;
      elliptic_data.Au = Aby_vec;
      elliptic_data.local_nodes = local_nodes;

      /* d4est_mesh_init_field(p4est, by_vec, boyen_york_model_analytic_solution, d4est_ops, d4est_geom, &boyen_york_model_input); */

      /* DEBUG_PRINT_ARR_DBL(by_vec, local_nodes); */
      
      boyen_york_model_build_residual
        (
         p4est,
         ghost,
         ghost_data,
         &elliptic_data,
         d4est_ops,
         d4est_geom,
         d4est_quad,
         &boyen_york_model_input
        );
      
      /* DEBUG_PRINT_2ARR_DBL(by_vec, Aby_vec, local_nodes); */
      
      P4EST_FREE(Aby_vec);

    
    d4est_amr_step
      (
       p4est,
       &ghost,
       &ghost_data,
       d4est_ops,
       d4est_amr,
       &by_vec,
       NULL
      );
    
  }


  P4EST_FREE(by_vec);

  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }
    
  d4est_poisson_flux_destroy(flux_data_for_jac);
  d4est_poisson_flux_destroy(flux_data_for_residual);
  d4est_mesh_geometry_storage_destroy(geometric_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_amr_destroy(d4est_amr);
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
