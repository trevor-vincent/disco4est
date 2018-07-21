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
#include <d4est_elliptic_eqns.h>
#include <d4est_solver_jacobian_tester.h>
#include <d4est_util.h>
#include <limits.h>
#include <ini.h>
#include "../Problems/Okendon/okendon_fcns.h"

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
  /* int deg = TEST_DEG_INIT; */
  /* double poly = pow(x,deg) + pow(y,deg); */
/* #if (P4EST_DIM)==3 */
  /* poly += ((P4EST_DIM)==3) ? pow(z,deg) : 0.; */
/* #endif */
  /* return poly; */

/*   double poly = x*(1.-x)*y*(1.-y); */
/* #if (P4EST_DIM)==3 */
/*   poly *= z*(1.-z); */
/* #endif */
/*   return poly; */

  return x*x - y*y;
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
laplacian_poly_vec_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
){
/*   int deg = TEST_DEG_INIT; */
/*   double poly = pow(x,deg-2) + pow(y,deg-2); */
/* #if (P4EST_DIM)==3 */
/*   poly += ((P4EST_DIM)==3) ? pow(z,deg-2) : 0.; */
/* #endif */
/*   double factor = deg; */
/*   while (factor != 0){ */
/*     poly *= factor; */
/*     factor -= 1; */
/*   } */
/*   return poly; */
  
/* #if (P4EST_DIM)==2 */
/*   return 2*(-x + pow(x,2) + (-1 + y)*y); */
/* #elif (P4EST_DIM)==3 */
/*   return -2*(-1 + x)*x*(-1 + y)*y - 2*(-1 + x)*x*(-1 + z)*z - 2*(-1 + y)*y*(-1 + z)*z; */
/* #else */
/*   D4EST_ABORT(""); */
/* #endif */

  /* return -2.; */
  return 0.;
}


/*  */

static double
testd4est_okendon_initial_guess
(
 double x,
 double y,
 double z,
 void* user
)
{
  return 1.0;
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
  
  d4est_geometry_t* d4est_geom = P4EST_ALLOC(d4est_geometry_t, 1);
  d4est_geom->geom_type = GEOM_BRICK;
  d4est_geom->DX_compute_method = GEOM_COMPUTE_ANALYTIC;
  d4est_geom->JAC_compute_method = GEOM_COMPUTE_NUMERICAL;

  const char* input_file = "testd4est_okendon_jacobian.input";
  d4est_geometry_brick_new(proc_rank, input_file, "geometry", "[Geometry]:", d4est_geom);
    
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

  okendon_params_t okendon_params = okendon_params_init(input_file);
  d4est_poisson_flux_data_t* flux_data_for_jac = d4est_poisson_flux_new(p4est, input_file, zero_fcn, NULL);
  d4est_poisson_flux_data_t* flux_data_for_residual = d4est_poisson_flux_new(p4est, input_file, poly_vec_fcn, NULL);
  okendon_params.flux_data_for_jac = flux_data_for_jac;
  okendon_params.flux_data_for_residual = flux_data_for_residual;
      
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);

  d4est_amr_t* d4est_amr = d4est_amr_init(
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
  
  double* poly_vec = P4EST_ALLOC_ZERO(double, local_nodes);
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

    
    printf("level = %d, elements = %d, nodes = %d\n", level, p4est->local_num_quadrants, local_nodes);


    d4est_elliptic_eqns_t elliptic_eqns;
    elliptic_eqns.apply_lhs = okendon_apply_jac;
    elliptic_eqns.build_residual = okendon_build_residual_strongbc;
    elliptic_eqns.user = &okendon_params;

    int num_vecs_to_try = 5;
    d4est_solver_jacobian_tester
      (
       p4est,
       ghost,
       ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &elliptic_eqns,
       local_nodes,
       testd4est_okendon_initial_guess,
       NULL,
       .01,
       JAC_TEST_FORWARD_DIFFERENCE,
       num_vecs_to_try
      );

    
    d4est_amr_step
      (
       p4est,
       &ghost,
       &ghost_data,
       d4est_ops,
       d4est_amr,
       &poly_vec,
       NULL
      );


  }


  P4EST_FREE(poly_vec);
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
