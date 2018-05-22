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
#include <zlog.h>
#include "../Problems/TwoPunctures/two_punctures_fcns.h"

#define D4EST_REAL_EPS 100*1e-15
#define NUM_OF_TRIALS 5

typedef struct {
  
  int deg;
  int deg_volume_quad;
  int deg_mortar_quad;
  
} test_d4est_puncture_fcns_t;


static int
get_deg_mortar_quad
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  test_d4est_puncture_fcns_t* data = user_ctx;
  return data->deg_mortar_quad;
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
  double r2 = (x*x + y*y + z*z);
  if (r2 == 0)
    return 1.;
  else
    return sin(r2)/sqrt(r2);
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
  return NAN;
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

/*  */

static double
test_d4est_twopunctures_initial_guess
(
 double x,
 double y,
 double z,
 void* user
)
{
  return poly_vec_fcn(x,y,z,user);
}

static void
problem_set_degrees_amr
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  test_d4est_puncture_fcns_t* data = user_ctx;
  elem_data->deg = data->deg;
  if (elem_data->tree != 12)
    elem_data->deg_quad = data->deg_volume_quad;
  else
    elem_data->deg_quad = data->deg;
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

  
  const char* input_file = "test_d4est_puncture_fcns.input";
  /*  */
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    input_file,
                                                    "geometry",
                                                    c_geom);



  /*  */
  p4est_t* p4est = problem_build_p4est
                   (
                    mpicomm,
                    d4est_geom->p4est_conn,
                    -1,
                    0,
                    1
                   );



  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_data_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t, 1);
  d4est_quad->quad_type = QUAD_TYPE_GAUSS_LEGENDRE;
  d4est_quadrature_legendre_new(d4est_quad, d4est_geom,"");

  test_d4est_puncture_fcns_t deg_data;
  deg_data.deg = atoi(argv[1]);
  double eps = atof(argv[2]);
  int num_vecs_to_try = atoi(argv[3]);
  D4EST_ASSERT(argc == 4);
  
  two_punctures_params_t two_punctures_params;
  init_two_punctures_data(&two_punctures_params);

  d4est_poisson_dirichlet_bc_t bc_data_for_jac;
  bc_data_for_jac.dirichlet_fcn = zero_fcn;

  d4est_poisson_dirichlet_bc_t bc_data_for_res;
  bc_data_for_res.dirichlet_fcn = boundary_fcn;

  
  d4est_poisson_flux_data_t* flux_data_for_jac =
    d4est_poisson_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_for_jac, get_deg_mortar_quad, &deg_data);
  d4est_poisson_flux_data_t* flux_data_for_res = d4est_poisson_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_for_res, get_deg_mortar_quad, &deg_data);
  
  problem_ctx_t ctx;
  ctx.flux_data_for_res = flux_data_for_res;
  ctx.flux_data_for_jac = flux_data_for_jac;
  ctx.two_punctures_params = &two_punctures_params;
  
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);

  d4est_amr_t* d4est_amr
    = d4est_amr_init(
                     p4est,
                     input_file,
                     NULL
    );


  double jac_last_node_error [NUM_OF_TRIALS];
  double jac_last_node [NUM_OF_TRIALS];

  double res_last_node_error [NUM_OF_TRIALS];
  double res_last_node [NUM_OF_TRIALS];
  
  for (int level = 0; level < d4est_amr->num_of_amr_steps; ++level){

        printf("\n LEVEL = %d, elements = %d\n", level, p4est->local_num_quadrants);
    
    for (int i = 0; i < NUM_OF_TRIALS; i++){
      deg_data.deg_volume_quad = deg_data.deg + i;
      deg_data.deg_mortar_quad = deg_data.deg + i;
      D4EST_ASSERT(deg_data.deg > 0 && deg_data.deg_volume_quad >= deg_data.deg && deg_data.deg_mortar_quad >= deg_data.deg);
     
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
                         problem_set_degrees_amr,
                         &deg_data
                        );
    

        d4est_elliptic_eqns_t elliptic_eqns;
        elliptic_eqns.apply_lhs = two_punctures_apply_jac;
        elliptic_eqns.build_residual = two_punctures_build_residual;
        elliptic_eqns.user = &ctx;

        double* poly_vec = P4EST_ALLOC(double, local_nodes);
        d4est_mesh_init_field(p4est, poly_vec, poly_vec_fcn, d4est_ops, d4est_geom, d4est_geom);
        double* Apoly_vec = P4EST_ALLOC(double, local_nodes);

        d4est_elliptic_data_t elliptic_data;
        elliptic_data.u = poly_vec;
        elliptic_data.u0 = poly_vec;
        elliptic_data.Au = Apoly_vec;
        elliptic_data.local_nodes = local_nodes;
 

        d4est_util_fill_array(Apoly_vec, 0., local_nodes);
        two_punctures_apply_jac_add_nonlinear_term
          (
           p4est,
           ghost,
           ghost_data,
           &elliptic_data,
           d4est_ops,
           d4est_geom,
           d4est_quad,
           &ctx
          );

        jac_last_node[i] = Apoly_vec[local_nodes-1];
        if(i == 0){
          jac_last_node_error[i] = -1.;
        }
        else{
          jac_last_node_error[i] = fabs(jac_last_node[i] - jac_last_node[i-1]);
        }
              
        d4est_util_fill_array(Apoly_vec, 0., local_nodes);
        two_punctures_build_residual_add_nonlinear_term
          (
           p4est,
           ghost,
           ghost_data,
           &elliptic_data,
           d4est_ops,
           d4est_geom,
           d4est_quad,
           &ctx
          );

        res_last_node[i] = Apoly_vec[local_nodes-1];
        if(i == 0){
          res_last_node_error[i] = -1.;
        }
        else{
          res_last_node_error[i] = fabs(res_last_node[i] - res_last_node[i-1]);
        }

        printf("elem %d deg %d deg_v %d deg_m %d jacl %.15f resl %.15f serr %.15f rerr %.15f \n", p4est->local_num_quadrants, deg_data.deg, deg_data.deg_volume_quad, deg_data.deg_mortar_quad, jac_last_node[i], res_last_node[i], jac_last_node_error[i],  res_last_node_error[i]);

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
           test_d4est_twopunctures_initial_guess,
           NULL,
           eps,
           JAC_TEST_FORWARD_DIFFERENCE,
           num_vecs_to_try
          );

        P4EST_FREE(Apoly_vec);
        P4EST_FREE(poly_vec);

    }

    double* tmp = NULL;
    d4est_amr_step
      (
       p4est,
       &ghost,
       &ghost_data,
       d4est_ops,
       d4est_amr,
       NULL,
       NULL
      );

  }

  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }
    
  d4est_poisson_flux_destroy(flux_data_for_jac);
  d4est_poisson_flux_destroy(flux_data_for_res);
  d4est_mesh_geometry_storage_destroy(geometric_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_amr_destroy(d4est_amr);
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);
  
  PetscFinalize();
  /* sc_finalize (); */
}
