#ifndef POISSON_SINX_FCNS_H
#define POISSON_SINX_FCNS_H 

#include <pXest.h>

#define PI 3.14159265358932384626433832795

typedef struct {

  d4est_poisson_flux_data_t* flux_data_for_apply_lhs;
  d4est_poisson_flux_data_t* flux_data_for_build_rhs;

} problem_ctx_t;


static double
poisson_sinx_analytic_solution
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
#if (P4EST_DIM)==3
  return sin((PI)*x)*sin((PI)*y)*sin((PI)*z);
#else
  return sin((PI)*x)*sin((PI)*y);
#endif
}

static double
poisson_sinx_rhs_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return (PI)*(PI)*(P4EST_DIM)*poisson_sinx_analytic_solution(x,y,z,user);
}

static double
poisson_sinx_initial_guess
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return 1.;
}

static double
poisson_sinx_boundary_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return poisson_sinx_analytic_solution(x,y,z,user);
}


static void
poisson_sinx_apply_lhs
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* user
)
{
  /* DEBUG_PRINT_ARR_DBL_SUM(prob_vecs->u, prob_vecs->local_nodes);  */
  /* DEBUG_PRINT_ARR_DBL_SUM(prob_vecs->rhs, prob_vecs->local_nodes);  */
  problem_ctx_t* ctx = user;
  d4est_poisson_flux_data_t* flux_fcn_data = ctx->flux_data_for_apply_lhs;
  d4est_poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, flux_fcn_data, d4est_ops, d4est_geom, d4est_quad, d4est_factors);
  /* DEBUG_PRINT_ARR_DBL_SUM(prob_vecs->Au, prob_vecs->local_nodes);  */
}


static void
poisson_sinx_build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* user
)
{
  poisson_sinx_apply_lhs(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom, d4est_quad, d4est_factors, user);
  d4est_linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);
}

#endif
