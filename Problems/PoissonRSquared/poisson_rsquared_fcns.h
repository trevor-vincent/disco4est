#ifndef POISSON_RSQUARED_FCNS_H
#define POISSON_RSQUARED_FCNS_H 


#include <pXest.h>

#define PI 3.14159265358932384626433832795

typedef struct {

  d4est_poisson_flux_data_t* flux_data_for_apply_lhs;
  d4est_poisson_flux_data_t* flux_data_for_build_rhs;

} problem_ctx_t;


static double
poisson_rsquared_analytic_solution
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  double ret = x*x + y*y;
#if (P4EST_DIM)==3
  ret += z*z;
#endif
}

static double
poisson_rsquared_rhs_fcn
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
  return -6.;
#else
  return -4.;
#endif
}

static double
poisson_rsquared_initial_guess
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
poisson_rsquared_boundary_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return poisson_rsquared_analytic_solution(x,y,z,user);
}


static void
poisson_rsquared_apply_lhs
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* user
)
{
  problem_ctx_t* ctx = user;
  d4est_poisson_flux_data_t* flux_fcn_data = ctx->flux_data_for_apply_lhs;
  d4est_poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, flux_fcn_data, d4est_ops, d4est_geom, d4est_quad);
}


static void
poisson_rsquared_build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* user
)
{
  poisson_rsquared_apply_lhs(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom, d4est_quad, user);
  d4est_linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);
}



#endif
