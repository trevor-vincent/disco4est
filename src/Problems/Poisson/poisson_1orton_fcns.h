#ifndef POISSON_1ORTON_FCNS_H
#define POISSON_1ORTON_FCNS_H 

#include <pXest.h>

#define PI 3.14159265358932384626433832795

typedef struct {

  int n;
  d4est_poisson_flux_data_t* flux_data_for_apply_lhs;
  d4est_poisson_flux_data_t* flux_data_for_build_rhs;

} problem_ctx_t;


double
poisson_1orton_robin_coeff_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user,
 d4est_poisson_flux_boundary_data_t* boundary_data,
 int mortar_node
)
{
  double n = ((problem_ctx_t*)user)->n;
  double r2 = x*x + y*y + z*z;
  return n/sqrt(r2);
}

double
poisson_1orton_robin_bc_rhs_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user,
 d4est_poisson_flux_boundary_data_t* boundary_data,
 int mortar_node
)
{
  return 0.;
}


static double
poisson_1orton_analytic_solution
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  int n = ((problem_ctx_t*)user)->n;
  double r2 = x*x + y*y + z*z;
  return 1./pow(r2, (double)n/2.);
}

static double
poisson_1orton_rhs_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  double n = ((problem_ctx_t*)user)->n;
  double r2 = x*x + y*y + z*z;

  /* return -1 * nabla^2  */
  return -1 * (-1 + n) * n * 1./(pow((r2),(1. + n/2.)));  
}


static double
poisson_1orton_boundary_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return poisson_1orton_analytic_solution
    (
     x,
     y,
#if (P4EST_DIM)==3
     z,
#endif
     user
    );
}


static double
poisson_1orton_initial_guess
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
poisson_1orton_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return poisson_1orton_analytic_solution(x,y,z,user);
}


static void
poisson_1orton_apply_lhs
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
  problem_ctx_t* ctx = user;
  d4est_poisson_flux_data_t* flux_fcn_data = ctx->flux_data_for_apply_lhs;
  d4est_poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, flux_fcn_data, d4est_ops, d4est_geom, d4est_quad, d4est_factors);
}


static void
poisson_1orton_build_residual
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
  poisson_1orton_apply_lhs(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom, d4est_quad, d4est_factors, user);
  d4est_linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);
}

#endif
