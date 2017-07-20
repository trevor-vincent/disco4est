#ifndef POISSON_SINX_FCNS_H
#define POISSON_SINX_FCNS_H 

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
  return sin((M_PI)*x)*sin((M_PI)*y)*sin((M_PI)*z);
#else
  return sin((M_PI)*x)*sin((M_PI)*y);
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
  return (M_PI)*(M_PI)*(P4EST_DIM)*poisson_sinx_analytic_solution(x,y,z,user);
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
  return 0.;
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
 void* user
)
{
  d4est_poisson_flux_data_t* flux_fcn_data = user;
  d4est_poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, flux_fcn_data, d4est_ops, d4est_geom, d4est_quad);
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
 void* user
)
{
  poisson_sinx_apply_lhs(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom, d4est_quad, user);
  d4est_linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);
}

#endif
