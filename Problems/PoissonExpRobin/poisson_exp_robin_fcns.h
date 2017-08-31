#ifndef POISSON_EXP_ROBIN_FCNS_H
#define POISSON_EXP_ROBIN_FCNS_H 


typedef struct {

  d4est_poisson_flux_data_t* flux_data_for_apply_lhs;
  d4est_poisson_flux_data_t* flux_data_for_build_rhs;

} problem_ctx_t;


double
poisson_exp_robin_coeff_fcn
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
  double nsum = 0.;
  for (int i = 0; i < (P4EST_DIM); i++){
    nsum += boundary_data->n_on_f_m_quad[i][mortar_node];
  }
  return -nsum;
}


double
poisson_exp_robin_bc_rhs_fcn
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
poisson_exp_robin_analytic_solution
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  double ret = exp(x)*exp(y);
#if (P4EST_DIM)==3
  ret *= exp(z);
#endif
  return ret;
}

static double
poisson_exp_robin_build_rhs_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  double ret = exp(x)*exp(y);
#if (P4EST_DIM)==3
  ret *= exp(z);
#endif
  return -(P4EST_DIM)*ret;
}

static double
poisson_exp_robin_initial_guess
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

static void
poisson_exp_robin_apply_lhs
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
poisson_exp_robin_build_residual
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
  poisson_exp_robin_apply_lhs(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom, d4est_quad, user);
  d4est_linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);
}



#endif
