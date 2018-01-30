#ifndef STAMM_FCNS_H
#define STAMM_FCNS_H 

typedef struct {

  double gamma_n;
  double gamma_h;
  double gamma_p;
  double sigma;
  double c2x;
  double c2y;
  double c2z;

} stamm_params_t;


typedef struct {

  d4est_amr_smooth_pred_params_t* smooth_pred_params;
  d4est_poisson_flux_data_t* flux_data_for_apply_lhs;
  d4est_poisson_flux_data_t* flux_data_for_build_rhs;
  stamm_params_t* stamm_params;

} problem_ctx_t;

static
int stamm_params_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  stamm_params_t* pconfig = (stamm_params_t*)user;
  if (d4est_util_match_couple(section,"problem",name,"c2x")) {
    D4EST_ASSERT(pconfig->c2x == -1);
    pconfig->c2x = atof(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"c2y")) {
    D4EST_ASSERT(pconfig->c2y == -1);
    pconfig->c2y = atof(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"c2z")) {
    D4EST_ASSERT(pconfig->c2z == -1);
    pconfig->c2z = atof(value);
  }  
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


static
stamm_params_t
stamm_params_input
(
 const char* input_file
)
{
  stamm_params_t input;
  input.c2x = -1;
  input.c2y = -1;
  input.c2z = -1;
  
  if (ini_parse(input_file, stamm_params_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("problem", input.c2x, -1);
  D4EST_CHECK_INPUT("problem", input.c2y, -1);
  D4EST_CHECK_INPUT("problem", input.c2z, -1);
  
  return input;
}


static
double stamm_analytic_solution
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void *user
)
{
  problem_ctx_t* ctx = user;
  stamm_params_t* params = ctx->stamm_params;
  
  double c2x = params->c2x;
  double c2y = params->c2y;
#if (P4EST_DIM)==3
  double c2z = params->c2z;
#endif
  double xp = x - c2x;
  double yp = y - c2y;
#if (P4EST_DIM)==3
  double zp = z - c2z;
#endif
  double rp = sqrt(xp*xp
                   + yp*yp
#if (P4EST_DIM)==3
                   + zp*zp
#endif   
                  );

#if (P4EST_DIM)==2  
  return x*(1.-x)*y*(1.-y)*d4est_util_dbl_pow_int(rp, 3);
#else
  return x*(1.-x)*y*(1.-y)*z*(1.-z)*d4est_util_dbl_pow_int(rp, 3);
#endif
}

static
double stamm_boundary_fcn
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

static
double stamm_boundary_fcn_sphere
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return stamm_analytic_solution(x,y,
#if(P4EST_DIM)==3
                                 z,
#endif
                                 user);
}


static
double stamm_rhs_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  problem_ctx_t* ctx = user;
  stamm_params_t* params = ctx->stamm_params;
  double c2x = params->c2x;
  double c2y = params->c2y;
#if (P4EST_DIM)==3
  double c2z = params->c2z;
#endif

  #if (P4EST_DIM)==2
  if (x == c2x && y == c2y){
    return 0.;
  }

  else{
    double ret = ((d4est_util_dbl_pow_int(c2x,2) + d4est_util_dbl_pow_int(c2y,2) - 2*c2x*x + d4est_util_dbl_pow_int(x,2) - 2*c2y*y + d4est_util_dbl_pow_int(y,2))*(-2*d4est_util_dbl_pow_int(c2x,2)*x - 6*c2y*x - 2*d4est_util_dbl_pow_int(c2y,2)*x + 4*c2x*d4est_util_dbl_pow_int(x,2) + 2*d4est_util_dbl_pow_int(c2x,2)*d4est_util_dbl_pow_int(x,2) + 6*c2y*d4est_util_dbl_pow_int(x,2) +
       2*d4est_util_dbl_pow_int(c2y,2)*d4est_util_dbl_pow_int(x,2) - 2*d4est_util_dbl_pow_int(x,3) - 4*c2x*d4est_util_dbl_pow_int(x,3) + 2*d4est_util_dbl_pow_int(x,4) - 6*c2x*y - 2*d4est_util_dbl_pow_int(c2x,2)*y - 2*d4est_util_dbl_pow_int(c2y,2)*y + 21*x*y + 16*c2x*x*y + 16*c2y*x*y - 29*d4est_util_dbl_pow_int(x,2)*y -
       16*c2y*d4est_util_dbl_pow_int(x,2)*y + 6*c2x*d4est_util_dbl_pow_int(y,2) + 2*d4est_util_dbl_pow_int(c2x,2)*d4est_util_dbl_pow_int(y,2) + 4*c2y*d4est_util_dbl_pow_int(y,2) + 2*d4est_util_dbl_pow_int(c2y,2)*d4est_util_dbl_pow_int(y,2) - 29*x*d4est_util_dbl_pow_int(y,2) - 16*c2x*x*d4est_util_dbl_pow_int(y,2) + 37*d4est_util_dbl_pow_int(x,2)*d4est_util_dbl_pow_int(y,2) -
                                                                                                                          2*d4est_util_dbl_pow_int(y,3) - 4*c2y*d4est_util_dbl_pow_int(y,3) + 2*d4est_util_dbl_pow_int(y,4)))/sqrt(d4est_util_dbl_pow_int(-c2x + x,2) + d4est_util_dbl_pow_int(-c2y + y,2));


  ret *= -1;
  return ret;
  }
#else
  if (x == c2x && y == c2y && z == c2z){
    return 0.;
  }
  else{
    double ret = (-2*(1 - x)*x*(1 - y)*y*d4est_util_dbl_pow_int(d4est_util_dbl_pow_int(c2x - x,2) + d4est_util_dbl_pow_int(c2y - y,2) + d4est_util_dbl_pow_int(c2z - z,2),2) + 
     9*(1 - x)*x*(1 - y)*y*(d4est_util_dbl_pow_int(c2x - x,2) + d4est_util_dbl_pow_int(c2y - y,2) + d4est_util_dbl_pow_int(c2z - z,2))*(1 - z)*
     z + 6*(1 - x)*(-c2x + x)*(1 - y)*y*
     (d4est_util_dbl_pow_int(c2x - x,2) + d4est_util_dbl_pow_int(c2y - y,2) + d4est_util_dbl_pow_int(c2z - z,2))*(1 - z)*z - 
     6*x*(-c2x + x)*(1 - y)*y*(d4est_util_dbl_pow_int(c2x - x,2) + d4est_util_dbl_pow_int(c2y - y,2) + d4est_util_dbl_pow_int(c2z - z,2))*
     (1 - z)*z + 6*(1 - x)*x*(1 - y)*(-c2y + y)*
     (d4est_util_dbl_pow_int(c2x - x,2) + d4est_util_dbl_pow_int(c2y - y,2) + d4est_util_dbl_pow_int(c2z - z,2))*(1 - z)*z - 
     6*(1 - x)*x*y*(-c2y + y)*(d4est_util_dbl_pow_int(c2x - x,2) + d4est_util_dbl_pow_int(c2y - y,2) + d4est_util_dbl_pow_int(c2z - z,2))*
     (1 - z)*z - 2*(1 - x)*x*d4est_util_dbl_pow_int(d4est_util_dbl_pow_int(c2x - x,2) + d4est_util_dbl_pow_int(c2y - y,2) + d4est_util_dbl_pow_int(c2z - z,2),
                                              2)*(1 - z)*z - 2*(1 - y)*y*d4est_util_dbl_pow_int(d4est_util_dbl_pow_int(c2x - x,2) + d4est_util_dbl_pow_int(c2y - y,2) + 
                                                                                          d4est_util_dbl_pow_int(c2z - z,2),2)*(1 - z)*z - 
     3*d4est_util_dbl_pow_int(c2x - x,2)*(-1 + x)*x*(-1 + y)*y*(-1 + z)*z - 
     3*(-1 + x)*x*d4est_util_dbl_pow_int(c2y - y,2)*(-1 + y)*y*(-1 + z)*z - 
     3*(-1 + x)*x*(-1 + y)*y*d4est_util_dbl_pow_int(c2z - z,2)*(-1 + z)*z + 
     6*(1 - x)*x*(1 - y)*y*(d4est_util_dbl_pow_int(c2x - x,2) + d4est_util_dbl_pow_int(c2y - y,2) + d4est_util_dbl_pow_int(c2z - z,2))*(1 - z)*
     (-c2z + z) - 6*(1 - x)*x*(1 - y)*y*
     (d4est_util_dbl_pow_int(c2x - x,2) + d4est_util_dbl_pow_int(c2y - y,2) + d4est_util_dbl_pow_int(c2z - z,2))*z*(-c2z + z))/
      sqrt(d4est_util_dbl_pow_int(c2x - x,2) + d4est_util_dbl_pow_int(c2y - y,2) + d4est_util_dbl_pow_int(c2z - z,2));
    ret *= -1;
    return ret;
  }
#endif
}

static void
stamm_apply_lhs
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
  stamm_params_t* params = ctx->stamm_params;
  d4est_poisson_flux_data_t* flux_fcn_data = ctx->flux_data_for_apply_lhs;
  d4est_poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, flux_fcn_data, d4est_ops, d4est_geom, d4est_quad);

  /* DEBUG_PRINT_ARR_DBL_SUM(prob_vecs->Au, prob_vecs->local_nodes); */
  
}


static void
stamm_build_residual
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
  stamm_apply_lhs(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom, d4est_quad, user);
  d4est_linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);
}

static double
stamm_initial_guess
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return 100.;
}


#endif
