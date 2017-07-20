#ifndef STAMM_FCNS_H
#define STAMM_FCNS_H 

static
double stamm_analytic_solution
(
 double x,
 double y
#if (P4EST_DIM)==3
 ,
 double z
#endif
)
{
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
 double y
#if (P4EST_DIM)==3
 ,
 double z
#endif
)
{
  return 0.;
}

static
double stamm_rhs_fcn
(
 double x,
 double y
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  stamm_params_t* params = user;
  double c2x = params->c2x;
  double c2y = params->c2y;
#if (P4EST_DIM)==3
  double c2z = params->c2z;
#endif
  
  #if (P4EST_DIM)==2
  if (x == c2x && y == c2y){
    return 0.;
  }

  else
    return ((d4est_util_dbl_pow_int(c2x,2) + d4est_util_dbl_pow_int(c2y,2) - 2*c2x*x + d4est_util_dbl_pow_int(x,2) - 2*c2y*y + d4est_util_dbl_pow_int(y,2))*(-2*d4est_util_dbl_pow_int(c2x,2)*x - 6*c2y*x - 2*d4est_util_dbl_pow_int(c2y,2)*x + 4*c2x*d4est_util_dbl_pow_int(x,2) + 2*d4est_util_dbl_pow_int(c2x,2)*d4est_util_dbl_pow_int(x,2) + 6*c2y*d4est_util_dbl_pow_int(x,2) +
       2*d4est_util_dbl_pow_int(c2y,2)*d4est_util_dbl_pow_int(x,2) - 2*d4est_util_dbl_pow_int(x,3) - 4*c2x*d4est_util_dbl_pow_int(x,3) + 2*d4est_util_dbl_pow_int(x,4) - 6*c2x*y - 2*d4est_util_dbl_pow_int(c2x,2)*y - 2*d4est_util_dbl_pow_int(c2y,2)*y + 21*x*y + 16*c2x*x*y + 16*c2y*x*y - 29*d4est_util_dbl_pow_int(x,2)*y -
       16*c2y*d4est_util_dbl_pow_int(x,2)*y + 6*c2x*d4est_util_dbl_pow_int(y,2) + 2*d4est_util_dbl_pow_int(c2x,2)*d4est_util_dbl_pow_int(y,2) + 4*c2y*d4est_util_dbl_pow_int(y,2) + 2*d4est_util_dbl_pow_int(c2y,2)*d4est_util_dbl_pow_int(y,2) - 29*x*d4est_util_dbl_pow_int(y,2) - 16*c2x*x*d4est_util_dbl_pow_int(y,2) + 37*d4est_util_dbl_pow_int(x,2)*d4est_util_dbl_pow_int(y,2) -
                                                                                                                          2*d4est_util_dbl_pow_int(y,3) - 4*c2y*d4est_util_dbl_pow_int(y,3) + 2*d4est_util_dbl_pow_int(y,4)))/sqrt(d4est_util_dbl_pow_int(-c2x + x,2) + d4est_util_dbl_pow_int(-c2y + y,2));

#else
  if (x == c2x && y == c2y && z == c2z){
    return 0.;
  }
  else{
    return (-2*(1 - x)*x*(1 - y)*y*d4est_util_dbl_pow_int(d4est_util_dbl_pow_int(c2x - x,2) + d4est_util_dbl_pow_int(c2y - y,2) + d4est_util_dbl_pow_int(c2z - z,2),2) + 
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
  }
#endif
}

static
void
stamm_build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
  poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops, NULL);
  d4est_linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);
}

static
void
stamm_apply_lhs
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
  poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops, NULL);
}

#endif
