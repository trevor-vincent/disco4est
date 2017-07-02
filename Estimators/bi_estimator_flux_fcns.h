#ifndef BI_ESTIMATOR_FLUX_FCNS_H
#define BI_ESTIMATOR_FLUX_FCNS_H 

#include "../Flux/compute_flux.h"
#include "../Utilities/util.h"

inline static double
bi_gradu_prefactor_maxp_minh
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double max_p = util_max_int(deg_m, deg_p);
  double min_h = (h_m < h_p) ? h_m : h_p;
  return sqrt(min_h/max_p);
}

inline static double
bi_u_prefactor_conforming_maxp_minh
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double max_p = (double)util_max_int(deg_m, deg_p);
  double min_h = (h_m < h_p) ? h_m : h_p;
  return sqrt(penalty_prefactor*max_p*max_p/min_h);
}

inline static double
bi_gradu_prefactor_max_h_over_p
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double minus_h_over_p = h_m/(double)deg_m;
  double plus_h_over_p = h_p/(double)deg_p;

  double max_h_over_p = (minus_h_over_p > plus_h_over_p) ?
                        minus_h_over_p :
                        plus_h_over_p;
  
  return sqrt(max_h_over_p);
}


inline static double
bi_u_prefactor_conforming_max_p2_over_h
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double minus_p2_over_h = ((double)deg_m)*((double)deg_m)/h_m;
  double plus_p2_over_h = ((double)deg_p)*((double)deg_p)/h_p;
  double max_p2_over_h = (minus_p2_over_h > plus_p2_over_h) ?
                         minus_p2_over_h :
                         plus_p2_over_h;
  return sqrt(penalty_prefactor*max_p2_over_h);
}



/* inline static double */
/* houston_gradu_prefactor */
/* ( */
/*  int deg_m, */
/*  int deg_p, */
/*  double h, */
/*  double penalty_prefactor */
/* ) */
/* { */
/*   double max_p = util_max_int(deg_m, deg_p); */
/*   return sqrt(.5*h/max_p); */
/* } */

/* inline static double */
/* houston_u_prefactor_conforming */
/* ( */
/*  int deg_m, */
/*  int deg_p, */
/*  double h, */
/*  double penalty_prefactor */
/* ) */
/* { */
/*   double max_p = util_max_int(deg_m, deg_p); */
/*   return sqrt(.5*penalty_prefactor*max_p*max_p/h); */
/* } */

/* inline static double */
/* houston_u_dirichlet_prefactor_conforming */
/* ( */
/*  int deg_m, */
/*  int deg_p, */
/*  double h, */
/*  double penalty_prefactor */
/* ) */
/* { */
/*   double max_p = util_max_int(deg_m, deg_p); */
/*   return sqrt(penalty_prefactor*max_p*max_p/h); */
/* } */




inline static double
houston_gradu_prefactor_max_h_over_p
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double minus_h_over_p = h_m/(double)deg_m;
  double plus_h_over_p = h_p/(double)deg_p;

  double max_h_over_p = (minus_h_over_p > plus_h_over_p) ?
                        minus_h_over_p :
                        plus_h_over_p;
  
  return sqrt(.5*max_h_over_p);
}


inline static double
houston_u_prefactor_max_p2_over_h
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double minus_p2_over_h = ((double)deg_m)*((double)deg_m)/h_m;
  double plus_p2_over_h = ((double)deg_p)*((double)deg_p)/h_p;
  double max_p2_over_h = (minus_p2_over_h > plus_p2_over_h) ?
                         minus_p2_over_h :
                         plus_p2_over_h;
  return sqrt(.5*penalty_prefactor*max_p2_over_h);
}

inline static double
houston_u_dirichlet_prefactor_max_p2_over_h
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double minus_p2_over_h = ((double)deg_m)*((double)deg_m)/h_m;
  double plus_p2_over_h = ((double)deg_p)*((double)deg_p)/h_p;
  double max_p2_over_h = (minus_p2_over_h > plus_p2_over_h) ?
                         minus_p2_over_h :
                         plus_p2_over_h;
  return sqrt(penalty_prefactor*max_p2_over_h);
}



inline static double
houston_gradu_prefactor_maxp_minh
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double max_p = util_max_int(deg_m, deg_p);
  double min_h = (h_m < h_p) ? h_m : h_p;
  return sqrt(.5*min_h/max_p);
}

inline static double
houston_u_prefactor_maxp_minh
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double max_p = (double)util_max_int(deg_m, deg_p);
  double min_h = (h_m < h_p) ? h_m : h_p;
  return sqrt(.5*penalty_prefactor*max_p*max_p/min_h);
}

inline static double
houston_u_dirichlet_prefactor_maxp_minh
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double max_p = (double)util_max_int(deg_m, deg_p);
  double min_h = (h_m < h_p) ? h_m : h_p;
  return sqrt(penalty_prefactor*max_p*max_p/min_h);
}




d4est_mortar_fcn_ptrs_t
bi_est_dirichlet_fetch_fcns
(
 d4est_grid_fcn_t bndry_fcn,
 penalty_calc_t u_penalty_fcn,
 penalty_calc_t u_dirichlet_penalty_fcn,
 penalty_calc_t gradu_penalty_fcn,
 double penalty_prefactor
);


#endif
