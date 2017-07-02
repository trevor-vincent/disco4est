#ifndef BI_ESTIMATOR_H
#define BI_ESTIMATOR_H 

#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/d4est_elliptic_eqns.h"
#include "../Estimators/estimator_stats.h"

double
bi_estimator_test_compute
(
 p4est_t* p4est,
 int level
);

void
bi_estimator_compute
(
 p4est_t* p4est,
 d4est_elliptic_problem_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 penalty_calc_t u_penalty_fcn,
 penalty_calc_t u_dirichlet_penalty_fcn,
 penalty_calc_t gradu_penalty_fcn,
 d4est_grid_fcn_t u_bndry_fcn,
 double penalty_prefactor,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_operators_t* d4est_ops
);

#endif
