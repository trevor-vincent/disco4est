#ifndef BI_ESTIMATOR_H
#define BI_ESTIMATOR_H 

#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"
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
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 penalty_calc_t u_penalty_fcn,
 penalty_calc_t u_dirichlet_penalty_fcn,
 penalty_calc_t gradu_penalty_fcn,
 grid_fcn_t u_bndry_fcn,
 double penalty_prefactor,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

#endif
