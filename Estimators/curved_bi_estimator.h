#ifndef CURVED_BI_ESTIMATOR_H
#define CURVED_BI_ESTIMATOR_H 

#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"

void
curved_bi_estimator_compute
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
 curved_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* geom
);

#endif
