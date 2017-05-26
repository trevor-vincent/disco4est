#ifndef JACOBIAN_TESTER_H
#define JACOBIAN_TESTER_H 

#include <problem_data.h>
#include <problem_weakeqn_ptrs.h>


void
jacobian_tester
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 curved_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 weakeqn_ptrs_t* prob_fcns,
 problem_data_t* prob_vecs
);


#endif
