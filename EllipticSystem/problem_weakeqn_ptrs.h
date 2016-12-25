#ifndef WEAKEQN_PTRS_H
#define WEAKEQN_PTRS_H 

#include "problem_data.h"
#include "../dGMath/dgmath.h"

typedef
void (*problem_apply_op_fcn_t)
(
 p4est_t*,
 p4est_ghost_t*,
 element_data_t*,
 problem_data_t*,
 dgmath_jit_dbase_t*
);


typedef
void (*curved_problem_apply_op_fcn_t)
(
 p4est_t*,
 p4est_ghost_t*,
 curved_element_data_t*,
 problem_data_t*,
 dgmath_jit_dbase_t*,
 d4est_geometry_t*
);

typedef struct {

  /* function pointer to apply LHS of weak equations */
  problem_apply_op_fcn_t apply_lhs;

  /* function pointer to build residual */
  problem_apply_op_fcn_t build_residual;

} weakeqn_ptrs_t;

typedef struct {

  /* function pointer to apply LHS of weak equations */
  curved_problem_apply_op_fcn_t apply_lhs;

  /* function pointer to build residual */
  curved_problem_apply_op_fcn_t build_residual;

} curved_weakeqn_ptrs_t;

#endif
