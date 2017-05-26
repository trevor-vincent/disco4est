#ifndef WEAKEQN_PTRS_H
#define WEAKEQN_PTRS_H 

#include "problem_data.h"
#include "../dGMath/d4est_operators.h"

typedef
void (*problem_apply_op_fcn_t)
(
 p4est_t*,
 p4est_ghost_t*,
 void*, /* element_data */
 problem_data_t*,
 d4est_operators_t*,
 d4est_geometry_t*
);

typedef struct {

  /* function pointer to apply LHS of weak equations */
  problem_apply_op_fcn_t apply_lhs;

  /* function pointer to build residual */
  problem_apply_op_fcn_t build_residual;

} weakeqn_ptrs_t;

#endif
