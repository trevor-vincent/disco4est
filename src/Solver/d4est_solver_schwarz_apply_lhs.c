#include <d4est_solver_schwarz_apply_lhs.h>

d4est_solver_schwarz_apply_lhs_t*
d4est_solver_schwarz_apply_lhs_init
(
 d4est_solver_schwarz_apply_lhs_fcn_t fcn,
 void* fcn_ctx
)
{
  d4est_solver_schwarz_apply_lhs_t* apply_lhs
    = P4EST_ALLOC(d4est_solver_schwarz_apply_lhs_t, 1);

  apply_lhs->apply_lhs_fcn = fcn;
  apply_lhs->apply_lhs_ctx = fcn_ctx;
  return apply_lhs;
}


void
d4est_solver_schwarz_apply_lhs_destroy
(
 d4est_solver_schwarz_apply_lhs_t* apply_lhs 
)
{
  P4EST_FREE(apply_lhs);
}
