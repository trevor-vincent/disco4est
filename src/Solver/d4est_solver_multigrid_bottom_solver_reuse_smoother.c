#include <pXest.h>
#include <d4est_solver_multigrid_bottom_solver_reuse_smoother.h>
#include <d4est_solver_multigrid.h>

d4est_solver_multigrid_bottom_solver_t*
d4est_solver_multigrid_bottom_solver_reuse_smoother_init
(
){
  d4est_solver_multigrid_bottom_solver_t* bottom_solver
    = P4EST_ALLOC(d4est_solver_multigrid_bottom_solver_t,
                  1);
  
  bottom_solver->solve = d4est_solver_multigrid_bottom_solver_reuse_smoother_solve;
  bottom_solver->update = NULL;
  bottom_solver->user = NULL;
  return bottom_solver;
  
}
void
d4est_solver_multigrid_bottom_solver_reuse_smoother_solve
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 double* r
)
{
  d4est_solver_multigrid_t* mg_data = p4est->user_pointer;
  d4est_solver_multigrid_smoother_t* smoother = mg_data->smoother;
  smoother->smooth
    (
     p4est,
     vecs,
     fcns,
     r,
     0
    );  
}
void
d4est_solver_multigrid_bottom_solver_reuse_smoother_destroy
(
 d4est_solver_multigrid_bottom_solver_t* bottom_solver
)
{
  P4EST_FREE(bottom_solver);
}
