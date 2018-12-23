#include <d4est_solver_multigrid_bottom_solver_krylov_petsc.h>
#include <d4est_solver_krylov_petsc.h>


static void
d4est_solver_multigrid_bottom_solver_krylov_petsc
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 double* r
)
{

  d4est_solver_multigrid_t* mg_data = p4est->user_pointer;
  d4est_solver_multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  d4est_solver_krylov_petsc_params_t* params = mg_data->bottom_solver->user;

  d4est_solver_krylov_petsc_solve(p4est,
                     vecs,
                     fcns,
                     &updater->current_d4est_ghost,
                     &updater->current_d4est_ghost_data,
                     mg_data->d4est_ops,
                     mg_data->d4est_geom,
                     mg_data->d4est_quad,
                     updater->current_d4est_factors,
                     params,
                     NULL,-1);
}


d4est_solver_multigrid_bottom_solver_t*
d4est_solver_multigrid_bottom_solver_krylov_petsc_init
(
 p4est_t* p4est,
 const char* input_file
)
{
  zlog_category_t *c_default = zlog_get_category("solver_d4est_solver_multigrid_bottom");
  if (p4est->mpirank == 0)
    zlog_info(c_default, "Initializing d4est_solver_multigrid bottom solver...");

  d4est_solver_multigrid_bottom_solver_t* bottom_solver = P4EST_ALLOC(d4est_solver_multigrid_bottom_solver_t, 1);
  d4est_solver_krylov_petsc_params_t* params = P4EST_ALLOC(d4est_solver_krylov_petsc_params_t, 1);

  d4est_solver_krylov_petsc_input
    (
     p4est,
     input_file,
     "mg_bottom_solver_krylov_petsc",
     params
    );
  
  bottom_solver->user = params;
  bottom_solver->solve = d4est_solver_multigrid_bottom_solver_krylov_petsc;
  bottom_solver->update = NULL;

  if (p4est->mpirank == 0)
    zlog_info(c_default, "Initialization of d4est_solver_multigrid bottom solver complete.");

  return bottom_solver;
}


void
d4est_solver_multigrid_bottom_solver_krylov_petsc_destroy(d4est_solver_multigrid_bottom_solver_t* solver){
  P4EST_FREE(solver->user);
  solver->solve = NULL;
  P4EST_FREE(solver);
}
