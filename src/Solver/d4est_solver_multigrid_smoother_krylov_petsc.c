#include <d4est_solver_multigrid_smoother_krylov_petsc.h>
#include <d4est_solver_krylov_petsc.h>
#include <d4est_linalg.h>

static void
d4est_solver_multigrid_smoother_krylov_petsc
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 double* r,
 int level
)
{

  d4est_solver_multigrid_data_t* mg_data = p4est->user_pointer;
  d4est_operators_t* d4est_ops = mg_data->d4est_ops;
  d4est_solver_multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  d4est_solver_krylov_petsc_params_t* params = mg_data->smoother->user;
  d4est_ghost_t* d4est_ghost = updater->current_d4est_ghost;
  d4est_ghost_data_t* d4est_ghost_data = updater->current_d4est_ghost_data;
  
  d4est_solver_krylov_petsc_solve(p4est,
                     vecs,
                     fcns,
                     &d4est_ghost,
                     &d4est_ghost_data,
                     mg_data->d4est_ops,
                     mg_data->d4est_geom,
                     mg_data->d4est_quad,
                     updater->current_d4est_factors,
                     params,
                     NULL,-1);

  double* Au;
  double* rhs;
  double local_nodes = vecs->local_nodes;
  Au = vecs->Au;
  rhs = vecs->rhs;

  d4est_elliptic_eqns_apply_lhs
    (
     p4est,
     d4est_ghost,
     d4est_ghost_data,
     fcns,
     vecs,
     mg_data->d4est_ops,
     mg_data->d4est_geom,
     mg_data->d4est_quad,
     updater->current_d4est_factors
    );
  
  d4est_util_copy_1st_to_2nd(Au, r, local_nodes);
  d4est_linalg_vec_xpby(rhs, -1., r, local_nodes);
}

d4est_solver_multigrid_smoother_t*
d4est_solver_multigrid_smoother_krylov_petsc_init
(
 p4est_t* p4est,
 const char* input_file
)
{
  zlog_category_t *c_default = zlog_get_category("solver_d4est_solver_multigrid_smoother");
  if (p4est->mpirank == 0)
    zlog_info(c_default, "Initializing d4est_solver_multigrid smoother solver...");

  d4est_solver_multigrid_smoother_t* smoother = P4EST_ALLOC(d4est_solver_multigrid_smoother_t, 1);
  d4est_solver_krylov_petsc_params_t* params = P4EST_ALLOC(d4est_solver_krylov_petsc_params_t, 1);

  d4est_solver_krylov_petsc_input
    (
     p4est,
     input_file,
     "mg_smoother_krylov_petsc",
     params
    );
  
  smoother->user = params;
  smoother->smooth = d4est_solver_multigrid_smoother_krylov_petsc;
  smoother->update = NULL;

  if (p4est->mpirank == 0)
    zlog_info(c_default, "Initialization of d4est_solver_multigrid smoother solver complete.");

  return smoother;
}


void
d4est_solver_multigrid_smoother_krylov_petsc_destroy(d4est_solver_multigrid_smoother_t* solver){
  P4EST_FREE(solver->user);
  solver->smooth = NULL;
  P4EST_FREE(solver);
}
