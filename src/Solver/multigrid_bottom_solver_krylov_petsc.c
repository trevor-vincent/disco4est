#include <multigrid_bottom_solver_krylov_petsc.h>
#include <krylov_petsc.h>


static void
multigrid_bottom_solver_krylov_petsc
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 double* r
)
{

  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  krylov_petsc_params_t* params = mg_data->bottom_solver->user;

  krylov_petsc_solve(p4est,
                     vecs,
                     fcns,
                     updater->ghost,
                     updater->ghost_data,
                     mg_data->d4est_ops,
                     mg_data->d4est_geom,
                     mg_data->d4est_quad,
                     updater->current_geometric_factors,
                     params,
                     NULL);

  /* if r is needed which I don't think it is */
  /* double* Au;  */
  /* double* rhs;  */
  /* local_nodes = vecs->local_nodes; */
  /* Au = vecs->Au; */
  /* rhs = vecs->rhs; */

  /* fcns->apply_lhs(p4est, ghost, ghost_data, vecs, d4est_ops, d4est_geom); */
  /* d4est_util_copy_1st_to_2nd(Au, r, local_nodes);   */
  /* d4est_linalg_vec_xpby(rhs, -1., r, local_nodes); */
}


multigrid_bottom_solver_t*
multigrid_bottom_solver_krylov_petsc_init
(
 p4est_t* p4est,
 const char* input_file
)
{
  zlog_category_t *c_default = zlog_get_category("solver_multigrid_bottom");
  if (p4est->mpirank == 0)
    zlog_info(c_default, "Initializing multigrid bottom solver...");

  multigrid_bottom_solver_t* bottom_solver = P4EST_ALLOC(multigrid_bottom_solver_t, 1);
  krylov_petsc_params_t* params = P4EST_ALLOC(krylov_petsc_params_t, 1);

  krylov_petsc_input
    (
     p4est,
     input_file,
     "mg_bottom_solver_krylov_petsc",
     params
    );
  
  bottom_solver->user = params;
  bottom_solver->solve = multigrid_bottom_solver_krylov_petsc;
  bottom_solver->update = NULL;

  if (p4est->mpirank == 0)
    zlog_info(c_default, "Initialization of multigrid bottom solver complete.");

  return bottom_solver;
}


void
multigrid_bottom_solver_krylov_petsc_destroy(multigrid_bottom_solver_t* solver){
  P4EST_FREE(solver->user);
  solver->solve = NULL;
  P4EST_FREE(solver);
}
