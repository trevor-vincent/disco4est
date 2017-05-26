#include <multigrid_bottom_solver_krylov_petsc.h>
#include <krylov_petsc.h>


static void 
multigrid_bottom_solver_krylov_petsc
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 double* r
)
{

  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  krylov_petsc_params_t* params = mg_data->bottom_solver->user;
  d4est_geometry_t* d4est_geom = updater->d4est_geom;

  krylov_petsc_solve(p4est,
                     vecs,
                     fcns,
                     updater->ghost,
                     updater->ghost_data,
                     mg_data->d4est_ops,
                     d4est_geom,
                     params,
                     NULL);

  /* if r is needed which I don't think it is */
  /* double* Au;  */
  /* double* rhs;  */
  /* local_nodes = vecs->local_nodes; */
  /* Au = vecs->Au; */
  /* rhs = vecs->rhs; */

  /* fcns->apply_lhs(p4est, ghost, ghost_data, vecs, d4est_ops, d4est_geom); */
  /* linalg_copy_1st_to_2nd(Au, r, local_nodes);   */
  /* linalg_vec_xpby(rhs, -1., r, local_nodes); */
}


multigrid_bottom_solver_t*
multigrid_bottom_solver_krylov_petsc_init
(
 p4est_t* p4est,
 const char* input_file
)
{
  multigrid_bottom_solver_t* bottom_solver = P4EST_ALLOC(multigrid_bottom_solver_t, 1);
  krylov_petsc_params_t* params = P4EST_ALLOC(krylov_petsc_params_t, 1);

  krylov_petsc_input
    (
     p4est,
     input_file,
     "mg_bottom_solver_krylov_petsc",
     "[MG_BOTTOM_SOLVER_KRYLOV_PETSC]",
     params
    );
  
  bottom_solver->user = params;
  bottom_solver->solve = multigrid_bottom_solver_krylov_petsc;
  bottom_solver->update = NULL;

  return bottom_solver;
}


void
multigrid_bottom_solver_krylov_petsc_destroy(multigrid_bottom_solver_t* solver){
  P4EST_FREE(solver->user);
  solver->solve = NULL;
  P4EST_FREE(solver);
}
