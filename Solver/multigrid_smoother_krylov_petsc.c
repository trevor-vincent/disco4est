#include <multigrid_smoother_krylov_petsc.h>
#include <krylov_petsc.h>
#include <d4est_linalg.h>

static void 
multigrid_smoother_krylov_petsc
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 double* r,
 int level
)
{

  multigrid_data_t* mg_data = p4est->user_pointer;
  d4est_operators_t* d4est_ops = mg_data->d4est_ops;
  multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  krylov_petsc_params_t* params = mg_data->smoother->user;
  p4est_ghost_t** ghost = (updater->ghost);
  void** ghost_data = (updater->ghost_data);
  d4est_geometry_t* d4est_geom = updater->d4est_geom;

  krylov_petsc_solve(p4est,
                     vecs,
                     fcns,
                     ghost,
                     ghost_data,
                     d4est_ops,
                     d4est_geom,
                     params,
                     NULL);

  double* Au;
  double* rhs;
  double local_nodes = vecs->local_nodes;
  Au = vecs->Au;
  rhs = vecs->rhs;

  fcns->apply_lhs(p4est, *ghost, *ghost_data, vecs, d4est_ops, d4est_geom);
  d4est_linalg_copy_1st_to_2nd(Au, r, local_nodes);
  d4est_linalg_vec_xpby(rhs, -1., r, local_nodes);
}

multigrid_smoother_t*
multigrid_smoother_krylov_petsc_init
(
 p4est_t* p4est,
 const char* input_file
)
{
  multigrid_smoother_t* smoother = P4EST_ALLOC(multigrid_smoother_t, 1);
  krylov_petsc_params_t* params = P4EST_ALLOC(krylov_petsc_params_t, 1);

  krylov_petsc_input
    (
     p4est,
     input_file,
     "mg_smoother_krylov_petsc",
     "[MG_SMOOTHER_KRYLOV_PETSC]",
     params
    );
  
  smoother->user = params;
  smoother->smooth = multigrid_smoother_krylov_petsc;
  smoother->update = NULL;

  return smoother;
}


void
multigrid_smoother_krylov_petsc_destroy(multigrid_smoother_t* solver){
  P4EST_FREE(solver->user);
  solver->smooth = NULL;
  P4EST_FREE(solver);
}
