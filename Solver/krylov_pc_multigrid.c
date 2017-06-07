#include <krylov_pc_multigrid.h>
#include <multigrid.h>
#include <d4est_linalg.h>

void
krylov_pc_multigrid_setup
(
 krylov_pc_t* kpc
)
{
  krylov_pc_multigrid_data_t* kpcmgdata = kpc->pc_data;
  if (kpcmgdata->user_setup_fcn != NULL)
    kpcmgdata->user_setup_fcn(kpc);
}



void
krylov_pc_multigrid_destroy(krylov_pc_t* kpc){
  kpc->pc_apply = NULL;
  kpc->pc_setup = NULL;
  kpc->pc_ctx = NULL;
  P4EST_FREE(kpc->pc_data);
  P4EST_FREE(kpc);
}


void
krylov_pc_multigrid_apply(krylov_pc_t* kpc, double* xp, double* yp)
{
  krylov_pc_multigrid_data_t* kpcmgdata = kpc->pc_data;
  multigrid_data_t* mg_data = kpcmgdata->mg_data;
  petsc_ctx_t* kct = kpc->pc_ctx;
  int local_nodes = kct->vecs->local_nodes;
  
  d4est_linalg_fill_vec(yp, 0., local_nodes);

  double* Au = P4EST_ALLOC(double, local_nodes);
  
  problem_data_t mg_prob_vecs;
  problem_data_copy_ptrs(kct->vecs, &mg_prob_vecs);
  mg_prob_vecs.u = yp;
  mg_prob_vecs.rhs = xp;
  mg_prob_vecs.Au = Au;

  multigrid_solve
    (
     kct->p4est,
     &mg_prob_vecs,
     kct->fcns,
     mg_data
    );
      
  P4EST_FREE(Au);  
}


krylov_pc_t*
krylov_pc_multigrid_create
(
 multigrid_data_t* mg_data,
 void(*user_setup_fcn)(krylov_pc_t* kpc)
){

  krylov_pc_t* pc = P4EST_ALLOC(krylov_pc_t, 1);
  krylov_pc_multigrid_data_t* kpcmgdata = P4EST_ALLOC(krylov_pc_multigrid_data_t,1);
  kpcmgdata->mg_data = mg_data;
  kpcmgdata->user_setup_fcn = user_setup_fcn;
  pc->pc_apply = krylov_pc_multigrid_apply;
  if (user_setup_fcn != NULL){
    pc->pc_setup = krylov_pc_multigrid_setup;
  }
  else{
    pc->pc_setup = NULL;
  }
  pc->pc_data = kpcmgdata;

  return pc;
}
