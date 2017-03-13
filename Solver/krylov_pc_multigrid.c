#include <krylov_pc_multigrid.h>
#include <multigrid.h>
#include <linalg.h>

krylov_pc_t*
krylov_pc_multigrid_create
(
 multigrid_data_t* mg_data
){

  krylov_pc_t* pc = P4EST_ALLOC(krylov_pc_t, 1);
  
  pc->pc_apply = krylov_pc_multigrid_apply;
  pc->pc_setup = NULL;
  pc->pc_data = mg_data;

  return pc;
}

void
krylov_pc_multigrid_destroy(krylov_pc_t* pc){
  pc->pc_apply = NULL;
  pc->pc_setup = NULL;
  pc->pc_ctx = NULL;
  pc->pc_data = NULL;
  P4EST_FREE(pc);
}


void
krylov_pc_multigrid_apply(void* kpc_in, double* xp, double* yp)
{
  krylov_pc_t* kpc = kpc_in;
  multigrid_data_t* mg_data = (multigrid_data_t*)kpc->pc_data;
  petsc_ctx_t* kct = kpc->pc_ctx;
  int local_nodes = kct->vecs->local_nodes;
  
  linalg_fill_vec(yp, 0., local_nodes);

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
     mg_data,
     kct->ghost,
     (element_data_t**)kct->ghost_data
    );
      
  P4EST_FREE(Au);  
}
