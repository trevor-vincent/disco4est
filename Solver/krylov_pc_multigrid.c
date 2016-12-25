
#include <krylov_pc_multigrid.h>
#include <multigrid.h>
#include <linalg.h>

krylov_pc_t*
krylov_pc_multigrid_create(krylov_pc_ctx_t* kct){

  krylov_pc_t* pc = P4EST_ALLOC(krylov_pc_t, 1);
  
  pc->pc_apply = krylov_pc_multigrid_apply;
  pc->pc_setup = krylov_pc_multigrid_setup;
  pc->pc_ctx = kct;

  return pc;
}

void
krylov_pc_multigrid_setup
(
 krylov_pc_ctx_t* pc_ctx
)
{
}

void
krylov_pc_multigrid_destroy(krylov_pc_t* pc){
  pc->pc_apply = NULL;
  pc->pc_setup = NULL;
  pc->pc_ctx = NULL;
  P4EST_FREE(pc);
}


void
krylov_pc_multigrid_apply(krylov_pc_ctx_t* kct, double* xp, double* yp)
{
  int local_nodes = kct->vecs->local_nodes;
  multigrid_data_t* mg_data = (multigrid_data_t*)kct->pc_data;
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
     kct->ghost_data
    );
      
  P4EST_FREE(Au);  
}





/* krylov_pc_multigrid.c ends here */
