#include <krylov_pc_jacobi.h>

krylov_pc_t*
krylov_pc_jacobi_create(krylov_pc_ctx_t* kct){

  krylov_pc_t* pc = P4EST_ALLOC(krylov_pc_t, 1);
  krylov_pc_jacobi_data_t* pc_data = P4EST_ALLOC(krylov_pc_jacobi_data_t, 1);
  
  pc_data->local_nodes = kct->vecs->local_nodes;
  pc_data->inv_aii = P4EST_ALLOC(double, pc_data->local_nodes);
  kct->pc_data = (void*)pc_data;
  
  pc->pc_apply = krylov_pc_jacobi_apply;
  pc->pc_setup = krylov_pc_jacobi_setup;
  pc->pc_ctx = kct;

  return pc;
}

void
krylov_pc_jacobi_setup
(
 krylov_pc_ctx_t* pc_ctx
)
{
  krylov_pc_jacobi_data_t* jacobi_data = (krylov_pc_jacobi_data_t*)pc_ctx->pc_data;

  int local_nodes = jacobi_data->local_nodes;
  double* u_temp = P4EST_ALLOC_ZERO(double, local_nodes);
  double* Au_temp = P4EST_ALLOC(double, local_nodes);

  double* tmp = pc_ctx->vecs->u;
  double* tmp1 = pc_ctx->vecs->Au;
  pc_ctx->vecs->u = u_temp;
  pc_ctx->vecs->Au = Au_temp;
  
  for (int i = 0; i < local_nodes; i++){
    pc_ctx->vecs->u[i] = 1.;
    if (pc_ctx->d4est_geom == NULL){
      ((weakeqn_ptrs_t*)(pc_ctx->fcns))->apply_lhs(pc_ctx->p4est, *(pc_ctx->ghost), *(element_data_t**)(pc_ctx->ghost_data), pc_ctx->vecs, pc_ctx->dgmath_jit_dbase);
    }
    else{
      ((curved_weakeqn_ptrs_t*)(pc_ctx->fcns))->apply_lhs(pc_ctx->p4est, *(pc_ctx->ghost), *(curved_element_data_t**)(pc_ctx->ghost_data), pc_ctx->vecs, pc_ctx->dgmath_jit_dbase, pc_ctx->d4est_geom);
    }
    
    jacobi_data->inv_aii[i] = 1./(pc_ctx->vecs->Au[i]);
    pc_ctx->vecs->u[i] = 0.;
  }

  pc_ctx->vecs->u = tmp;
  pc_ctx->vecs->Au = tmp1;
  P4EST_FREE(u_temp);
  P4EST_FREE(Au_temp);
}

void
krylov_pc_jacobi_destroy(krylov_pc_t* pc){
  pc->pc_apply = NULL;
  pc->pc_setup = NULL;
  
  krylov_pc_jacobi_data_t* pc_data = (krylov_pc_jacobi_data_t*)pc->pc_ctx->pc_data;
  P4EST_FREE(pc_data->inv_aii);
  P4EST_FREE(pc_data);
  pc->pc_ctx = NULL;
  
  P4EST_FREE(pc);
}


void
krylov_pc_jacobi_apply(krylov_pc_ctx_t* kct, double* xp, double* yp)
{
  krylov_pc_jacobi_data_t* pc_data = (krylov_pc_jacobi_data_t*)kct->pc_data;
  for (int i = 0; i < pc_data->local_nodes; i++){
    yp[i] = pc_data->inv_aii[i]*xp[i];
  }
}





/* krylov_pc_jacobi.c ends here */
