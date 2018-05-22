#include <krylov_pc_multigrid.h>
#include <multigrid.h>
#include <d4est_linalg.h>

void
krylov_pc_multigrid_setup
(
 krylov_pc_t* kpc
)
{
  /* krylov_pc_multigrid_data_t* kpcmgdata = kpc->pc_data; */
  /* if (kpcmgdata->user_setup_fcn != NULL) */
    /* kpcmgdata->user_setup_fcn(kpc); */
}



void
krylov_pc_multigrid_destroy(krylov_pc_t* kpc){
  kpc->pc_apply = NULL;
  kpc->pc_setup = NULL;
  kpc->pc_ctx = NULL;
  kpc->pc_data = NULL;
  /* P4EST_FREE(kpc->pc_data); */
  P4EST_FREE(kpc);
}


void
krylov_pc_multigrid_apply(krylov_pc_t* kpc, double* xp, double* yp)
{
  //  krylov_pc_multigrid_data_t* kpcmgdata = kpc->pc_data;
  multigrid_data_t* mg_data = kpc->pc_data;//kpcmgdata->mg_data;
  krylov_ctx_t* kct = kpc->pc_ctx;
  int local_nodes = kct->vecs->local_nodes;
  /* if(kct->residual_history_size !=0){ */

  int na;
  double* residual_history;
  KSPGetResidualHistory(*kct->ksp,&residual_history,&na);
 
  if(na >=2)
    {
      double ratio = residual_history[na-1]/residual_history[na-2];      
      if (ratio > .95){
        printf("ratio is getting bad, it is = %f\n",ratio);
        d4est_util_copy_1st_to_2nd(xp,yp,local_nodes);
        return;
      }
    }
  
    d4est_util_fill_array(yp, 0., local_nodes);
    double* Au = P4EST_ALLOC(double, local_nodes); 
    d4est_elliptic_data_t mg_prob_vecs;
    d4est_elliptic_data_copy_ptrs(kct->vecs, &mg_prob_vecs);
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
  /* krylov_pc_multigrid_data_t* kpcmgdata = P4EST_ALLOC(krylov_pc_multigrid_data_t,1); */
  /* kpcmgdata->mg_data = mg_data; */
  /* kpcmgdata->user_setup_fcn = user_setup_fcn; */
  
  pc->pc_apply = krylov_pc_multigrid_apply;

  /* if (user_setup_fcn != NULL){ */
    /* pc->pc_setup = krylov_pc_multigrid_setup; */
  /* } */
  /* else{ */
    /* pc->pc_setup = NULL; */
  /* } */
  pc->pc_setup = user_setup_fcn;
  pc->pc_data = mg_data;

  return pc;
}
