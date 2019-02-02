#include <d4est_krylov_pc_multigrid.h>
#include <d4est_solver_multigrid.h>
#include <d4est_linalg.h>
#include <zlog.h>

void
d4est_krylov_pc_multigrid_setup
(
 d4est_krylov_pc_t* kpc
)
{
  d4est_krylov_pc_multigrid_data_t* kpcmgdata = kpc->pc_data;
  d4est_solver_multigrid_t* mg_data = kpcmgdata->mg_data;
  mg_data->linear_operator_updates++;
  krylov_ctx_t* kct = kpc->pc_ctx;
  
  if(kct->p4est->mpirank == 0){
    zlog_category_t *c_def = zlog_get_category("d4est_krylov_pc_multigrid");
    zlog_info(c_def, "Operator is changing, running d4est_krylov_pc_multigrid_setup");
    zlog_info(c_def, "mg_data->linear_operator_updates before setup = %d\n",
              mg_data->linear_operator_updates);
  }
  mg_data->linear_operator_updates++;
  /* add to mg_data setup int every time this is called */
  /* d4est_krylov_pc_multigrid_data_t* kpcmgdata = kpc->pc_data; */
  
  if (kpcmgdata->user_setup_fcn != NULL)
    kpcmgdata->user_setup_fcn(kpc);
}



void
d4est_krylov_pc_multigrid_destroy(d4est_krylov_pc_t* kpc){
  kpc->pc_apply = NULL;
  kpc->pc_setup = NULL;
  kpc->pc_ctx = NULL;
  P4EST_FREE(kpc->pc_data);
  kpc->pc_data = NULL;
  P4EST_FREE(kpc);
}


void
d4est_krylov_pc_multigrid_apply(d4est_krylov_pc_t* kpc, double* xp, double* yp)
{
  d4est_krylov_pc_multigrid_data_t* kpcmgdata = kpc->pc_data;
 d4est_solver_multigrid_t* mg_data = kpcmgdata->mg_data;
  krylov_ctx_t* kct = kpc->pc_ctx;
  int local_nodes = kct->vecs->local_nodes;
  int na;
  double* residual_history;
 
  d4est_util_fill_array(yp, 0., local_nodes);
  double* Au = P4EST_ALLOC(double, local_nodes); 
  d4est_elliptic_data_t mg_prob_vecs;
  d4est_elliptic_data_copy_ptrs(kct->vecs, &mg_prob_vecs);
  mg_prob_vecs.u = yp;
  mg_prob_vecs.rhs = xp;
  mg_prob_vecs.Au = Au;

  d4est_solver_multigrid_solve
    (
     kct->p4est,
     &mg_prob_vecs,
     kct->fcns,
     mg_data
    );
      
  P4EST_FREE(Au);
}


d4est_krylov_pc_t*
d4est_krylov_pc_multigrid_create
(
d4est_solver_multigrid_t* mg_data,
 void(*user_setup_fcn)(d4est_krylov_pc_t* kpc)
){

  d4est_krylov_pc_t* pc = P4EST_ALLOC(d4est_krylov_pc_t, 1);
  pc->pc_apply = d4est_krylov_pc_multigrid_apply;

  d4est_krylov_pc_multigrid_data_t* kpcmgdata = P4EST_ALLOC(d4est_krylov_pc_multigrid_data_t,1);
  kpcmgdata->mg_data = mg_data;
  kpcmgdata->ratio_is_getting_bad_counts = 0;
  mg_data->linear_operator_updates = 0;
  kpcmgdata->user_setup_fcn = user_setup_fcn;
  pc->pc_setup = d4est_krylov_pc_multigrid_setup;
  pc->pc_data = kpcmgdata;

  return pc;
}
