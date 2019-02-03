#include <pXest.h>
#include <d4est_krylov_pc_schwarz.h>
#include <d4est_solver_schwarz.h>
#include <d4est_linalg.h>
#include <zlog.h>
#include <ini.h>

static
int d4est_krylov_pc_schwarz_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_krylov_pc_schwarz_data_t* pconfig = (d4est_krylov_pc_schwarz_data_t*)user;
  const char* input_section = "schwarz_pc";
  
  if (
      d4est_util_match_couple(section,input_section,name,"schwarz_pc_iterations")) {
    D4EST_ASSERT(pconfig->iterations == -1);
    pconfig->iterations = atoi(value);
  }
  else if (
           d4est_util_match_couple(section,input_section,name,"schwarz_pc_verbose")) {
    D4EST_ASSERT(pconfig->verbose == -1);
    pconfig->verbose = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


void
d4est_krylov_pc_schwarz_setup
(
 d4est_krylov_pc_t* kpc
)
{
  d4est_krylov_pc_schwarz_data_t* kpcschwarz = kpc->pc_data;
  
  if (kpcschwarz->user_setup_fcn != NULL)
    kpcschwarz->user_setup_fcn(kpc);
}

void
d4est_krylov_pc_schwarz_destroy(d4est_krylov_pc_t* kpc){
  kpc->pc_apply = NULL;
  kpc->pc_setup = NULL;
  kpc->pc_ctx = NULL;
  P4EST_FREE(kpc->pc_data);
  kpc->pc_data = NULL;
  P4EST_FREE(kpc);
}

void
d4est_krylov_pc_schwarz_apply(d4est_krylov_pc_t* kpc, double* xp, double* yp)
{
  d4est_krylov_pc_schwarz_data_t* kpcschwarz = kpc->pc_data;
  d4est_solver_schwarz_t* schwarz = kpcschwarz->schwarz;
  krylov_ctx_t* kct = kpc->pc_ctx;
  int local_nodes = kct->vecs->local_nodes;
  int na;
  double* residual_history;
 
  d4est_util_fill_array(yp, 0., local_nodes);
  double* Au = P4EST_ALLOC(double, local_nodes); 
  double* residual = P4EST_ALLOC(double, local_nodes); 

  d4est_elliptic_data_t schwarz_prob_vecs;
  d4est_elliptic_data_copy_ptrs(kct->vecs, &schwarz_prob_vecs);
  
  schwarz_prob_vecs.u = yp;
  schwarz_prob_vecs.rhs = xp;
  schwarz_prob_vecs.Au = Au;

  for (int i = 0; i < kpcschwarz->iterations; i++){
  
    d4est_elliptic_eqns_apply_lhs
      (
       kct->p4est,
       *kct->ghost,
       *kct->ghost_data,
       kct->fcns,
       &schwarz_prob_vecs,
       kct->d4est_ops,
       kct->d4est_geom,
       kct->d4est_quad,
       kct->d4est_factors
      );
  

    d4est_linalg_vec_axpyeqz(-1., schwarz_prob_vecs.Au, schwarz_prob_vecs.rhs, residual, schwarz_prob_vecs.local_nodes);

    d4est_solver_schwarz_iterate
      (
       kct->p4est,
       kct->d4est_geom,
       kct->d4est_quad,
       kct->d4est_factors,
       *kct->ghost,
       schwarz,
       schwarz_prob_vecs.u,
       residual
      );

  }
       
  P4EST_FREE(residual);
  P4EST_FREE(Au);
}


d4est_krylov_pc_t*
d4est_krylov_pc_schwarz_create
(
 const char* input_file,
 d4est_solver_schwarz_t* schwarz,
 void(*user_setup_fcn)(d4est_krylov_pc_t* kpc)
){

  d4est_krylov_pc_t* pc = P4EST_ALLOC(d4est_krylov_pc_t, 1);
  pc->pc_apply = d4est_krylov_pc_schwarz_apply;
  d4est_krylov_pc_schwarz_data_t* kpcschwarz = P4EST_ALLOC(d4est_krylov_pc_schwarz_data_t,1);

  kpcschwarz->iterations = -1;
  kpcschwarz->verbose = -1;
  
  if(
     ini_parse(input_file,
               d4est_krylov_pc_schwarz_input_handler,
               kpcschwarz) < 0
  ){
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("schwarz_pc", kpcschwarz->iterations, -1);
  D4EST_CHECK_INPUT("schwarz_pc", kpcschwarz->verbose, -1);
  
  kpcschwarz->schwarz = schwarz;
  kpcschwarz->user_setup_fcn = user_setup_fcn;
  pc->pc_setup = d4est_krylov_pc_schwarz_setup;
  pc->pc_data = kpcschwarz;

  return pc;
}
