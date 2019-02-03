#include <pXest.h>
#include <d4est_krylov_pc_cheby.h>
#include <d4est_solver_multigrid_smoother_cheby.h>
#include <d4est_solver_cg_eigs.h>
#include <d4est_linalg.h>
#include <zlog.h>
#include <ini.h>

static
int d4est_krylov_pc_cheby_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_krylov_pc_cheby_data_t* pconfig = (d4est_krylov_pc_cheby_data_t*)user;
  const char* input_section = "d4est_krylov_pc_cheby";
  
  if (d4est_util_match_couple(section,input_section,name,"cheby_imax")) {
    D4EST_ASSERT(pconfig->cheby_imax == -1);
    pconfig->cheby_imax = atoi(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"cheby_reuse_eig")) {
    D4EST_ASSERT(pconfig->cheby_reuse_eig == -1);
    pconfig->cheby_reuse_eig = atoi(value);
  }  
  else if (d4est_util_match_couple(section,input_section,name,"cheby_eigs_cg_imax")) {
    D4EST_ASSERT(pconfig->cheby_eigs_cg_imax == -1);
    pconfig->cheby_eigs_cg_imax = atoi(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"cheby_eigs_lmax_lmin_ratio")) {
    D4EST_ASSERT(pconfig->cheby_eigs_lmax_lmin_ratio == -1);
    pconfig->cheby_eigs_lmax_lmin_ratio = atof(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"cheby_eigs_max_multiplier")) {
    D4EST_ASSERT(pconfig->cheby_eigs_max_multiplier == -1);
    pconfig->cheby_eigs_max_multiplier = atof(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"cheby_print_residual_norm")) {
    D4EST_ASSERT(pconfig->cheby_print_residual_norm == 0);
    pconfig->cheby_print_residual_norm = atoi(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"cheby_print_spectral_bound")) {
    D4EST_ASSERT(pconfig->cheby_print_spectral_bound == 0);
    pconfig->cheby_print_spectral_bound = atoi(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"cheby_print_spectral_bound_iterations")) {
    D4EST_ASSERT(pconfig->cheby_print_spectral_bound_iterations == 0);
    pconfig->cheby_print_spectral_bound_iterations = atoi(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"cheby_use_zero_guess_for_eigs")) {
    D4EST_ASSERT(pconfig->cheby_use_zero_guess_for_eigs == 0);
    pconfig->cheby_use_zero_guess_for_eigs = atoi(value);
  }
  
  else if (d4est_util_match_couple(section,input_section,name,"cheby_use_new_cg_eigs")) {
    D4EST_ASSERT(pconfig->cheby_use_new_cg_eigs == 0);
    pconfig->cheby_use_new_cg_eigs = atoi(value);
  }  
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


void
d4est_krylov_pc_cheby_setup
(
 d4est_krylov_pc_t* kpc
)
{
  d4est_krylov_pc_cheby_data_t* kpccheby = kpc->pc_data;
  
  if (kpccheby->user_setup_fcn != NULL)
    kpccheby->user_setup_fcn(kpc);
}

void
d4est_krylov_pc_cheby_destroy(d4est_krylov_pc_t* kpc){
  kpc->pc_apply = NULL;
  kpc->pc_setup = NULL;
  kpc->pc_ctx = NULL;
  P4EST_FREE(kpc->pc_data);
  kpc->pc_data = NULL;
  P4EST_FREE(kpc);
}

void
d4est_krylov_pc_cheby_apply
(
 d4est_krylov_pc_t* kpc,
 double* xp,
 double* yp
)
{
  d4est_krylov_pc_cheby_data_t* kpccheby = kpc->pc_data;
  krylov_ctx_t* kct = kpc->pc_ctx;
  int local_nodes = kct->vecs->local_nodes;
  double* residual_history;
 
  d4est_util_fill_array(yp, 0., local_nodes);
  double* Au = P4EST_ALLOC(double, local_nodes); 

  d4est_elliptic_data_t cheby_prob_vecs;
  d4est_elliptic_data_copy_ptrs(kct->vecs, &cheby_prob_vecs);
  
  cheby_prob_vecs.u = yp;
  cheby_prob_vecs.rhs = xp;
  cheby_prob_vecs.Au = Au;

  if(!kpccheby->cheby_reuse_eig
     || kpccheby->eig == -1){

    double* temp = NULL;
    double* zero_vec = NULL;
    if (kpccheby->cheby_use_zero_guess_for_eigs){
      zero_vec = P4EST_ALLOC_ZERO(double, cheby_prob_vecs.local_nodes);
      cheby_prob_vecs.u = zero_vec;
    }
    cg_eigs
      (
       kct->p4est,
       &cheby_prob_vecs,
       kct->fcns,
       *kct->ghost,
       *kct->ghost_data,
       kct->d4est_ops,
       kct->d4est_geom,
       kct->d4est_quad,
       kct->d4est_factors,
       kpccheby->cheby_eigs_cg_imax,
       kpccheby->cheby_print_spectral_bound_iterations,
       kpccheby->cheby_use_new_cg_eigs,
       &kpccheby->eig
      );

   if (kpccheby->cheby_use_zero_guess_for_eigs){
     P4EST_FREE(zero_vec);
     cheby_prob_vecs.u = yp;
   }
   
   kpccheby->eig *= kpccheby->cheby_eigs_max_multiplier;
  }
  
  double* r = P4EST_ALLOC(double, cheby_prob_vecs.local_nodes);
  d4est_solver_multigrid_smoother_cheby_iterate_aux
    (
     kct->p4est,
     kct->d4est_ops,
     kct->d4est_geom,
     kct->d4est_quad,
     kct->d4est_factors,
     *kct->ghost,
     *kct->ghost_data,
     &cheby_prob_vecs,
     kct->fcns,
     r,
     kpccheby->cheby_imax,
     kpccheby->eig/kpccheby->cheby_eigs_lmax_lmin_ratio,
     kpccheby->eig,
     kpccheby->cheby_print_residual_norm,
     -1,
     0
    );

  P4EST_FREE(Au);
  P4EST_FREE(r);
}


d4est_krylov_pc_t*
d4est_krylov_pc_cheby_create
(
 const char* input_file,
 void(*user_setup_fcn)(d4est_krylov_pc_t* kpc)
){

  d4est_krylov_pc_t* pc = P4EST_ALLOC(d4est_krylov_pc_t, 1);
  pc->pc_apply = d4est_krylov_pc_cheby_apply;
  d4est_krylov_pc_cheby_data_t* kpccheby = P4EST_ALLOC(d4est_krylov_pc_cheby_data_t,1);
  kpccheby->user_setup_fcn = user_setup_fcn;
  /* set externally in input file */
  kpccheby->cheby_imax = -1;
  kpccheby->cheby_eigs_cg_imax = -1;
  kpccheby->cheby_eigs_lmax_lmin_ratio = -1;
  kpccheby->cheby_eigs_max_multiplier = -1;
  kpccheby->cheby_print_residual_norm = 0;
  kpccheby->cheby_print_spectral_bound = 0;
  kpccheby->cheby_print_spectral_bound_iterations = 0;
  kpccheby->cheby_use_new_cg_eigs = 0;
  kpccheby->cheby_use_zero_guess_for_eigs = 0;
  kpccheby->cheby_reuse_eig = -1;
  kpccheby->eig = -1;
  
  if(
     ini_parse(input_file,
               d4est_krylov_pc_cheby_input_handler,
               kpccheby) < 0
  ){
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("d4est_krylov_pc_cheby", kpccheby->cheby_imax, -1);
  D4EST_CHECK_INPUT("d4est_krylov_pc_cheby", kpccheby->cheby_eigs_cg_imax, -1);
  D4EST_CHECK_INPUT("d4est_krylov_pc_cheby", kpccheby->cheby_eigs_lmax_lmin_ratio, -1);
  D4EST_CHECK_INPUT("d4est_krylov_pc_cheby", kpccheby->cheby_eigs_max_multiplier, -1);
  D4EST_CHECK_INPUT("d4est_krylov_pc_cheby", kpccheby->cheby_reuse_eig, -1);
  
  /* kpccheby->cheby = cheby; */
  kpccheby->user_setup_fcn = user_setup_fcn;
  pc->pc_setup = d4est_krylov_pc_cheby_setup;
  pc->pc_data = kpccheby;

  return pc;
}
