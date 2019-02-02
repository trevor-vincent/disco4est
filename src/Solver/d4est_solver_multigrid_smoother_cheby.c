#include <pXest.h>
#include <d4est_solver_multigrid_smoother_cheby.h>
#include <d4est_linalg.h>
#include <ini.h>
#include <d4est_util.h>
#include <d4est_solver_cg_eigs.h>

static
int
d4est_solver_multigrid_smoother_cheby_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{  
  d4est_solver_multigrid_smoother_cheby_t* pconfig = ((d4est_solver_multigrid_smoother_cheby_t*)user);

  if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_imax")) {
    D4EST_ASSERT(pconfig->cheby_imax == -1);
    pconfig->cheby_imax = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_eigs_cg_imax")) {
    D4EST_ASSERT(pconfig->cheby_eigs_cg_imax == -1);
    pconfig->cheby_eigs_cg_imax = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_eigs_lmax_lmin_ratio")) {
    D4EST_ASSERT(pconfig->cheby_eigs_lmax_lmin_ratio == -1);
    pconfig->cheby_eigs_lmax_lmin_ratio = atof(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_eigs_max_multiplier")) {
    D4EST_ASSERT(pconfig->cheby_eigs_max_multiplier == -1);
    pconfig->cheby_eigs_max_multiplier = atof(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_eigs_reuse_fromdownvcycle")) {
    D4EST_ASSERT(pconfig->cheby_eigs_reuse_fromdownvcycle == -1);
    pconfig->cheby_eigs_reuse_fromdownvcycle = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_eigs_reuse_fromlastvcycle")) {
    D4EST_ASSERT(pconfig->cheby_eigs_reuse_fromlastvcycle == -1);
    pconfig->cheby_eigs_reuse_fromlastvcycle = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_print_residual_norm")) {
    D4EST_ASSERT(pconfig->cheby_print_residual_norm == 0);
    pconfig->cheby_print_residual_norm = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_print_spectral_bound")) {
    D4EST_ASSERT(pconfig->cheby_print_spectral_bound == 0);
    pconfig->cheby_print_spectral_bound = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_print_spectral_bound_iterations")) {
    D4EST_ASSERT(pconfig->cheby_print_spectral_bound_iterations == 0);
    pconfig->cheby_print_spectral_bound_iterations = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"make_cheby_symmetric")) {
    D4EST_ASSERT(pconfig->make_cheby_symmetric == 0);
    pconfig->make_cheby_symmetric = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_use_new_cg_eigs")) {
    D4EST_ASSERT(pconfig->cheby_use_new_cg_eigs == 0);
    pconfig->cheby_use_new_cg_eigs = atoi(value);
  }
  
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void
d4est_solver_multigrid_smoother_cheby_destroy(d4est_solver_multigrid_smoother_t* smoother)
{
  d4est_solver_multigrid_smoother_cheby_t* cheby = smoother->user;
  P4EST_FREE(cheby->eigs);
  P4EST_FREE(cheby);
  P4EST_FREE(smoother);
}

void 
d4est_solver_multigrid_smoother_cheby_iterate
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 double* r,
 int iter,
 double lmin,
 double lmax,
 int print_residual_norm,
 int mg_level
 /* d4est_solver_multigrid_cheby_params_t* cheby_params */
)
{
  d4est_solver_multigrid_t* mg_data = (d4est_solver_multigrid_t*) p4est->user_pointer;
  d4est_solver_multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  zlog_category_t *c_default = zlog_get_category("d4est_solver_multigrid_smoother_cheby");   
  d4est_ghost_t* d4est_ghost = updater->current_d4est_ghost;
  d4est_ghost_data_t* d4est_ghost_data = updater->current_d4est_ghost_data;
  
  int i;
  double d = (lmax + lmin)*.5;
  double c = (lmax - lmin)*.5;

  int local_nodes = vecs->local_nodes;
  double* Au = vecs->Au;
  double* u = vecs->u;
  double* rhs = vecs->rhs;
  
  double* p;
  double alpha,beta;
  p = P4EST_ALLOC(double, local_nodes);
  
  d4est_util_fill_array(p, 0., local_nodes);
  for (i = 0; i < iter; i++){
    /* calculate residual r = rhs - Au */

  d4est_elliptic_eqns_apply_lhs
    (
     p4est,
     d4est_ghost,
     d4est_ghost_data,
     fcns,
     vecs,
     mg_data->d4est_ops,
     mg_data->d4est_geom,
     mg_data->d4est_quad,
     updater->current_d4est_factors
    );
  
    d4est_util_copy_1st_to_2nd(Au, r, local_nodes);
    d4est_linalg_vec_xpby(rhs, -1., r, local_nodes);

    if(print_residual_norm && p4est->mpirank == 0){
      zlog_info(c_default, "mg_level, iter, residual = %d,%d,%.25f", mg_level, i, d4est_linalg_vec_dot(r,r,local_nodes));
    }
    
    if (i == 0)
      alpha = 1./d;
    else if (i == 1)
      alpha = 2.*d/(2*d*d - c*c);
    else
      alpha = 1./(d-(alpha*c*c/4.));

    beta = alpha*d - 1.;
   
    d4est_linalg_vec_scale(alpha,r,local_nodes);
    d4est_linalg_vec_xpby(&r[0], beta, &p[0], local_nodes);   
    d4est_linalg_vec_axpy(1., p, u, local_nodes);
  }

  /* calculate the residual */
  d4est_elliptic_eqns_apply_lhs
    (
     p4est,
     d4est_ghost,
     d4est_ghost_data,
     fcns,
     vecs,
     mg_data->d4est_ops,
     mg_data->d4est_geom,
     mg_data->d4est_quad,
     updater->current_d4est_factors
    );
  
  d4est_util_copy_1st_to_2nd(Au, r, local_nodes);
  d4est_linalg_vec_xpby(rhs, -1., r, local_nodes);

  P4EST_FREE(p);
}

static void
d4est_solver_multigrid_smoother_cheby_update
(
 p4est_t* p4est,
 int level,
 d4est_elliptic_data_t* vecs
)
{
  d4est_solver_multigrid_t* mg_data = (d4est_solver_multigrid_t*) p4est->user_pointer;
  d4est_solver_multigrid_smoother_cheby_t* cheby = mg_data->smoother->user;  
  int vcycle = mg_data->vcycle_num_finished;

  if(mg_data->mg_state == PRE_V){
    if(
       cheby->cheby_eigs_reuse_fromlastvcycle == 1 &&
       vcycle != 0
    ){
      cheby->cheby_eigs_compute = 0;
    }
    else {
      cheby->cheby_eigs_compute = 1;
    }
  }
  
  else if (mg_data->mg_state == UPV_PRE_SMOOTH){
    if (
        cheby->cheby_eigs_reuse_fromdownvcycle == 1
        || (cheby->cheby_eigs_reuse_fromlastvcycle == 1 &&
            vcycle != 0)
    ){
      cheby->cheby_eigs_compute = 0;
    }
    else {
      cheby->cheby_eigs_compute = 1;
    }
  }
    
  else {
    return;
  }
}

static void
d4est_solver_multigrid_smoother_cheby
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 double* r,
 int level
)
{  
  d4est_solver_multigrid_t* mg_data = p4est->user_pointer;
  d4est_solver_multigrid_smoother_cheby_t* cheby = mg_data->smoother->user;
  d4est_solver_multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;

  double* zero_vector = NULL;
  double* temp_vector = NULL;
  if (cheby->cheby_eigs_compute){

    if(cheby->make_cheby_symmetric){
      temp_vector = vecs->u;
      zero_vector = P4EST_ALLOC_ZERO(double, vecs->local_nodes);
      vecs->u = zero_vector;
    }
    
      cg_eigs
        (
         p4est,
         vecs,
         fcns,
         updater->current_d4est_ghost,
         updater->current_d4est_ghost_data,
         mg_data->d4est_ops,
         mg_data->d4est_geom,
         mg_data->d4est_quad,
         updater->current_d4est_factors,
         cheby->cheby_eigs_cg_imax,
         cheby->cheby_print_spectral_bound_iterations,
         cheby->cheby_use_new_cg_eigs,
         &cheby->eigs[level]
        );

      if(cheby->make_cheby_symmetric){
        vecs->u = temp_vector;
        P4EST_FREE(zero_vector);
      }
      
      cheby->eigs[level] *= cheby->cheby_eigs_max_multiplier;
  }

  if(cheby->make_cheby_symmetric == 1){
    if (cheby->cheby_eigs_reuse_fromdownvcycle != 1){
      printf("If you set make_cheby_symmetric == 1, please set cheby_eigs_reuse_fromdownvcycle = 1\n");
      D4EST_ABORT("");
    }
  }
  
  if (cheby->cheby_eigs_compute == 0 &&
      cheby->make_cheby_symmetric == 1){
    double dummy;

    if(cheby->make_cheby_symmetric){
      temp_vector = vecs->u;
      zero_vector = P4EST_ALLOC_ZERO(double, vecs->local_nodes);
      vecs->u = zero_vector;
    }

    cg_eigs
        (
         p4est,
         vecs,
         fcns,
         updater->current_d4est_ghost,
         updater->current_d4est_ghost_data,
         mg_data->d4est_ops,
         mg_data->d4est_geom,
         mg_data->d4est_quad,
         updater->current_d4est_factors,
         cheby->cheby_eigs_cg_imax,
         cheby->cheby_print_spectral_bound_iterations,
         cheby->cheby_use_new_cg_eigs,
         &dummy
        );

      if(cheby->make_cheby_symmetric){
        vecs->u = temp_vector;
        P4EST_FREE(zero_vector);
      }

      /* cheby->eigs[level] *= cheby->cheby_eigs_max_multiplier; */
  }
  
  int iter = cheby->cheby_imax;
  double lmin = cheby->eigs[level]/cheby->cheby_eigs_lmax_lmin_ratio;
  double lmax = cheby->eigs[level];

  if (cheby->cheby_print_spectral_bound && p4est->mpirank == 0){
    zlog_category_t *c_default = zlog_get_category("d4est_solver_multigrid_smoother_cheby");    
    zlog_info(c_default, "mg_level,spectral_bound,multiplier = %d,%f,%f", level, cheby->eigs[level], cheby->cheby_eigs_max_multiplier);
  }
  
  d4est_solver_multigrid_smoother_cheby_iterate
    (
     p4est,
     vecs,
     fcns,
     r,
     iter,
     lmin,
     lmax,
     cheby->cheby_print_residual_norm,
     level
    );
}



d4est_solver_multigrid_smoother_t*
d4est_solver_multigrid_smoother_cheby_init
(
 p4est_t* p4est,
 int num_of_levels,
 const char* input_file
)
{
  d4est_solver_multigrid_smoother_t* smoother = P4EST_ALLOC(d4est_solver_multigrid_smoother_t, 1);
  d4est_solver_multigrid_smoother_cheby_t* cheby_data = P4EST_ALLOC(d4est_solver_multigrid_smoother_cheby_t, 1);
  
  cheby_data->eigs = P4EST_ALLOC(double, num_of_levels);  
  
  /* set externally in input file */
  cheby_data->cheby_imax = -1;
  cheby_data->cheby_eigs_cg_imax = -1;
  cheby_data->cheby_eigs_lmax_lmin_ratio = -1;
  cheby_data->cheby_eigs_max_multiplier = -1;
  cheby_data->cheby_eigs_reuse_fromdownvcycle = -1;
  cheby_data->cheby_eigs_reuse_fromlastvcycle = -1;
  cheby_data->cheby_print_residual_norm = 0;
  cheby_data->cheby_print_spectral_bound = 0;
  cheby_data->cheby_print_spectral_bound_iterations = 0;
  cheby_data->cheby_use_new_cg_eigs = 0;
  cheby_data->make_cheby_symmetric = 0;

  /* set internally */
  cheby_data->mpirank = p4est->mpirank;
  cheby_data->cheby_eigs_compute = -1;
  
  if (ini_parse(input_file, d4est_solver_multigrid_smoother_cheby_input_handler, cheby_data) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  if (cheby_data->make_cheby_symmetric){
    if(!cheby_data->cheby_eigs_reuse_fromdownvcycle){
      D4EST_ABORT("Do not use make_cheby_symmetric without cheby_eigs_reuse_fromdownvycle");
    }
  }
  
  D4EST_CHECK_INPUT("mg_smoother_cheby", cheby_data->cheby_imax, -1);
  D4EST_CHECK_INPUT("mg_smoother_cheby", cheby_data->cheby_eigs_cg_imax, -1);
  D4EST_CHECK_INPUT("mg_smoother_cheby", cheby_data->cheby_eigs_lmax_lmin_ratio, -1);
  D4EST_CHECK_INPUT("mg_smoother_cheby", cheby_data->cheby_eigs_max_multiplier, -1);
  D4EST_CHECK_INPUT("mg_smoother_cheby", cheby_data->cheby_eigs_reuse_fromdownvcycle, -1);
  D4EST_CHECK_INPUT("mg_smoother_cheby", cheby_data->cheby_eigs_reuse_fromlastvcycle, -1);
  /* D4EST_CHECK_INPUT("mg_smoother_cheby", cheby_data->cheby_use_new_cg_eigs, -1); */
  /* D4EST_CHECK_INPUT("mg_smoother_cheby", cheby_data->cheby_print_residual_norm, -1); */
  /* D4EST_CHECK_INPUT("mg_smoother_cheby", cheby_data->cheby_print_spectral_bound, -1); */
  /* D4EST_CHECK_INPUT("mg_smoother_cheby", cheby_data->cheby_print_spectral_bound_iterations, -1); */
  
  
  if(p4est->mpirank == 0){
    zlog_category_t *c_default = zlog_get_category("d4est_solver_multigrid_smoother_cheby");
    zlog_debug(c_default, " Smoother imax = %d", cheby_data->cheby_imax);
    zlog_debug(c_default, " Smoother eigs cg max = %d", cheby_data->cheby_eigs_cg_imax);
    zlog_debug(c_default, " Smoother eigs lmax_lmin_ratio = %f", cheby_data->cheby_eigs_lmax_lmin_ratio);
    zlog_debug(c_default, " Smoother eigs max multiplier = %.25f", cheby_data->cheby_eigs_max_multiplier);
    zlog_debug(c_default, " Smoother eigs reuse up vcycle = %d", cheby_data->cheby_eigs_reuse_fromdownvcycle);
    zlog_debug(c_default, " Smoother eigs reuse from last vcycle = %d", cheby_data->cheby_eigs_reuse_fromlastvcycle);
    zlog_debug(c_default, " Smoother print residual norm = %d", cheby_data->cheby_print_residual_norm);
    zlog_debug(c_default, " Smoother print eigs = %d", cheby_data->cheby_print_spectral_bound);
    zlog_debug(c_default, " Smoother use new cg eigs scheme = %d", cheby_data->cheby_use_new_cg_eigs);
    zlog_debug(c_default, " Smoother make_cheby_symmetric = %d", cheby_data->make_cheby_symmetric);
  }

  smoother->user = cheby_data;
  smoother->smooth = d4est_solver_multigrid_smoother_cheby;
  smoother->update = d4est_solver_multigrid_smoother_cheby_update;

  return smoother;
}


