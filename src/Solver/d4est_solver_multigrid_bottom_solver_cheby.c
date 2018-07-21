#include <d4est_solver_multigrid_bottom_solver_cheby.h>
#include <d4est_solver_multigrid_smoother_cheby.h>
#include <d4est_linalg.h>
#include <ini.h>
#include <d4est_util.h>
#include <d4est_solver_cg_eigs.h>

static
int
d4est_solver_multigrid_bottom_solver_cheby_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{  
  d4est_solver_multigrid_bottom_solver_cheby_t* pconfig = ((d4est_solver_multigrid_bottom_solver_cheby_t*)user);

  if (d4est_util_match_couple(section,"mg_bottom_solver_cheby",name,"cheby_imax")) {
    D4EST_ASSERT(pconfig->cheby_imax == -1);
    pconfig->cheby_imax = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_bottom_solver_cheby",name,"cheby_eigs_cg_imax")) {
    D4EST_ASSERT(pconfig->cheby_eigs_cg_imax == -1);
    pconfig->cheby_eigs_cg_imax = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_bottom_solver_cheby",name,"cheby_eigs_lmax_lmin_ratio")) {
    D4EST_ASSERT(pconfig->cheby_eigs_lmax_lmin_ratio == -1);
    pconfig->cheby_eigs_lmax_lmin_ratio = atof(value);
  }
  else if (d4est_util_match_couple(section,"mg_bottom_solver_cheby",name,"cheby_eigs_max_multiplier")) {
    D4EST_ASSERT(pconfig->cheby_eigs_max_multiplier == -1);
    pconfig->cheby_eigs_max_multiplier = atof(value);
  }
  else if (d4est_util_match_couple(section,"mg_bottom_solver_cheby",name,"cheby_print_residual_norm")) {
    D4EST_ASSERT(pconfig->cheby_print_residual_norm == 0);
    pconfig->cheby_print_residual_norm = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_bottom_solver_cheby",name,"cheby_print_spectral_bound")) {
    D4EST_ASSERT(pconfig->cheby_print_spectral_bound == 0);
    pconfig->cheby_print_spectral_bound = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_bottom_solver_cheby",name,"cheby_print_spectral_bound_iterations")) {
    D4EST_ASSERT(pconfig->cheby_print_spectral_bound_iterations == 0);
    pconfig->cheby_print_spectral_bound_iterations = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


static void
d4est_solver_multigrid_bottom_solver_cheby
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 double* r
)
{
  d4est_solver_multigrid_data_t* mg_data = p4est->user_pointer;
  d4est_solver_multigrid_bottom_solver_cheby_t* cheby = mg_data->bottom_solver->user;
  d4est_solver_multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  int level = 0;
  
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
     &cheby->eig
    );

  cheby->eig *= cheby->cheby_eigs_max_multiplier;

  int iter = cheby->cheby_imax;
  double lmin = cheby->eig/cheby->cheby_eigs_lmax_lmin_ratio;
  double lmax = cheby->eig;

  if (cheby->cheby_print_spectral_bound){
    zlog_category_t *c_default = zlog_get_category("d4est_d4est_solver_multigrid_bottom_solver_cheby");    
    zlog_info(c_default, "Lev %d Max_eig %f Multiplier %f\n", level, cheby->eig, cheby->cheby_eigs_max_multiplier);
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

void
d4est_solver_multigrid_bottom_solver_cheby_destroy
(
 d4est_solver_multigrid_bottom_solver_t* bottom_solver
)
{
  P4EST_FREE(bottom_solver->user);
  bottom_solver->solve = NULL;
  bottom_solver->update = NULL;
  P4EST_FREE(bottom_solver);
}

d4est_solver_multigrid_bottom_solver_t*
d4est_solver_multigrid_bottom_solver_cheby_init
(
 p4est_t* p4est,
 int num_of_levels,
 const char* input_file
)
{
  d4est_solver_multigrid_bottom_solver_t* bottom_solver = P4EST_ALLOC(d4est_solver_multigrid_bottom_solver_t, 1);
  d4est_solver_multigrid_bottom_solver_cheby_t* cheby_data = P4EST_ALLOC(d4est_solver_multigrid_bottom_solver_cheby_t, 1);
  
  /* set externally in input file */
  cheby_data->cheby_imax = -1;
  cheby_data->cheby_eigs_cg_imax = -1;
  cheby_data->cheby_eigs_lmax_lmin_ratio = -1;
  cheby_data->cheby_eigs_max_multiplier = -1;
  cheby_data->cheby_print_residual_norm = 0;
  cheby_data->cheby_print_spectral_bound = 0;
  cheby_data->cheby_print_spectral_bound_iterations = 0;
  
  if (ini_parse(input_file, d4est_solver_multigrid_bottom_solver_cheby_input_handler, cheby_data) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("mg_bottom_solver_cheby", cheby_data->cheby_imax, -1);
  D4EST_CHECK_INPUT("mg_bottom_solver_cheby", cheby_data->cheby_eigs_cg_imax, -1);
  D4EST_CHECK_INPUT("mg_bottom_solver_cheby", cheby_data->cheby_eigs_lmax_lmin_ratio, -1);
  D4EST_CHECK_INPUT("mg_bottom_solver_cheby", cheby_data->cheby_eigs_max_multiplier, -1);
  /* D4EST_CHECK_INPUT("mg_bottom_solver_cheby", cheby_data->cheby_print_residual_norm, -1); */
  /* D4EST_CHECK_INPUT("mg_bottom_solver_cheby", cheby_data->cheby_print_spectral_bound, -1); */
  /* D4EST_CHECK_INPUT("mg_bottom_solver_cheby", cheby_data->cheby_print_spectral_bound_iterations, -1); */
  
  if(p4est->mpirank == 0){
    printf("[D4EST_INFO]: D4est_Solver_Multigrid Bottom Solver Cheby Parameters\n");
    printf("[D4EST_INFO]: Cheby imax = %d\n", cheby_data->cheby_imax);
    printf("[D4EST_INFO]: Cheby eigs cg max = %d\n", cheby_data->cheby_eigs_cg_imax);
    printf("[D4EST_INFO]: Cheby eigs lmax_lmin_ratio = %f\n", cheby_data->cheby_eigs_lmax_lmin_ratio);
    printf("[D4EST_INFO]: Cheby eigs max multiplier = %.25f\n", cheby_data->cheby_eigs_max_multiplier);
    printf("[D4EST_INFO]: Cheby print residual norm = %d\n", cheby_data->cheby_print_residual_norm);
    printf("[D4EST_INFO]: Cheby print spectral bound = %d\n", cheby_data->cheby_print_spectral_bound);
    printf("[D4EST_INFO]: Cheby print spectral bound iterations = %d\n", cheby_data->cheby_print_spectral_bound_iterations);
  }

  bottom_solver->user = cheby_data;
  bottom_solver->solve = d4est_solver_multigrid_bottom_solver_cheby;
  bottom_solver->update = NULL;

  return bottom_solver;
}
