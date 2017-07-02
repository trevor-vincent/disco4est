#include <multigrid_bottom_solver_cheby_d4est.h>
#include <multigrid_smoother_cheby_d4est.h>
#include <d4est_linalg.h>
#include <ini.h>
#include <util.h>
#include <cg_eigs.h>

static
int
multigrid_bottom_solver_cheby_d4est_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{  
  multigrid_bottom_solver_cheby_d4est_t* pconfig = ((multigrid_bottom_solver_cheby_d4est_t*)user);

  if (util_match_couple(section,"mg_bottom_solver_cheby_d4est",name,"cheby_imax")) {
    mpi_assert(pconfig->cheby_imax == -1);
    pconfig->cheby_imax = atoi(value);
  }
  else if (util_match_couple(section,"mg_bottom_solver_cheby_d4est",name,"cheby_eigs_cg_imax")) {
    mpi_assert(pconfig->cheby_eigs_cg_imax == -1);
    pconfig->cheby_eigs_cg_imax = atoi(value);
  }
  else if (util_match_couple(section,"mg_bottom_solver_cheby_d4est",name,"cheby_eigs_lmax_lmin_ratio")) {
    mpi_assert(pconfig->cheby_eigs_lmax_lmin_ratio == -1);
    pconfig->cheby_eigs_lmax_lmin_ratio = atof(value);
  }
  else if (util_match_couple(section,"mg_bottom_solver_cheby_d4est",name,"cheby_eigs_max_multiplier")) {
    mpi_assert(pconfig->cheby_eigs_max_multiplier == -1);
    pconfig->cheby_eigs_max_multiplier = atof(value);
  }
  else if (util_match_couple(section,"mg_bottom_solver_cheby_d4est",name,"cheby_print_residual_norm")) {
    mpi_assert(pconfig->cheby_print_residual_norm == -1);
    pconfig->cheby_print_residual_norm = atoi(value);
  }
  else if (util_match_couple(section,"mg_bottom_solver_cheby_d4est",name,"cheby_print_eig")) {
    mpi_assert(pconfig->cheby_print_eig == -1);
    pconfig->cheby_print_eig = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


static void
multigrid_bottom_solver_cheby_d4est
(
 p4est_t* p4est,
 d4est_elliptic_problem_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 double* r
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_bottom_solver_cheby_d4est_t* cheby = mg_data->bottom_solver->user;
  multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  int level = 0;
  
  cg_eigs
    (
     p4est,
     vecs,
     fcns,
     *(updater->ghost),
     *(updater->ghost_data),
     mg_data->d4est_ops,
     updater->d4est_geom,
     cheby->cheby_eigs_cg_imax,
     &cheby->eig
    );

  cheby->eig *= cheby->cheby_eigs_max_multiplier;

  int iter = cheby->cheby_imax;
  double lmin = cheby->eig/cheby->cheby_eigs_lmax_lmin_ratio;
  double lmax = cheby->eig;

  if (cheby->cheby_print_eig){
    printf("[MG_BOTTOM_SOLVER_CHEBY]: Lev %d Max_eig %.25f\n", level, cheby->eig);
  }
  
  multigrid_smoother_cheby_d4est_iterate
    (
     p4est,
     vecs,
     fcns,
     r,
     iter,
     lmin,
     lmax,
     cheby->cheby_print_residual_norm
    );
}

void
multigrid_bottom_solver_cheby_d4est_destroy
(
 multigrid_bottom_solver_t* bottom_solver
)
{
  P4EST_FREE(bottom_solver->user);
  bottom_solver->solve = NULL;
  bottom_solver->update = NULL;
  P4EST_FREE(bottom_solver);
}

multigrid_bottom_solver_t*
multigrid_bottom_solver_cheby_d4est_init
(
 p4est_t* p4est,
 int num_of_levels,
 const char* input_file
)
{
  multigrid_bottom_solver_t* bottom_solver = P4EST_ALLOC(multigrid_bottom_solver_t, 1);
  multigrid_bottom_solver_cheby_d4est_t* cheby_data = P4EST_ALLOC(multigrid_bottom_solver_cheby_d4est_t, 1);
  
  /* set externally in input file */
  cheby_data->cheby_imax = -1;
  cheby_data->cheby_eigs_cg_imax = -1;
  cheby_data->cheby_eigs_lmax_lmin_ratio = -1;
  cheby_data->cheby_eigs_max_multiplier = -1;
  cheby_data->cheby_print_residual_norm = -1;
  cheby_data->cheby_print_eig = -1;
  
  if (ini_parse(input_file, multigrid_bottom_solver_cheby_d4est_input_handler, cheby_data) < 0) {
    mpi_abort("Can't load input file");
  }

  D4EST_CHECK_INPUT("multigrid_bottom_solver", cheby_data->cheby_imax, -1);
  D4EST_CHECK_INPUT("multigrid_bottom_solver", cheby_data->cheby_eigs_cg_imax, -1);
  D4EST_CHECK_INPUT("multigrid_bottom_solver", cheby_data->cheby_eigs_lmax_lmin_ratio, -1);
  D4EST_CHECK_INPUT("multigrid_bottom_solver", cheby_data->cheby_eigs_max_multiplier, -1);
  D4EST_CHECK_INPUT("multigrid_bottom_solver", cheby_data->cheby_print_residual_norm, -1);
  D4EST_CHECK_INPUT("multigrid_bottom_solver", cheby_data->cheby_print_eig, -1);
  
  if(p4est->mpirank == 0){
    printf("[D4EST_INFO]: Multigrid Bottom Solver Cheby Parameters\n");
    printf("[D4EST_INFO]: Cheby imax = %d\n", cheby_data->cheby_imax);
    printf("[D4EST_INFO]: Cheby eigs cg max = %d\n", cheby_data->cheby_eigs_cg_imax);
    printf("[D4EST_INFO]: Cheby eigs lmax_lmin_ratio = %f\n", cheby_data->cheby_eigs_lmax_lmin_ratio);
    printf("[D4EST_INFO]: Cheby eigs max multiplier = %.25f\n", cheby_data->cheby_eigs_max_multiplier);
    printf("[D4EST_INFO]: Cheby print residual norm = %d\n", cheby_data->cheby_print_residual_norm);
    printf("[D4EST_INFO]: Cheby print eig = %d\n", cheby_data->cheby_print_eig);
  }

  bottom_solver->user = cheby_data;
  bottom_solver->solve = multigrid_bottom_solver_cheby_d4est;
  bottom_solver->update = NULL;

  return bottom_solver;
}
