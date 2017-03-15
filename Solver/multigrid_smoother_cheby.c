#include <multigrid_smoother_cheby.h>
#include <linalg.h>
#include <ini.h>
#include <util.h>
#include <cg_eigs.h>

static
int
multigrid_smoother_cheby_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{  
  multigrid_smoother_cheby_t* pconfig = ((multigrid_smoother_cheby_t*)user);

  if (util_match_couple(section,"multigrid",name,"smoother_imax")) {
    mpi_assert(pconfig->smoother_imax == -1);
    pconfig->smoother_imax = atoi(value);
  }
  else if (util_match_couple(section,"multigrid",name,"smoother_eigs_cg_imax")) {
    mpi_assert(pconfig->smoother_eigs_cg_imax == -1);
    pconfig->smoother_eigs_cg_imax = atoi(value);
  }
  else if (util_match_couple(section,"multigrid",name,"smoother_eigs_lmax_lmin_ratio")) {
    mpi_assert(pconfig->smoother_eigs_lmax_lmin_ratio == -1);
    pconfig->smoother_eigs_lmax_lmin_ratio = atof(value);
  }
  else if (util_match_couple(section,"multigrid",name,"smoother_eigs_max_multiplier")) {
    mpi_assert(pconfig->smoother_eigs_max_multiplier == -1);
    pconfig->smoother_eigs_max_multiplier = atof(value);
  }
  else if (util_match_couple(section,"multigrid",name,"smoother_eigs_reuse_fromdownvcycle")) {
    mpi_assert(pconfig->smoother_eigs_reuse_fromdownvcycle == -1);
    pconfig->smoother_eigs_reuse_fromdownvcycle = atoi(value);
  }
  else if (util_match_couple(section,"multigrid",name,"smoother_eigs_reuse_fromlastvcycle")) {
    mpi_assert(pconfig->smoother_eigs_reuse_fromlastvcycle == -1);
    pconfig->smoother_eigs_reuse_fromlastvcycle = atoi(value);
  }
  else if (util_match_couple(section,"multigrid",name,"smoother_print_residual_norm")) {
    mpi_assert(pconfig->smoother_print_residual_norm == -1);
    pconfig->smoother_print_residual_norm = atoi(value);
  }
  else if (util_match_couple(section,"multigrid",name,"smoother_print_eigs")) {
    mpi_assert(pconfig->smoother_print_eigs == -1);
    pconfig->smoother_print_eigs = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void
multigrid_smoother_cheby_destroy(multigrid_smoother_t* smoother)
{

  multigrid_smoother_cheby_t* cheby = smoother->user;
  P4EST_FREE(cheby->eigs);
  P4EST_FREE(cheby);
  P4EST_FREE(smoother);
}


static void 
multigrid_smoother_cheby_iterate
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 double* r,
 int iter,
 double lmin,
 double lmax,
 int print_residual_norm
 /* multigrid_cheby_params_t* cheby_params */
)
{
  /* const int iter = cheby_params->iter; */
  /* const double lmin = cheby_params->lmin; */
  /* const double lmax = cheby_params->lmax; */

  multigrid_data_t* mg_data = (multigrid_data_t*) p4est->user_pointer;
  dgmath_jit_dbase_t* dgmath_jit_dbase = mg_data->dgmath_jit_dbase;
 
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
  
  linalg_fill_vec(p, 0., local_nodes);
  for (i = 0; i < iter; i++){
    /* calculate residual r = rhs - Au */
    fcns->apply_lhs(p4est, ghost, ghost_data, vecs, dgmath_jit_dbase);
    linalg_copy_1st_to_2nd(Au, r, local_nodes);
    linalg_vec_xpby(rhs, -1., r, local_nodes);

    if(print_residual_norm && p4est->mpirank == 0){
      printf("[MG_CHEBY_SMOOTHER]: iter, residual = %d, %.25f\n", i, linalg_vec_dot(r,r,local_nodes));
    }
    
    if (i == 0)
      alpha = 1./d;
    else if (i == 1)
      alpha = 2.*d/(2*d*d - c*c);
    else
      alpha = 1./(d-(alpha*c*c/4.));

    beta = alpha*d - 1.;
   
    linalg_vec_scale(alpha,r,local_nodes);
    linalg_vec_xpby(&r[0], beta, &p[0], local_nodes);   
    linalg_vec_axpy(1., p, u, local_nodes);
  }

  /* calculate the residual */
  fcns->apply_lhs(p4est, ghost, ghost_data, vecs, dgmath_jit_dbase);
  linalg_copy_1st_to_2nd(Au, r, local_nodes);
  linalg_vec_xpby(rhs, -1., r, local_nodes);

  P4EST_FREE(p);
  /* P4EST_FREE(ghost_data); */
  /* p4est_ghost_destroy (ghost); */
}

static void
multigrid_smoother_cheby_update
(
 p4est_t* p4est,
 int level,
 problem_data_t* vecs
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_smoother_cheby_t* cheby = mg_data->smoother->user;
    
  int vcycle = mg_data->vcycle_num_finished;

  if(mg_data->mg_state == PRE_V){
    if(
       cheby->smoother_eigs_reuse_fromlastvcycle == 1 &&
       vcycle != 0
    ){
      cheby->smoother_eigs_compute == 0;
    }
    else {
      cheby->smoother_eigs_compute == 1;
    }
  }
  
  else if (mg_data->mg_state == UPV_PRE_SMOOTH){
    if (
        cheby->smoother_eigs_reuse_fromdownvcycle == 1
    ){
      cheby->smoother_eigs_compute == 0;
    }
    else {
      cheby->smoother_eigs_compute == 1;
    }
  }
  
  else {
    return;
  }
}
 

static void
multigrid_smoother_cheby
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 double* r,
 int level
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_smoother_cheby_t* cheby = mg_data->smoother->user;
  
  if (cheby->smoother_eigs_compute){
      cg_eigs
        (
         p4est,
         vecs,
         fcns,
         ghost,
         ghost_data,
         mg_data->dgmath_jit_dbase,
         cheby->smoother_eigs_cg_imax,
         &cheby->eigs[level]
        );

      cheby->eigs[level] *= cheby->smoother_eigs_max_multiplier;
  }

  int iter = cheby->smoother_imax;
  double lmin = cheby->eigs[level]/cheby->smoother_eigs_lmax_lmin_ratio;
  double lmax = cheby->eigs[level];

  if (cheby->smoother_print_eigs){
    printf("[MG_CHEBY_SMOOTHER]: Lev %d Max_eig %.25f\n", level, cheby->eigs[level]);
  }
  
  multigrid_smoother_cheby_iterate
    (
     p4est,
     vecs,
     fcns,
     ghost,
     ghost_data,
     r,
     iter,
     lmin,
     lmax,
     cheby->smoother_print_residual_norm
    );
}


multigrid_smoother_t*
multigrid_smoother_cheby_init
(
 p4est_t* p4est,
 int num_of_levels,
 const char* input_file
)
{
  multigrid_smoother_t* smoother = P4EST_ALLOC(multigrid_smoother_t, 1);
  multigrid_smoother_cheby_t* smoother_data = P4EST_ALLOC(multigrid_smoother_cheby_t, 1);
  
  smoother_data->eigs = P4EST_ALLOC(double, num_of_levels);  
  
  /* set externally in input file */
  smoother_data->smoother_imax = -1;
  smoother_data->smoother_eigs_cg_imax = -1;
  smoother_data->smoother_eigs_lmax_lmin_ratio = -1;
  smoother_data->smoother_eigs_max_multiplier = -1;
  smoother_data->smoother_eigs_reuse_fromdownvcycle = -1;
  smoother_data->smoother_eigs_reuse_fromlastvcycle = -1;
  smoother_data->smoother_print_residual_norm = -1;
  smoother_data->smoother_print_eigs = -1;

  /* set internally */
  smoother_data->mpirank = p4est->mpirank;
  smoother_data->smoother_eigs_compute = -1;
  
  if (ini_parse(input_file, multigrid_smoother_cheby_input_handler, smoother_data) < 0) {
    mpi_abort("Can't load input file");
  }
  if(smoother_data->smoother_imax == -1){
    mpi_abort("[D4EST_ERROR]: smoother_imax not set in multigrid input");
  }
  if(smoother_data->smoother_eigs_cg_imax == -1){
    mpi_abort("[D4EST_ERROR]: smoother_eigs_cg_imax not set in multigrid input");
  }
  if(smoother_data->smoother_eigs_lmax_lmin_ratio == -1){
    mpi_abort("[D4EST_ERROR]: smoother_eigs_lmax_lmin_ratio not set in multigrid input");
  }
  if(smoother_data->smoother_eigs_max_multiplier == -1){
    mpi_abort("[D4EST_ERROR]: smoother_eigs_max_multiplier not set in multigrid input");
  }
  if(smoother_data->smoother_eigs_reuse_fromdownvcycle == -1){
    mpi_abort("[D4EST_ERROR]: smoother_eigs_reuse_fromdownvcycle not set in multigrid input");
  }
  if(smoother_data->smoother_eigs_reuse_fromlastvcycle == -1){
    mpi_abort("[D4EST_ERROR]: smoother_eigs_reuse_fromlastvcycle not set in multigrid input");
  }
  if(smoother_data->smoother_print_residual_norm == -1){
    mpi_abort("[D4EST_ERROR]: smoother_print_residual_norm  not set in multigrid input");
  }
  if(smoother_data->smoother_print_eigs == -1){
    mpi_abort("[D4EST_ERROR]: smoother_print_eigs not set in multigrid input");
  }
  
  if(p4est->mpirank == 0){
    printf("[D4EST_INFO]: Multigrid_Smoother_CHEBY Parameters\n");
    printf("[D4EST_INFO]: Smoother imax = %d\n", smoother_data->smoother_imax);
    printf("[D4EST_INFO]: Smoother eigs cg max = %d\n", smoother_data->smoother_eigs_cg_imax);
    printf("[D4EST_INFO]: Smoother eigs lmax_lmin_ratio = %f\n", smoother_data->smoother_eigs_lmax_lmin_ratio);
    printf("[D4EST_INFO]: Smoother eigs max multiplier = %.25f\n", smoother_data->smoother_eigs_max_multiplier);
    printf("[D4EST_INFO]: Smoother eigs reuse up vcycle = %d\n", smoother_data->smoother_eigs_reuse_fromdownvcycle);
    printf("[D4EST_INFO]: Smoother eigs reuse from last vcycle = %d\n", smoother_data->smoother_eigs_reuse_fromlastvcycle);
    printf("[D4EST_INFO]: Smoother print residual norm = %d\n", smoother_data->smoother_print_residual_norm);
    printf("[D4EST_INFO]: Smoother print eigs = %d\n", smoother_data->smoother_print_eigs);
  }

  smoother->user = smoother_data;
  smoother->smooth = multigrid_smoother_cheby;
  smoother->update = multigrid_smoother_cheby_update;

  return smoother;
}
