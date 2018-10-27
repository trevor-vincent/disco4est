#include <d4est_util.h>
#include <pXest.h>
#include <d4est_linalg.h>
#include <d4est_mesh.h>
#include <d4est_elliptic_data.h>
#include <d4est_solver_krylov_petsc.h>
#include <petscsnes.h>
#include <d4est_krylov_pc.h>
#include <d4est_checkpoint.h>
#include <ini.h>
#include <time.h>

static
int d4est_solver_krylov_petsc_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_solver_krylov_petsc_params_t* pconfig = (d4est_solver_krylov_petsc_params_t*)user;
  if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_atol")) {
    D4EST_ASSERT(pconfig->ksp_atol[0] == '*');
    snprintf (pconfig->ksp_atol, sizeof(pconfig->ksp_atol), "%s", value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_rtol")) {
    D4EST_ASSERT(pconfig->ksp_rtol[0] == '*');
    snprintf (pconfig->ksp_rtol, sizeof(pconfig->ksp_rtol), "%s", value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_max_it")) {
    D4EST_ASSERT(pconfig->ksp_max_it[0] == '*');
    snprintf (pconfig->ksp_max_it, sizeof(pconfig->ksp_max_it), "%s", value);
    D4EST_ASSERT(atoi(value) > -1);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_view")) {
    D4EST_ASSERT(pconfig->ksp_view == -1);
    pconfig->ksp_view = atoi(value);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_monitor")) {
    D4EST_ASSERT(pconfig->ksp_monitor == -1);
    pconfig->ksp_monitor = atoi(value);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_converged_reason")) {
    D4EST_ASSERT(pconfig->ksp_converged_reason == -1);
    pconfig->ksp_converged_reason = atoi(value);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_monitor_singular_value")) {
    pconfig->ksp_monitor_singular_value = atoi(value);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_initial_guess_nonzero")) {
    D4EST_ASSERT(pconfig->ksp_initial_guess_nonzero == -1);
    pconfig->ksp_initial_guess_nonzero = atoi(value);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_type")) {
    D4EST_ASSERT(pconfig->ksp_type[0] == '*');
    snprintf (pconfig->ksp_type, sizeof(pconfig->ksp_type), "%s", value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_chebyshev_esteig_steps")) {
    D4EST_ASSERT(pconfig->ksp_chebyshev_esteig_steps[0] == '*');
    snprintf (pconfig->ksp_chebyshev_esteig_steps, sizeof(pconfig->ksp_chebyshev_esteig_steps), "%s", value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_chebyshev_esteig")) {
    D4EST_ASSERT(pconfig->ksp_chebyshev_esteig[0] == '*');
    snprintf (pconfig->ksp_chebyshev_esteig, sizeof(pconfig->ksp_chebyshev_esteig), "%s", value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_chebyshev_esteig_random")) {
    D4EST_ASSERT(pconfig->ksp_chebyshev_esteig_random == -1);
    pconfig->ksp_chebyshev_esteig_random = atoi(value);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_do_not_use_preconditioner")) {
    D4EST_ASSERT(pconfig->ksp_do_not_use_preconditioner == 0);
    pconfig->ksp_do_not_use_preconditioner = atoi(value);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"checkpoint_every_n_krylov_its")) {
    D4EST_ASSERT(pconfig->checkpoint_every_n_krylov_its == 0);
    pconfig->checkpoint_every_n_krylov_its = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void
d4est_solver_krylov_petsc_set_options_database_from_params
(
 d4est_solver_krylov_petsc_params_t* input
)
{
  if(input->ksp_monitor)
    PetscOptionsSetValue(NULL,"-ksp_monitor","");
  else
    PetscOptionsClearValue(NULL,"-ksp_monitor");

  if(input->ksp_view)
    PetscOptionsSetValue(NULL,"-ksp_view","");
  else
    PetscOptionsClearValue(NULL,"-ksp_view");

  if(input->ksp_converged_reason)
     PetscOptionsSetValue(NULL,"-ksp_converged_reason","");
  else
    PetscOptionsClearValue(NULL,"-ksp_converged_reason");

  if(input->ksp_initial_guess_nonzero)
    PetscOptionsSetValue(NULL,"-ksp_initial_guess_nonzero","");
  else
    PetscOptionsClearValue(NULL,"-ksp_initial_guess_nonzero");

  if(input->ksp_monitor_singular_value == 1)
    PetscOptionsSetValue(NULL,"-ksp_monitor_singular_value","");
  else
    PetscOptionsClearValue(NULL,"-ksp_monitor_singular_value");
  
  PetscOptionsClearValue(NULL,"-ksp_type");
  PetscOptionsSetValue(NULL,"-ksp_type",input->ksp_type);
  
  PetscOptionsClearValue(NULL,"-ksp_atol");
  PetscOptionsSetValue(NULL,"-ksp_atol",input->ksp_atol);
  
  PetscOptionsClearValue(NULL,"-ksp_rtol");
  PetscOptionsSetValue(NULL,"-ksp_rtol",input->ksp_rtol);
  
  PetscOptionsClearValue(NULL,"-ksp_max_it");
  PetscOptionsSetValue(NULL,"-ksp_max_it", input->ksp_max_it);

  if(d4est_util_match(input->ksp_type,"chebyshev")){
    PetscOptionsClearValue(NULL,"-ksp_chebyshev_esteig_steps");
    PetscOptionsSetValue(NULL,"-ksp_chebyshev_esteig_steps",input->ksp_chebyshev_esteig_steps);
    PetscOptionsClearValue(NULL,"-ksp_chebyshev_esteig");
    PetscOptionsSetValue(NULL,"-ksp_chebyshev_esteig",input->ksp_chebyshev_esteig);
    if(input->ksp_chebyshev_esteig_random)
      PetscOptionsSetValue(NULL,"-ksp_chebyshev_esteig_random","");
    else
      PetscOptionsClearValue(NULL,"-ksp_chebyshev_esteig_random");
  }
}

void
d4est_solver_krylov_petsc_input
(
 p4est_t* p4est,
 const char* input_file,
 const char* input_section,
 d4est_solver_krylov_petsc_params_t* input
)
{
  input->ksp_view = -1;
  input->ksp_monitor = -1;
  input->ksp_converged_reason = -1;
  input->ksp_initial_guess_nonzero = -1;
  input->ksp_type[0] = '*';
  input->ksp_atol[0] = '*';
  input->ksp_rtol[0] = '*';
  input->ksp_max_it[0] = '*';
  input->ksp_chebyshev_esteig_steps[0] = '*';
  input->ksp_chebyshev_esteig[0] = '*';
  input->ksp_chebyshev_esteig_random = -1;
  input->ksp_monitor_singular_value = 0;
  input->ksp_do_not_use_preconditioner = 0;
  input->checkpoint_every_n_krylov_its = 0;
  
  D4EST_ASSERT(sizeof(input->input_section) <= 50);
  snprintf (input->input_section, sizeof(input->input_section), "%s", input_section);
  if (ini_parse(input_file, d4est_solver_krylov_petsc_input_handler, input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, input->ksp_view, -1);
  D4EST_CHECK_INPUT(input_section, input->ksp_monitor, -1);
  D4EST_CHECK_INPUT(input_section, input->ksp_converged_reason, -1);
  D4EST_CHECK_INPUT(input_section, input->ksp_initial_guess_nonzero, -1);
  D4EST_CHECK_INPUT(input_section, input->ksp_type[0], '*');
  D4EST_CHECK_INPUT(input_section, input->ksp_atol[0], '*');
  D4EST_CHECK_INPUT(input_section, input->ksp_rtol[0], '*');
  D4EST_CHECK_INPUT(input_section, input->ksp_max_it[0], '*');

  if(d4est_util_match(input->ksp_type,"chebyshev")){
    D4EST_CHECK_INPUT(input_section, input->ksp_chebyshev_esteig_steps[0], '*');
    D4EST_CHECK_INPUT(input_section, input->ksp_chebyshev_esteig[0], '*');
    D4EST_CHECK_INPUT(input_section, input->ksp_chebyshev_esteig_random, -1);
  }
    
  if(p4est->mpirank == 0){
    zlog_category_t *c_default = zlog_get_category("d4est_solver_krylov_petsc");
    zlog_debug(c_default, "ksp_type = %s", input->ksp_type);
    zlog_debug(c_default, "ksp_view = %d", input->ksp_view);
    zlog_debug(c_default, "ksp_monitor = %d", input->ksp_monitor);
    zlog_debug(c_default, "ksp_atol = %s", input->ksp_atol);
    zlog_debug(c_default, "ksp_rtol = %s", input->ksp_rtol);
    zlog_debug(c_default, "ksp_maxit = %s", input->ksp_max_it);
    zlog_debug(c_default, "ksp_converged_reason = %d", input->ksp_converged_reason);
    zlog_debug(c_default, "ksp_initial_guess_nonzero = %d", input->ksp_initial_guess_nonzero);
    zlog_debug(c_default, "ksp_do_not_use_preconditioner = %d", input->ksp_do_not_use_preconditioner);
    if(d4est_util_match(input->ksp_type,"chebyshev")){
      zlog_debug(c_default, "ksp_chebyshev_esteig_steps = %s", input->ksp_chebyshev_esteig_steps);
      zlog_debug(c_default, "ksp_chebyshev_esteig = %s", input->ksp_chebyshev_esteig);
      zlog_debug(c_default, "ksp_chebyshev_esteig_random = %d", input->ksp_chebyshev_esteig_random);
    }
  }  
}

static
PetscErrorCode d4est_solver_krylov_petsc_monitor(KSP ksp,PetscInt it, PetscReal norm, void *ctx)
{
  krylov_ctx_t* petsc_ctx = (krylov_ctx_t*) ctx;
  zlog_category_t* c_default = zlog_get_category("d4est_solver_krylov_petsc");
  zlog_category_t* its_output = zlog_get_category("d4est_solver_iteration_info");

  if (petsc_ctx->p4est->mpirank == 0){
    double duration_seconds = ((double)(clock() - petsc_ctx->time_start)) / CLOCKS_PER_SEC;
    zlog_info(its_output, "AMR_IT KSP_IT KSP_NORM TIME: %d %d %.25f %f", petsc_ctx->amr_level, it, norm, duration_seconds);
  }  
  if (petsc_ctx->checkpoint_every_n_krylov_its > 0 && petsc_ctx->amr_level >= 0){
    int it;
    KSPGetIterationNumber(*(petsc_ctx->ksp),&it);
    if ( (it - petsc_ctx->last_krylov_checkpoint_it) >= petsc_ctx->checkpoint_every_n_krylov_its ){
      char* output = NULL;
      asprintf(&output,"checkpoint_krylov_%d", it);

      d4est_elliptic_data_t* vecs = petsc_ctx->vecs;
      d4est_checkpoint_save
        (
         petsc_ctx->amr_level,
         output,
         petsc_ctx->p4est,
         NULL,
         NULL,
         (const char * []){"u", NULL},
         (hid_t []){H5T_NATIVE_DOUBLE},
         (int []){vecs->local_nodes},
         (void* []){vecs->u}
        );

      petsc_ctx->last_krylov_checkpoint_it = it;
      free(output);
    }
    else {
      if (petsc_ctx->p4est->mpirank == 0){
        zlog_info(c_default, "No checkpoint at krylov iteration %d", it);
        zlog_info(c_default, "petsc_ctx->last_krylov_checkpoint_it = %d", petsc_ctx->last_krylov_checkpoint_it);
        zlog_info(c_default,"petsc_ctx->checkpoint_every_n_krylov_its = %d",petsc_ctx->checkpoint_every_n_krylov_its);
      }
    }
  }

  return 0;
}

static
PetscErrorCode d4est_solver_krylov_petsc_apply_aij( Mat A, Vec x, Vec y )
{
  void           *ctx;
  PetscErrorCode ierr;

  krylov_ctx_t* petsc_ctx;
  const double* px;
  double* py;

  /* PetscFunctionBegin; */
  ierr = MatShellGetContext( A, &ctx ); CHKERRQ(ierr);
  petsc_ctx = (krylov_ctx_t *)ctx;
  ierr = VecGetArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecGetArray( y, &py ); CHKERRQ(ierr);

  /* d4est_elliptic_eqns_t* fcns = petsc_ctx->fcns; */
  p4est_t* p4est = petsc_ctx->p4est;
  d4est_ghost_t* ghost = *petsc_ctx->ghost;

  d4est_elliptic_data_t vecs_for_aij;
  d4est_elliptic_data_copy_ptrs(petsc_ctx->vecs, &vecs_for_aij);

  vecs_for_aij.u = (double*)px;
  vecs_for_aij.Au = py;

  d4est_elliptic_eqns_apply_lhs
    (
     petsc_ctx->p4est,
     *petsc_ctx->ghost,
     *petsc_ctx->ghost_data,
     petsc_ctx->fcns,
     &vecs_for_aij,
     petsc_ctx->d4est_ops,
     petsc_ctx->d4est_geom,
     petsc_ctx->d4est_quad,
     petsc_ctx->d4est_factors
    );

  ierr = VecRestoreArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecRestoreArray( y, &py ); CHKERRQ(ierr);
  return ierr;
}

d4est_solver_krylov_petsc_info_t
d4est_solver_krylov_petsc_solve
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 d4est_ghost_t** ghost,
 d4est_ghost_data_t** ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_krylov_petsc_params_t* d4est_solver_krylov_petsc_params,
 d4est_krylov_pc_t* d4est_krylov_pc,
 int amr_level
)
{
  zlog_category_t *c_default = zlog_get_category("d4est_solver_krylov_petsc");
  clock_t start = clock();
  if (p4est->mpirank == 0) {
    zlog_info(c_default, "Performing Krylov PETSc solve...");
  }

  d4est_solver_krylov_petsc_set_options_database_from_params(d4est_solver_krylov_petsc_params);

  d4est_solver_krylov_petsc_info_t info;
  KSP ksp;
  Vec x,b;
  PC             pc;
  /* double* u_temp; */
  /* double* rhs_temp; */

  krylov_ctx_t petsc_ctx;
  petsc_ctx.p4est = p4est;
  petsc_ctx.vecs = vecs;
  petsc_ctx.fcns = fcns;
  petsc_ctx.ghost = ghost;
  petsc_ctx.ghost_data = ghost_data;
  petsc_ctx.d4est_ops = d4est_ops;
  petsc_ctx.d4est_geom = d4est_geom;
  petsc_ctx.d4est_quad = d4est_quad;
  petsc_ctx.d4est_factors = d4est_factors;
  petsc_ctx.checkpoint_every_n_krylov_its = d4est_solver_krylov_petsc_params->checkpoint_every_n_krylov_its;
  petsc_ctx.last_krylov_checkpoint_it = 0;
  petsc_ctx.amr_level = amr_level;
  petsc_ctx.time_start = start;

  petsc_ctx.ksp = &ksp;

  
  int local_nodes = vecs->local_nodes;
  double* u = vecs->u;
  double* rhs = vecs->rhs;

  /* printf("PRE KRYLOV SOLVE REDUCTIONS\n"); */
  /* DEBUG_PRINT_ARR_DBL_SUM(u, local_nodes); */
  /* DEBUG_PRINT_ARR_DBL_SUM(rhs, local_nodes); */
  
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  VecCreate(PETSC_COMM_WORLD,&x);//CHKERRQ(ierr);
  VecSetSizes(x, local_nodes, PETSC_DECIDE);//CHKERRQ(ierr);
  VecSetFromOptions(x);//CHKERRQ(ierr);
  VecDuplicate(x,&b);//CHKERRQ(ierr);

  
  KSPGetPC(ksp,&pc);
  if (d4est_krylov_pc != NULL && d4est_solver_krylov_petsc_params->ksp_do_not_use_preconditioner == 0) {
    PCSetType(pc,PCSHELL);//CHKERRQ(ierr);
    d4est_krylov_pc->pc_ctx = &petsc_ctx;
    PCShellSetApply(pc, d4est_solver_krylov_petsc_pc_apply);//CHKERRQ(ierr);
    if (d4est_krylov_pc->pc_setup != NULL){
      PCShellSetSetUp(pc, d4est_solver_krylov_petsc_pc_setup);
    }
    PCShellSetContext(pc, d4est_krylov_pc);//CHKERRQ(ierr);
  }
  else {
    PCSetType(pc,PCNONE);//CHKERRQ(ierr);
  }

  KSPSetResidualHistory(ksp,
                        PETSC_NULL,   // pointer to the array which holds the history
                        PETSC_DECIDE, // size of the array holding the history
                        PETSC_TRUE);  // Whether or not to reset the history for each solve.
  
  KSPSetFromOptions(ksp);
  /* KSPSetUp(ksp); */
  /* Create matrix-free shell for Aij */
  Mat A;
  MatCreateShell
    (
     PETSC_COMM_WORLD,
     local_nodes,
     local_nodes,
     PETSC_DETERMINE,
     PETSC_DETERMINE,
     (void*)&petsc_ctx,
     &A
    );
  MatShellSetOperation(A,MATOP_MULT,(void(*)())d4est_solver_krylov_petsc_apply_aij);

  /* Set Amat and Pmat, where Pmat is the matrix the Preconditioner needs */
  KSPSetOperators(ksp,A,A);
  VecPlaceArray(b, rhs);
  VecPlaceArray(x, u);
  KSPMonitorSet(ksp, d4est_solver_krylov_petsc_monitor, &petsc_ctx, NULL);

  
  KSPSolve(ksp,b,x);
  
  KSPGetIterationNumber(ksp, &(info.total_krylov_iterations));
  KSPGetResidualNorm(ksp, &(info.residual_norm));

  MatDestroy(&A);
  VecResetArray(b);
  VecResetArray(x);
  VecDestroy(&x);
  VecDestroy(&b);
  KSPDestroy(&ksp);
  
  if (p4est->mpirank == 0) {
    double duration_seconds = ((double)(clock() - start)) / CLOCKS_PER_SEC;
    zlog_info(c_default, "Krylov PETSc solve complete in %.2f seconds (%d iterations). Residual norm: %.2e", duration_seconds, info.total_krylov_iterations, info.residual_norm);
  }

  return info;
}
