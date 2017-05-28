#include "../pXest/pXest.h"

#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../Solver/newton_petsc.h"
#include "../Solver/krylov_petsc.h"
#include "petscsnes.h"
#include <ini.h>

static
int newton_petsc_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  newton_petsc_params_t* pconfig = (newton_petsc_params_t*)user;

  if (util_match_couple(section,"newton_petsc",name,"snes_type")) {
    mpi_assert(pconfig->snes_type[0] == '*');
    snprintf (pconfig->snes_type, sizeof(pconfig->snes_type), "%s", value);
  }
  else if (util_match_couple(section,"newton_petsc",name,"snes_atol")) {
    mpi_assert(pconfig->snes_atol[0] == '*');
    snprintf (pconfig->snes_atol, sizeof(pconfig->snes_atol), "%s", value);
  }
  else if (util_match_couple(section,"newton_petsc",name,"snes_rtol")) {
    mpi_assert(pconfig->snes_rtol[0] == '*');
    snprintf (pconfig->snes_rtol, sizeof(pconfig->snes_rtol), "%s", value);
  }
  else if (util_match_couple(section,"newton_petsc",name,"snes_stol")) {
    mpi_assert(pconfig->snes_stol[0] == '*');
    snprintf (pconfig->snes_stol, sizeof(pconfig->snes_stol), "%s", value);
  }
 else if (util_match_couple(section,"newton_petsc",name,"snes_trtol")) {
    mpi_assert(pconfig->snes_trtol[0] == '*');
    snprintf (pconfig->snes_trtol, sizeof(pconfig->snes_trtol), "%s", value);
  } 
  else if (util_match_couple(section,"newton_petsc",name,"snes_max_funcs")) {
    mpi_assert(pconfig->snes_max_funcs[0] == '*');
    snprintf (pconfig->snes_max_funcs, sizeof(pconfig->snes_max_funcs), "%s", value);
  }
  else if (util_match_couple(section,"newton_petsc",name,"snes_max_it")) {
    mpi_assert(pconfig->snes_max_it[0] == '*');
    snprintf (pconfig->snes_max_it, sizeof(pconfig->snes_max_it), "%s", value);
  }
  else if (util_match_couple(section,"newton_petsc",name,"snes_linesearch_monitor")) {
    mpi_assert(pconfig->snes_linesearch_monitor == -1);
    pconfig->snes_linesearch_monitor = atoi(value);
    mpi_assert(atoi(value) == 0 || atoi(value) == 1);
  }  
  else if (util_match_couple(section,"newton_petsc",name,"snes_linesearch_order")) {
    mpi_assert(pconfig->snes_linesearch_order[0] == '*');
    snprintf (pconfig->snes_linesearch_order, sizeof(pconfig->snes_linesearch_order), "%s", value);
  }
  else if (util_match_couple(section,"newton_petsc",name,"snes_monitor")) {
    mpi_assert(pconfig->snes_monitor == -1);
    pconfig->snes_monitor = atoi(value);
    mpi_assert(atoi(value) == 0 || atoi(value) == 1);
  }
  else if (util_match_couple(section,"newton_petsc",name,"snes_view")) {
    mpi_assert(pconfig->snes_view == -1);
    pconfig->snes_view = atoi(value);
    mpi_assert(atoi(value) == 0 || atoi(value) == 1);
  }    
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


void
newton_petsc_input
(
 p4est_t* p4est,
 const char* input_file,
 const char* printf_prefix,
 newton_petsc_params_t* input
)
{
  input->snes_atol[0] = '*';
  input->snes_rtol[0] = '*';
  input->snes_stol[0] = '*';
  input->snes_max_funcs[0] = '*';
  input->snes_type[0] = '*';
  input->snes_max_it[0] = '*';
  input->snes_trtol[0] = '*';
  input->snes_linesearch_order[0] = '*';
  input->snes_monitor = -1;
  input->snes_view = -1;
  input->snes_linesearch_monitor = -1;
  input->snes_converged_reason = -1;


  if (ini_parse(input_file, newton_petsc_input_handler, input) < 0) {
    mpi_abort("Can't load input file");
  }  

  D4EST_CHECK_INPUT("newton_petsc", input->snes_type[0], '*');
  D4EST_CHECK_INPUT("newton_petsc", input->snes_atol[0], '*');
  D4EST_CHECK_INPUT("newton_petsc", input->snes_rtol[0], '*');
  D4EST_CHECK_INPUT("newton_petsc", input->snes_max_it[0], '*');
  D4EST_CHECK_INPUT("newton_petsc", input->snes_max_funcs[0], '*');

  if(util_match(input->snes_type,"newtonls")){
    D4EST_CHECK_INPUT("newton_petsc", input->snes_linesearch_order[0], '*');
    D4EST_CHECK_INPUT("newton_petsc", input->snes_linesearch_monitor, -1);
  }
  else if(util_match(input->snes_type,"newtontr")){
    D4EST_CHECK_INPUT("newton_petsc", input->snes_trtol[0], '*');
  }  
  
  D4EST_CHECK_INPUT("newton_petsc", input->snes_view, -1);
  D4EST_CHECK_INPUT("newton_petsc", input->snes_monitor, -1);

  if(p4est->mpirank == 0){
    printf("%s: snes_type = %s\n",printf_prefix, input->snes_type);
    printf("%s: snes_view = %d\n",printf_prefix, input->snes_view);
    printf("%s: snes_monitor = %d\n",printf_prefix, input->snes_monitor);
    printf("%s: snes_atol = %s\n",printf_prefix, input->snes_atol);
    printf("%s: snes_rtol = %s\n",printf_prefix, input->snes_rtol);
    printf("%s: snes_stol = %s\n",printf_prefix, input->snes_stol);
    printf("%s: snes_maxit = %s\n",printf_prefix, input->snes_max_it);
    printf("%s: snes_maxfuncs = %s\n",printf_prefix, input->snes_max_funcs);
    if(util_match(input->snes_type,"newtonls")){
      printf("%s: snes_linesearch_order = %s\n",printf_prefix, input->snes_max_funcs);
      printf("%s: snes_linesearch_monitor = %d\n",printf_prefix, input->snes_linesearch_monitor);
    }
    if(util_match(input->snes_type,"newtontr")){
      printf("%s: snes_trtol = %s\n",printf_prefix, input->snes_trtol);
    }
  } 
}


/** 
 * Required for grabbing x0 
 * in J(x0)dx = -F(x0)
 * 
 * @param MojSnes 
 * @param X 
 * @param Jac 
 * @param Precond 
 * @param flags 
 * @param ctx 
 * 
 * @return 
 */
static
PetscErrorCode newton_petsc_save_x0
(
 SNES MojSnes,
 Vec x0,
 Mat Jac, 
 Mat Precond,
 void* ctx
) 
{

  PetscErrorCode ierr;
  petsc_ctx_t* petsc_ctx = (petsc_ctx_t*) ctx;
  problem_data_t* vecs = petsc_ctx->vecs;
  const double* px0;
  ierr = VecGetArrayRead( x0, &px0 ); CHKERRQ(ierr);

  int i;
  for (i = 0; i < vecs->local_nodes; i++){
    vecs->u0[i] = px0[i];
    /* printf("u0[%d] = %f\n",i, vecs->u0[i]); */
  }
    
  VecRestoreArrayRead
    (
     x0,
     &px0
    );
  
  return 0; 
}


static
PetscErrorCode newton_petsc_get_residual(SNES snes, Vec x, Vec f, void *ctx){

  const double* xx;
  petsc_ctx_t* petsc_ctx = (petsc_ctx_t*) ctx;

  double* ftemp;

  VecGetArray(f,&ftemp); 
  VecGetArrayRead(x,&xx);

  problem_data_t* vecs = petsc_ctx->vecs;
  p4est_t* p4est = petsc_ctx->p4est;
  p4est_ghost_t* ghost = *petsc_ctx->ghost;
  d4est_operators_t* d4est_ops = petsc_ctx->d4est_ops;
  d4est_geometry_t* d4est_geom = petsc_ctx->d4est_geom;
  d4est_quadrature_t* d4est_quad = petsc_ctx->d4est_quad;
  
  problem_data_t vecs_for_res_build;
  problem_data_copy_ptrs(vecs, &vecs_for_res_build);
  vecs_for_res_build.u = (double*)xx;//x_temp;
  vecs_for_res_build.u0 = (double*)xx;//x_temp;
  vecs_for_res_build.Au = ftemp;

  /* if (d4est_geom != NULL){ */
  /*   ((curved_weakeqn_ptrs_t*)(petsc_ctx->fcns))->build_residual(p4est, ghost, (d4est_element_data_t*)(petsc_ctx->ghost_data), &vecs_for_res_build, d4est_ops,d4est_geom); */
  /* } */
  /* else { */
  ((weakeqn_ptrs_t*)(petsc_ctx->fcns))->build_residual(p4est, ghost, (petsc_ctx->ghost_data), &vecs_for_res_build, d4est_ops, d4est_geom, d4est_quad);
  /* } */
  VecRestoreArray(f,&ftemp);//CHKERRQ(ierr);
  VecRestoreArrayRead(x,&xx);

  /* P4EST_FREE(x_temp); */
  
  return 0;
}

static
PetscErrorCode newton_petsc_apply_jacobian( Mat jac, Vec x, Vec y )
{
  void           *ctx;
  PetscErrorCode ierr;

  petsc_ctx_t* petsc_ctx;
  const double* px;
  double* py;

  /* PetscFunctionBegin; */
  ierr = MatShellGetContext( jac, &ctx ); CHKERRQ(ierr);  
  petsc_ctx = (petsc_ctx_t *)ctx;
  ierr = VecGetArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecGetArray( y, &py ); CHKERRQ(ierr);

  problem_data_t* vecs = petsc_ctx->vecs;
  p4est_t* p4est = petsc_ctx->p4est;
  p4est_ghost_t* ghost = *petsc_ctx->ghost;
  d4est_operators_t* d4est_ops = petsc_ctx->d4est_ops;
  d4est_geometry_t* d4est_geom = petsc_ctx->d4est_geom;
  d4est_quadrature_t* d4est_quad = petsc_ctx->d4est_quad;


  /* int i; */
  /* for (i = 0; i < vecs->local_nodes; i++){ */
    /* vecs->u[i] = px[i]; */
  /* } */

  problem_data_t vecs_for_jac;
  problem_data_copy_ptrs(vecs, &vecs_for_jac);

  vecs_for_jac.u = (double*)px;
  vecs_for_jac.Au = py;
  /* vecs->u is already set above */
  /* vecs->u0 is already set by newton_petsc_save_x0 */


  /* if (d4est_geom != NULL){ */
    /* ((curved_weakeqn_ptrs_t*)(petsc_ctx->fcns))->apply_lhs(p4est, ghost, (d4est_element_data_t*)(petsc_ctx->ghost_data), &vecs_for_jac, d4est_ops,d4est_geom); */
  /* } */
  /* else { */
((weakeqn_ptrs_t*)(petsc_ctx->fcns))->apply_lhs(p4est,
                                                ghost,
                                                (petsc_ctx->ghost_data),
                                                &vecs_for_jac,
                                                d4est_ops,
                                                d4est_geom,
                                                d4est_quad
                                               );
  /* } */
  
  ierr = VecRestoreArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecRestoreArray( y, &py ); CHKERRQ(ierr);
  return ierr;
}

void
newton_petsc_set_options_database_from_params
(
 newton_petsc_params_t* input
)
{
  if(input->snes_monitor)
    PetscOptionsSetValue(NULL,"-snes_monitor","");
  else
    PetscOptionsClearValue(NULL,"-snes_monitor");

  if(input->snes_view)
    PetscOptionsSetValue(NULL,"-snes_view","");
  else
    PetscOptionsClearValue(NULL,"-snes_view");

  if(input->snes_converged_reason)
     PetscOptionsSetValue(NULL,"-snes_converged_reason","");
  else
    PetscOptionsClearValue(NULL,"-snes_converged_reason");

  PetscOptionsClearValue(NULL,"-snes_type");
  PetscOptionsSetValue(NULL,"-snes_type",input->snes_type);
  
  PetscOptionsClearValue(NULL,"-snes_atol");
  PetscOptionsSetValue(NULL,"-snes_atol",input->snes_atol);
  
  PetscOptionsClearValue(NULL,"-snes_rtol");
  PetscOptionsSetValue(NULL,"-snes_rtol",input->snes_rtol);
  
  PetscOptionsClearValue(NULL,"-snes_max_it");
  PetscOptionsSetValue(NULL,"-snes_max_it", input->snes_max_it);

  PetscOptionsClearValue(NULL,"-snes_max_funcs");
  PetscOptionsSetValue(NULL,"-snes_max_funcs", input->snes_max_funcs);

  if(input->snes_stol[0] != '*'){
    PetscOptionsClearValue(NULL,"-snes_stol");
    PetscOptionsSetValue(NULL,"-snes_stol", input->snes_stol);
  }
    
  if(util_match(input->snes_type,"newtonls")){
    PetscOptionsClearValue(NULL,"-snes_linesearch_order");
    PetscOptionsSetValue(NULL,"-snes_linesearch_order", input->snes_linesearch_order);
    if(input->snes_linesearch_monitor){
      PetscOptionsSetValue(NULL,"-snes_linesearch_monitor", "");
    }
    else {
      PetscOptionsClearValue(NULL,"-snes_linesearch_monitor");
    }
  }
  else if(util_match(input->snes_type,"newtontr")){
    PetscOptionsClearValue(NULL,"-snes_trtol");
    PetscOptionsSetValue(NULL,"-snes_trtol", input->snes_trtol);
  }  
  
}

void newton_petsc_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 void* fcns,
 p4est_ghost_t** ghost,
 void** ghost_data, 
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 krylov_petsc_params_t* krylov_options,
 newton_petsc_params_t* newton_options,
 krylov_pc_t* krylov_pc
)
{
  krylov_petsc_set_options_database_from_params(krylov_options);
  newton_petsc_set_options_database_from_params(newton_options);
  
  SNES snes;
  KSP ksp;
  Vec x,r;
  /* PetscMPIInt size; */
  /* double* u_temp; */
  SNESCreate(PETSC_COMM_WORLD,&snes);//CHKERRQ(ierr);
  VecCreate(PETSC_COMM_WORLD,&x);//CHKERRQ(ierr);
  VecSetSizes(x, vecs->local_nodes, PETSC_DECIDE);//CHKERRQ(ierr);
  VecSetFromOptions(x);//CHKERRQ(ierr);
  VecDuplicate(x,&r);//CHKERRQ(ierr);
  
  petsc_ctx_t petsc_ctx;
  petsc_ctx.p4est = p4est;
  petsc_ctx.vecs = vecs;
  petsc_ctx.fcns = fcns;
  petsc_ctx.ghost = ghost;
  petsc_ctx.ghost_data = ghost_data;
  petsc_ctx.d4est_ops = d4est_ops;
  petsc_ctx.d4est_geom = d4est_geom;
  
  SNESSetFunction(snes,r,newton_petsc_get_residual,(void*)&petsc_ctx);//CHKERRQ(ierr);
  SNESGetKSP(snes,&ksp);
  /* PetscOptionsSetValue(NULL,"-snes_mf",""); */
  /* PetscOptionsSetValue(NULL,"-snes_converged_reason",""); */
  /* PetscOptionsSetValue(NULL,"-ksp_converged_reason",""); */
  /* PetscOptionsSetValue(NULL,"-ksp_monitor_true_residual",""); */
  /* PetscOptionsSetValue(NULL,"-snes_monitor",""); */
  /* PetscOptionsSetValue("-snes_monitor_solution",""); */

  
  /* PetscOptionsSetValue(NULL,"-ksp_atol","1e-50"); */
  /* PetscOptionsSetValue(NULL,"-ksp_rtol","1e-5"); */
  /* PetscOptionsSetValue(NULL,"-ksp_max_it","1000000"); */

  /* PetscOptionsSetValue(NULL,"-snes_rtol","1e-5"); */
  /* PetscOptionsSetValue(NULL,"-snes_atol","1e-50"); */
  /* PetscOptionsSetValue(NULL,"-snes_linesearch_monitor",""); */
  /* PetscOptionsSetValue(NULL,"-ksp_view",""); */
  /* PetscOptionsSetValue(NULL,"-info","100000"); */

  /* absurdly large number */
  /* PetscOptionsSetValue(NULL,"-snes_max_funcs","10000000000000000"); */
  /* PetscOptionsSetValue(NULL,"-ksp_type",krylov_params->ksp_type); */

  /* PetscOptionsSetValue(NULL,"-snes_linesearch_order","3"); */
  /* if (krylov_params->ksp_monitor == 1) */
    /* PetscOptionsSetValue(NULL,"-ksp_monitor",""); */


  /* PetscOptionsSetValue(NULL,"-ksp_monitor",""); */
  /* PetscOptionsSetValue(NULL,"-ksp_converged_reason",""); */
  /* PetscOptionsSetValue(NULL,"-ksp_atol","1e-20"); */
  /* PetscOptionsSetValue(NULL,"-with-debugging","1"); */
  /* PetscOptionsSetValue(NULL,"-ksp_rtol","1e-5"); */
  /* PetscOptionsSetValue(NULL,"-ksp_max_it","1000000"); */
  PC pc;
  KSPGetPC(ksp,&pc);
  if (krylov_pc != NULL) {
    PCSetType(pc,PCSHELL);//CHKERRQ(ierr);
    krylov_pc->pc_ctx = &petsc_ctx;
    PCShellSetApply(pc, krylov_petsc_pc_apply);//CHKERRQ(ierr);
    if (krylov_pc->pc_setup != NULL){
      PCShellSetSetUp(pc, krylov_petsc_pc_setup);
    }
    PCShellSetContext(pc, krylov_pc);//CHKERRQ(ierr);
  }
  else {
    PCSetType(pc,PCNONE);//CHKERRQ(ierr);
  }
  
  /* KSPSetType(ksp, krylov_params.ksp_type); */
  KSPSetFromOptions(ksp);
  
  /* END SET PRECONDITIONER */
  /* END SET PRECONDITIONER */
  /* END SET PRECONDITIONER */

  /* PetscOptionsSetValue(NULL,"-pc_type","none"); */
  /* PetscOptionsSetValue(NULL,"-ksp_type","cg");   */
  /* PetscOptionsSetValue(NULL,"-ksp_view",""); */
  /* PetscOptionsSetValue(NULL,"-snes_mf_operator",""); */
  /* PetscOptionsSetValue(NULL,"-snes_view",""); */
  SNESSetFromOptions(snes);//CHKERRQ(ierr);
    
  double* u0 = P4EST_ALLOC(double, vecs->local_nodes);
  vecs->u0 = u0;
  
  Mat J;
  MatCreateShell
    (
     PETSC_COMM_WORLD,
     vecs->local_nodes,
     vecs->local_nodes,
     PETSC_DETERMINE,
     PETSC_DETERMINE,
     (void*)&petsc_ctx,
     &J
    );
  

  /* MatSetFromOptions(J); */
  MatShellSetOperation(J,MATOP_MULT,(void(*)())newton_petsc_apply_jacobian);
  /* MatCreateSNESMF(snes, &J); */
  SNESSetJacobian(snes,J,J,newton_petsc_save_x0,(void*)&petsc_ctx);

  /* linalg_copy_1st_to_2nd(vecs->u, u_temp, vecs->local_nodes); */

  VecPlaceArray(x, vecs->u);
 
  
  SNESSolve(snes,NULL,x);

  /* linalg_copy_1st_to_2nd(u_temp, vecs->u, vecs->local_nodes); */

  P4EST_FREE(u0);
  /* VecRestoreArray(x,&u_temp); */
  VecResetArray(x);
  VecDestroy(&x);//CHKERRQ(ierr);  
  VecDestroy(&r);//CHKERRQ(ierr);
  SNESDestroy(&snes);//CHKERRQ(ierr);
}
