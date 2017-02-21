#include "../pXest/pXest.h"

#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../Solver/newton_petsc.h"
#include "petscsnes.h"
#include <ini.h>
#include <krylov_petsc_pc.h>


typedef struct {

  int snes_monitor;
  int snes_linesearch_monitor;
  int snes_linesearch_order;
  double snes_atol;
  double snes_rtol;
  int snes_max_funcs;

  int count;
  
} newton_petsc_params_t;


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
  
  if (util_match_couple(section,"solver",name,"snes_atol")) {
    mpi_assert(pconfig->snes_atol == -1);
    pconfig->snes_atol = atof(value);
    PetscOptionsSetValue(NULL,"-snes_atol",value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"solver",name,"snes_rtol")) {
    mpi_assert(pconfig->snes_rtol == -1);
    pconfig->snes_rtol = atof(value);
    PetscOptionsSetValue(NULL,"-snes_rtol",value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"solver",name,"snes_max_funcs")) {
    mpi_assert(pconfig->snes_max_funcs == -1);
    pconfig->snes_max_funcs = atoi(value);
    PetscOptionsSetValue(NULL,"-snes_max_funcs",value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"solver",name,"snes_linesearch_monitor")) {
    mpi_assert(pconfig->snes_linesearch_monitor == -1);
    pconfig->snes_linesearch_monitor = atoi(value);
    mpi_assert(atoi(value) == 0 || atoi(value) == 1);
    if (atoi(value) == 1){
      PetscOptionsSetValue(NULL,"-snes_linesearch_monitor","");
    }
    pconfig->count += 1;
  }  
  else if (util_match_couple(section,"solver",name,"snes_linesearch_order")) {
    mpi_assert(pconfig->snes_linesearch_order == -1);
    pconfig->snes_linesearch_order = atoi(value);
    PetscOptionsSetValue(NULL,"-snes_linesearch_order",value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"solver",name,"snes_monitor")) {
    mpi_assert(pconfig->snes_monitor == -1);
    pconfig->snes_monitor = atoi(value);
    mpi_assert(atoi(value) == 0 || atoi(value) == 1);
    if (atoi(value) == 1){
      PetscOptionsSetValue(NULL,"-snes_monitor","");
    }
    pconfig->count += 1;
  }    
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


newton_petsc_params_t
newton_petsc_input
(
 const char* input_file
)
{
  int num_of_options = 6;
  
  newton_petsc_params_t input;
  input.count = 0;
  input.snes_atol = -1;
  input.snes_rtol = -1;
  input.snes_max_funcs = -1;
  input.snes_monitor = -1;
  input.snes_linesearch_order = -1;
  input.snes_linesearch_monitor = -1;
  
  if (ini_parse(input_file, newton_petsc_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  /* printf("[D4EST_INFO]: snes_atol = %.25f\n", input.snes_atol); */
  /* printf("[D4EST_INFO]: snes_rtol = %.25f\n", input.snes_rtol); */
  /* printf("[D4EST_INFO]: snes_max_funcs = %d\n", input.snes_max_funcs); */
  /* printf("[D4EST_INFO]: snes_monitor = %d\n", input.snes_monitor); */
  /* printf("[D4EST_INFO]: snes_linesearch_order = %d\n", input.snes_linesearch_order); */
  /* printf("[D4EST_INFO]: snes_linesearch_monitor = %d\n", input.snes_linesearch_monitor); */
  
  if (input.count != num_of_options){
    mpi_abort("[D4EST_ERROR]: input.count != num_of_options in newton_petsc_params");
  }
  return input;
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
  newton_petsc_ctx_t* rct = (newton_petsc_ctx_t*) ctx;
  problem_data_t* vecs = rct->vecs;
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
  newton_petsc_ctx_t* rct = (newton_petsc_ctx_t*) ctx;

  double* ftemp;

  VecGetArray(f,&ftemp); 
  VecGetArrayRead(x,&xx);

  problem_data_t* vecs = rct->vecs;
  p4est_t* p4est = rct->p4est;
  p4est_ghost_t* ghost = rct->ghost;
  dgmath_jit_dbase_t* dgmath_jit_dbase = rct->dgmath_jit_dbase;
  d4est_geometry_t* d4est_geom = rct->d4est_geom;
  
  problem_data_t vecs_for_res_build;
  problem_data_copy_ptrs(vecs, &vecs_for_res_build);
  vecs_for_res_build.u = (double*)xx;//x_temp;
  vecs_for_res_build.u0 = (double*)xx;//x_temp;
  vecs_for_res_build.Au = ftemp;

  if (d4est_geom != NULL){
    ((curved_weakeqn_ptrs_t*)(rct->fcns))->build_residual(p4est, ghost, (curved_element_data_t*)(rct->ghost_data), &vecs_for_res_build, dgmath_jit_dbase,d4est_geom);
  }
  else {
    ((weakeqn_ptrs_t*)(rct->fcns))->build_residual(p4est, ghost, (element_data_t*)(rct->ghost_data), &vecs_for_res_build, dgmath_jit_dbase);
  }
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

  newton_petsc_ctx_t* rct;
  const double* px;
  double* py;

  /* PetscFunctionBegin; */
  ierr = MatShellGetContext( jac, &ctx ); CHKERRQ(ierr);  
  rct = (newton_petsc_ctx_t *)ctx;
  ierr = VecGetArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecGetArray( y, &py ); CHKERRQ(ierr);

  problem_data_t* vecs = rct->vecs;
  p4est_t* p4est = rct->p4est;
  p4est_ghost_t* ghost = rct->ghost;
  dgmath_jit_dbase_t* dgmath_jit_dbase = rct->dgmath_jit_dbase;
  d4est_geometry_t* d4est_geom = rct->d4est_geom;


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


  if (d4est_geom != NULL){
    ((curved_weakeqn_ptrs_t*)(rct->fcns))->apply_lhs(p4est, ghost, (curved_element_data_t*)(rct->ghost_data), &vecs_for_jac, dgmath_jit_dbase,d4est_geom);
  }
  else {
    ((weakeqn_ptrs_t*)(rct->fcns))->apply_lhs(p4est, ghost, (element_data_t*)(rct->ghost_data), &vecs_for_jac, dgmath_jit_dbase);
  }
  
  ierr = VecRestoreArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecRestoreArray( y, &py ); CHKERRQ(ierr);
  return ierr;
}



void newton_petsc_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 void* fcns,
 p4est_ghost_t** ghost,
 void** ghost_data, 
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geom,
 const char* input_file,
 krylov_pc_create_fcn_t pc_create,
 krylov_pc_destroy_fcn_t pc_destroy,
 void* pc_data
)
{
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
  /* VecGetArray(x,&u_temp); */
  krylov_petsc_params_t krylov_params = krylov_petsc_input(input_file);
  krylov_params.pc_create = (pc_create == NULL) ? NULL : pc_create;
  krylov_params.pc_destroy = (pc_destroy == NULL) ? NULL : pc_destroy;
  krylov_params.pc_data = pc_data;
  newton_petsc_params_t newton_params = newton_petsc_input(input_file);
  
  newton_petsc_ctx_t rct;
  rct.p4est = p4est;
  rct.vecs = vecs;
  rct.fcns = fcns;
  rct.ghost = *ghost;
  rct.ghost_data = *ghost_data;
  rct.dgmath_jit_dbase = dgmath_jit_dbase;
  rct.d4est_geom = d4est_geom;
  SNESSetFunction(snes,r,newton_petsc_get_residual,(void*)&rct);//CHKERRQ(ierr);
  SNESGetKSP(snes,&ksp);
  /* PetscOptionsSetValue(NULL,"-snes_mf",""); */
  PetscOptionsSetValue(NULL,"-snes_converged_reason","");
  PetscOptionsSetValue(NULL,"-ksp_converged_reason","");
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


  /* SET PRECONDITIONER */
  /* SET PRECONDITIONER */
  /* SET PRECONDITIONER */
  krylov_pc_ctx_t kct;
  kct.p4est = p4est;
  kct.vecs = vecs;
  kct.fcns = fcns;
  kct.ghost = ghost;
  kct.ghost_data = ghost_data;
  kct.dgmath_jit_dbase = dgmath_jit_dbase;
  kct.d4est_geom = d4est_geom;
  /* PetscOptionsSetValue(NULL,"-ksp_monitor",""); */
  /* PetscOptionsSetValue(NULL,"-ksp_converged_reason",""); */
  /* PetscOptionsSetValue(NULL,"-ksp_atol","1e-20"); */
  /* PetscOptionsSetValue(NULL,"-with-debugging","1"); */
  /* PetscOptionsSetValue(NULL,"-ksp_rtol","1e-5"); */
  /* PetscOptionsSetValue(NULL,"-ksp_max_it","1000000"); */
  PC pc;
  KSPGetPC(ksp,&pc);
  
  krylov_pc_t* kp = NULL;
  if (krylov_params.ksp_user_defined_pc) {
    PCSetType(pc,PCSHELL);//CHKERRQ(ierr);
    kp = krylov_params.pc_create(&kct);
    PCShellSetApply(pc, krylov_petsc_pc_apply);//CHKERRQ(ierr);
    PCShellSetSetUp(pc, krylov_petsc_pc_setup);
    PCShellSetContext(pc, kp);//CHKERRQ(ierr);
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
     (void*)&rct,
     &J
    );
  

  /* MatSetFromOptions(J); */
  MatShellSetOperation(J,MATOP_MULT,(void(*)())newton_petsc_apply_jacobian);
  /* MatCreateSNESMF(snes, &J); */
  SNESSetJacobian(snes,J,J,newton_petsc_save_x0,(void*)&rct);

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
