#include "../pXest/pXest.h"

#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../Solver/newton_petsc.h"
#include "petscsnes.h"
/* #include <krylov_petsc.h> */
#include <krylov_petsc_pc.h>

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
  weakeqn_ptrs_t* fcns = rct->fcns;
  p4est_t* p4est = rct->p4est;
  p4est_ghost_t* ghost = rct->ghost;
  element_data_t* ghost_data = rct->ghost_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = rct->dgmath_jit_dbase;

  /* double* x_temp = P4EST_ALLOC(double, vecs->local_nodes); */
  /* int i; */
  /* for(i = 0; i < vecs->local_nodes; i++){ */
    /* x_temp[i] = xx[i]; */
  /* } */
  
  problem_data_t vecs_for_res_build;
  problem_data_copy_ptrs(vecs, &vecs_for_res_build);
  vecs_for_res_build.u = (double*)xx;//x_temp;
  vecs_for_res_build.u0 = (double*)xx;//x_temp;
  vecs_for_res_build.Au = ftemp;
  
  fcns->build_residual(p4est, ghost, ghost_data, &vecs_for_res_build, dgmath_jit_dbase);

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
  weakeqn_ptrs_t* fcns = rct->fcns;
  p4est_t* p4est = rct->p4est;
  p4est_ghost_t* ghost = rct->ghost;
  element_data_t* ghost_data = rct->ghost_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = rct->dgmath_jit_dbase;


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

  /* for (i = 0; i < vecs->local_nodes; i++){ */
  
    /* printf("\n \n %f, %f\n \n", px[0], py[0]); */
  /* } */
  
  fcns->apply_lhs(p4est, ghost, ghost_data, &vecs_for_jac, dgmath_jit_dbase);

  ierr = VecRestoreArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecRestoreArray( y, &py ); CHKERRQ(ierr);
  return ierr;
}


void newton_petsc_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 p4est_ghost_t** ghost,
 element_data_t** ghost_data, 
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 krylov_petsc_params_t* krylov_params
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

  
  newton_petsc_ctx_t rct;
  rct.p4est = p4est;
  rct.vecs = vecs;
  rct.fcns = fcns;
  rct.ghost = *ghost;
  rct.ghost_data = *ghost_data;
  rct.dgmath_jit_dbase = dgmath_jit_dbase;
  SNESSetFunction(snes,r,newton_petsc_get_residual,(void*)&rct);//CHKERRQ(ierr);
  SNESGetKSP(snes,&ksp);
  /* PetscOptionsSetValue(NULL,"-snes_mf",""); */
  PetscOptionsSetValue(NULL,"-snes_converged_reason","");
  PetscOptionsSetValue(NULL,"-ksp_converged_reason","");
  /* PetscOptionsSetValue(NULL,"-ksp_monitor_true_residual",""); */
  PetscOptionsSetValue(NULL,"-snes_monitor","");
  /* PetscOptionsSetValue("-snes_monitor_solution",""); */
  /* PetscOptionsSetValue(NULL,"-pc_type","none"); */
  PetscOptionsSetValue(NULL,"-snes_atol","1e-50");
  PetscOptionsSetValue(NULL,"-ksp_atol","1e-50");
  PetscOptionsSetValue(NULL,"-ksp_rtol","1e-5");
  PetscOptionsSetValue(NULL,"-snes_rtol","1e-5");
  PetscOptionsSetValue(NULL,"-ksp_max_it","1000000");
  PetscOptionsSetValue(NULL,"-snes_linesearch_monitor","");
  PetscOptionsSetValue(NULL,"-ksp_view","");
  /* PetscOptionsSetValue(NULL,"-info","100000"); */

  /* absurdly large number */
  PetscOptionsSetValue(NULL,"-snes_max_funcs","10000000000000000");

  PetscOptionsSetValue(NULL,"-snes_linesearch_order","3");
  if (krylov_params->ksp_monitor == 1)
    PetscOptionsSetValue(NULL,"-ksp_monitor","");


  /* SET PRECONDITIONER */
  /* SET PRECONDITIONER */
  /* SET PRECONDITIONER */
  krylov_pc_ctx_t kct;
  kct.p4est = p4est;
  kct.vecs = vecs;
  kct.fcns = (void*)fcns;
  kct.ghost = ghost;
  kct.ghost_data = (void**)ghost_data;
  kct.dgmath_jit_dbase = dgmath_jit_dbase;

  /* PetscOptionsSetValue(NULL,"-ksp_monitor",""); */
  /* PetscOptionsSetValue(NULL,"-ksp_converged_reason",""); */
  /* PetscOptionsSetValue(NULL,"-ksp_atol","1e-20"); */
  /* PetscOptionsSetValue(NULL,"-with-debugging","1"); */
  /* PetscOptionsSetValue(NULL,"-ksp_rtol","1e-5"); */
  /* PetscOptionsSetValue(NULL,"-ksp_max_it","1000000"); */
  PC pc;
  KSPGetPC(ksp,&pc);
  
  krylov_pc_t* kp = NULL;
  if (krylov_params != NULL && krylov_params->user_defined_pc) {
    PCSetType(pc,PCSHELL);//CHKERRQ(ierr);
    kp = krylov_params->pc_create(&kct);
    PCShellSetApply(pc, krylov_petsc_pc_apply);//CHKERRQ(ierr);
    PCShellSetSetUp(pc, krylov_petsc_pc_setup);
    PCShellSetContext(pc, kp);//CHKERRQ(ierr);
  }
  else {
    PCSetType(pc,PCNONE);//CHKERRQ(ierr);
  }

  KSPSetType(ksp, krylov_params->ksp_type);
  KSPSetFromOptions(ksp);
  
  /* END SET PRECONDITIONER */
  /* END SET PRECONDITIONER */
  /* END SET PRECONDITIONER */

  
  /* PetscOptionsSetValue(NULL,"-ksp_view",""); */
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
