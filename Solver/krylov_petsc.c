#include "../pXest/pXest.h"
#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../Solver/krylov_petsc.h"
#include "petscsnes.h"
#include <krylov_pc.h>
#include <krylov_petsc_pc.h>

static
PetscErrorCode krylov_petsc_apply_aij( Mat A, Vec x, Vec y )
{
  void           *ctx;
  PetscErrorCode ierr;

  krylov_pc_ctx_t* kct;
  const double* px;
  double* py;

  /* PetscFunctionBegin; */
  ierr = MatShellGetContext( A, &ctx ); CHKERRQ(ierr);  
  kct = (krylov_pc_ctx_t *)ctx;
  ierr = VecGetArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecGetArray( y, &py ); CHKERRQ(ierr);

  /* weakeqn_ptrs_t* fcns = kct->fcns; */
  p4est_t* p4est = kct->p4est;
  p4est_ghost_t* ghost = *kct->ghost;
  dgmath_jit_dbase_t* dgmath_jit_dbase = kct->dgmath_jit_dbase;

  problem_data_t vecs_for_aij;
  problem_data_copy_ptrs(kct->vecs, &vecs_for_aij);

  vecs_for_aij.u = (double*)px;
  vecs_for_aij.Au = py;

  if (kct->d4est_geom == NULL){
    ((weakeqn_ptrs_t*)(kct->fcns))->apply_lhs(p4est,
                                              ghost,
                                              (element_data_t*)(*kct->ghost_data),
                                              &vecs_for_aij,
                                              dgmath_jit_dbase);
  }
  else {
    ((curved_weakeqn_ptrs_t*)(kct->fcns))->apply_lhs(p4est,
                                                     ghost,
                                                     (curved_element_data_t*)(*kct->ghost_data),
                                                     &vecs_for_aij,
                                                     dgmath_jit_dbase,
                                                     kct->d4est_geom);
  }
  ierr = VecRestoreArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecRestoreArray( y, &py ); CHKERRQ(ierr);
  return ierr;
}
krylov_petsc_info_t
krylov_petsc_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 void* fcns,
 p4est_ghost_t** ghost,
 void** ghost_data, 
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geom,
 krylov_petsc_params_t* krylov_params
)
{
  krylov_petsc_info_t info;
  KSP ksp;
  Vec x,b;
  PC             pc;
  /* double* u_temp; */
  /* double* rhs_temp; */

  krylov_pc_ctx_t kct;
  kct.p4est = p4est;
  kct.vecs = vecs;
  kct.fcns = fcns;
  kct.ghost = ghost;
  kct.ghost_data = ghost_data;
  kct.dgmath_jit_dbase = dgmath_jit_dbase;
  kct.d4est_geom = d4est_geom;
  kct.pc_data = krylov_params->pc_data;

  int local_nodes = vecs->local_nodes;
  double* u = vecs->u;
  double* rhs = vecs->rhs;

  
  KSPCreate(PETSC_COMM_WORLD,&ksp);  
  VecCreate(PETSC_COMM_WORLD,&x);//CHKERRQ(ierr);
  VecSetSizes(x, local_nodes, PETSC_DECIDE);//CHKERRQ(ierr);
  VecSetFromOptions(x);//CHKERRQ(ierr);
  VecDuplicate(x,&b);//CHKERRQ(ierr);


  
  if (krylov_params->ksp_monitor)
    PetscOptionsSetValue(NULL,"-ksp_monitor","");
  if (krylov_params->ksp_monitor)
    PetscOptionsSetValue(NULL,"-ksp_view","");
    /* KSPMonitorSet(ksp, KSPMonitorDefault, NULL, NULL); */

  /* PetscOptionsSetValue(NULL,"-ksp_converged_reason",""); */
  PetscOptionsSetValue(NULL,"-ksp_atol","1e-18");
  /* PetscOptionsSetValue(NULL,"-with-debugging","1"); */
  PetscOptionsSetValue(NULL,"-ksp_rtol","1e-20");
  PetscOptionsSetValue(NULL,"-ksp_max_it","1000000");

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

  KSPSetType(ksp, krylov_params->krylov_type);
  KSPSetFromOptions(ksp);

  /* Create matrix-free shell for Aij */
  Mat A;
  MatCreateShell
    (
     PETSC_COMM_WORLD,
     local_nodes,
     local_nodes,
     PETSC_DETERMINE,
     PETSC_DETERMINE,
     (void*)&kct,
     &A
    ); 
  MatShellSetOperation(A,MATOP_MULT,(void(*)())krylov_petsc_apply_aij);

  /* Set Amat and Pmat, where Pmat is the matrix the Preconditioner needs */
  KSPSetOperators(ksp,A,A);

  /* linalg_copy_1st_to_2nd(vecs->u, u_temp, local_nodes); */
  /* linalg_copy_1st_to_2nd(vecs->rhs, rhs_temp, local_nodes); */

  VecPlaceArray(b, rhs);
  VecPlaceArray(x, u);
  
  KSPSolve(ksp,b,x);

  if (krylov_params != NULL && krylov_params->user_defined_pc) {
    krylov_params->pc_destroy(kp);
  }
  
  KSPGetIterationNumber(ksp, &(info.iterations));
  KSPGetResidualNorm(ksp, &(info.residual_norm));
  
  /* linalg_copy_1st_to_2nd(u_temp, vecs->u, local_nodes); */

  /* VecRestoreArray(x,&u_temp); */
  /* VecRestoreArray(b,&rhs_temp); */
  VecResetArray(b);
  VecResetArray(x);
  VecDestroy(&x);
  VecDestroy(&b);
  KSPDestroy(&ksp);

  return info;
}
