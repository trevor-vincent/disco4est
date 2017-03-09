#ifndef KRYLOV_PC_H
#define KRYLOV_PC_H 

#include "../pXest/pXest.h"
#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"
#include <d4est_petsc.h>

typedef
void
(*krylov_pc_apply_fcn_t)
(
 void*,
 double*,
 double*
);

typedef
void
(*krylov_pc_setup_fcn_t)
(
 void*
);

typedef struct {

  /* Start EXTERNAL Parameters */
  /* THESE MUST BE SET */

  /* Cannot be NULL */
  krylov_pc_apply_fcn_t pc_apply;

  /* can be NULL */
  krylov_pc_setup_fcn_t pc_setup;

  /* can be NULL */
  void* pc_data;

  /* End EXTERNAL Parameters */
  
  /* Start INTERNAL Parameters */
  /* DO NOT SET THESE */
  
  petsc_ctx_t* pc_ctx;

  /* End INTERNAL Parameters */
  
} krylov_pc_t;


static
PetscErrorCode
krylov_petsc_pc_setup
(
 PC pc
)
{
  PetscErrorCode ierr;
  krylov_pc_t* kp;
  ierr = PCShellGetContext(pc,(void**)&kp);CHKERRQ(ierr);
  if (kp->pc_setup != NULL){
    kp->pc_setup(kp->pc_ctx);
  }
  return ierr;
}

static
PetscErrorCode
krylov_petsc_pc_apply(PC pc,Vec x,Vec y)
{
  PetscErrorCode ierr;
  krylov_pc_t* kp;
  ierr = PCShellGetContext(pc,(void**)&kp);CHKERRQ(ierr);

  double* yp;
  const double* xp;
  
  VecGetArray(y,&yp); 
  VecGetArrayRead(x,&xp);

  kp->pc_apply(kp, (double*)xp, yp);
  
  VecRestoreArray(y,&yp);//CHKERRQ(ierr);
  VecRestoreArrayRead(x,&xp);
  
  return ierr;
}


#endif
