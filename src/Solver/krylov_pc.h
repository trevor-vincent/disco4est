#ifndef KRYLOV_PC_H
#define KRYLOV_PC_H 

#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_petsc.h>

typedef struct krylov_pc krylov_pc_t;

typedef
void
(*krylov_pc_apply_fcn_t)
(
 krylov_pc_t*,
 double*,
 double*
);

typedef
void
(*krylov_pc_setup_fcn_t)
(
 krylov_pc_t*
);

struct krylov_pc {

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
  
  krylov_ctx_t* pc_ctx;

  /* End INTERNAL Parameters */
  
};


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
    kp->pc_setup(kp);
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
