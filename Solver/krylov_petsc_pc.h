#ifndef KRYLOV_PETSC_PC_H
#define KRYLOV_PETSC_PC_H 

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
  kp->pc_setup(kp->pc_ctx);
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

  kp->pc_apply(kp->pc_ctx, (double*)xp, yp);
  
  VecRestoreArray(y,&yp);//CHKERRQ(ierr);
  VecRestoreArrayRead(x,&xp);
  
  return ierr;
}


#endif
