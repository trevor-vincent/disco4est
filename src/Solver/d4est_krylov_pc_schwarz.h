#ifndef D4EST_KRYLOV_PC_SCHWARZ_H
#define D4EST_KRYLOV_PC_SCHWARZ_H 

#include <pXest.h>
#include <d4est_solver_petsc.h>
#include <d4est_solver_schwarz.h>
#include <d4est_krylov_pc.h>

typedef struct {

  d4est_solver_schwarz_t* schwarz;
  int iterations;
  int verbose;
  void(*user_setup_fcn)(d4est_krylov_pc_t* kpc);
  
} d4est_krylov_pc_schwarz_data_t;



/* This file was automatically generated.  Do not edit! */
d4est_krylov_pc_t *d4est_krylov_pc_schwarz_create(const char *input_file,d4est_solver_schwarz_t *schwarz,void(*user_setup_fcn)(d4est_krylov_pc_t *kpc));
void d4est_krylov_pc_schwarz_apply(d4est_krylov_pc_t *kpc,double *xp,double *yp);
void d4est_krylov_pc_schwarz_destroy(d4est_krylov_pc_t *kpc);
void d4est_krylov_pc_schwarz_setup(d4est_krylov_pc_t *kpc);

#endif
