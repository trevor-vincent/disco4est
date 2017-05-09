#ifndef D4EST_CG_H
#define D4EST_CG_H 

#include <problem_data.h>
#include <problem_weakeqn_ptrs.h>

typedef struct {

  char input_section [50]; /* useful when using cg in multiple contexts */
  
  int monitor;
  int iter;
  double rtol;
  double atol;
  
} d4est_cg_params_t;
/* This file was automatically generated.  Do not edit! */

void d4est_cg_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 p4est_ghost_t** ghost,
 void** ghost_data, 
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geom,
 d4est_cg_params_t* params
);
void d4est_cg_input(p4est_t *p4est,const char *input_file,const char *input_section,const char *printf_prefix,d4est_cg_params_t *input);


#endif
