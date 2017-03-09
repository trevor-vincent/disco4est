#ifndef KRYLOV_H
#define KRYLOV_H 

typedef struct {

  p4est_t* p4est;
  problem_data_t* vecs;
  void* fcns;
  p4est_ghost_t** ghost;
  void** ghost_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase;
  d4est_geometry_t* d4est_geom;
  
} petsc_ctx_t;

#endif
