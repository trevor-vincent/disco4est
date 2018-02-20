#ifndef KRYLOV_H
#define KRYLOV_H 

typedef struct {

  p4est_t* p4est;
  d4est_elliptic_data_t* vecs;
  void* fcns;
  p4est_ghost_t** ghost;
  void** ghost_data;
  d4est_operators_t* d4est_ops;
  d4est_geometry_t* d4est_geom;
  
} krylov_ctx_t;

#endif
