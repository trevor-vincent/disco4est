#ifndef NEWTON_D4EST_H
#define NEWTON_D4EST_H 

typedef
void
(*krylov_solver_fcn_t)
(
 p4est_t*,
 d4est_elliptic_data_t*,
 d4est_elliptic_eqns_t*,
 p4est_ghost_t**,
 void**,
 d4est_operators_t*,
 d4est_geometry_t*,
 const char*,
 krylov_pc_t* 
);

int newton_d4est_solve
(
 p4est_t *p4est,
 d4est_elliptic_data_t *vecs,
 void *fcns,
 p4est_ghost_t **ghost,
 void **ghost_data,
 d4est_operators_t *d4est_ops,
 d4est_geometry_t *d4est_geom,
 const char *input_file,
 krylov_pc_t *krylov_pc
);


#endif
