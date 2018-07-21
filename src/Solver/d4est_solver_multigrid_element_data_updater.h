#ifndef D4EST_SOLVER_MULTIGRID_ELEMENT_DATA_UPDATER_H
#define D4EST_SOLVER_MULTIGRID_ELEMENT_DATA_UPDATER_H 

#include <d4est_solver_multigrid.h>

/* This file was automatically generated.  Do not edit! */
void d4est_solver_multigrid_element_data_updater_destroy(d4est_solver_multigrid_element_data_updater_t *updater,int num_of_levels);
void d4est_solver_multigrid_element_data_updater_update(p4est_t *p4est,int level,d4est_elliptic_data_t *vecs);
d4est_solver_multigrid_element_data_updater_t *d4est_solver_multigrid_element_data_updater_init(int num_of_levels,d4est_ghost_t *d4est_ghost_toplevel,d4est_ghost_data_t *d4est_ghost_data_toplevel,d4est_mesh_data_t *d4est_mesh_data_toplevel,void(*element_data_init_user_fcn)(d4est_element_data_t *,void *),d4est_mesh_initial_extents_t *initial_extents);


#endif
