#ifndef MULTIGRID_ELEMENT_DATA_UPDATER_H
#define MULTIGRID_ELEMENT_DATA_UPDATER_H 

#include <multigrid.h>

/* This file was automatically generated.  Do not edit! */
void multigrid_element_data_updater_destroy(multigrid_element_data_updater_t *updater,int num_of_levels);
void multigrid_element_data_updater_update(p4est_t *p4est,int level,d4est_elliptic_data_t *vecs);
multigrid_element_data_updater_t *multigrid_element_data_updater_init(int num_of_levels,p4est_ghost_t **ghost,d4est_element_data_t **ghost_data,d4est_mesh_data_t *d4est_mesh_data_toplevel,void(*element_data_init_user_fcn)(d4est_element_data_t *,void *),void *user);



#endif
