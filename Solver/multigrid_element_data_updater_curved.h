#ifndef MULTIGRID_ELEMENT_DATA_UPDATER_CURVED_H
#define MULTIGRID_ELEMENT_DATA_UPDATER_CURVED_H 

#include <multigrid.h>

void multigrid_element_data_updater_curved_update(p4est_t *p4est,int level,problem_data_t *vecs);
void multigrid_element_data_updater_curved_destroy(multigrid_element_data_updater_t *updater,int num_of_levels);
multigrid_element_data_updater_t *multigrid_element_data_updater_curved_init(int num_of_levels,p4est_ghost_t **ghost,d4est_element_data_t **ghost_data,d4est_mesh_geometry_storage_t *d4est_mesh_geometry_storage_toplevel,d4est_geometry_t *d4est_geom,void(*element_data_init_user_fcn)(void *,void *),void *user);



#endif
