#ifndef MULTIGRID_ELEMENT_DATA_UPDATE_NOCURVED_H
#define MULTIGRID_ELEMENT_DATA_UPDATE_NOCURVED_H 

#include <multigrid.h>

/* This file was automatically generated.  Do not edit! */
void multigrid_element_data_updater_nocurved_destroy(multigrid_element_data_updater_t *updater);
void multigrid_element_data_updater_nocurved_update(p4est_t *p4est,int level,d4est_elliptic_problem_data_t *vecs);
multigrid_element_data_updater_t *multigrid_element_data_updater_nocurved_init(p4est_ghost_t **ghost,element_data_t **ghost_data);

#endif
