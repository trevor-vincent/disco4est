#ifndef D4EST_SOLVER_SCHWARZ_TRANSFER_GHOST_DATA_H
#define D4EST_SOLVER_SCHWARZ_TRANSFER_GHOST_DATA_H 

#include <pXest.h>
#include <d4est_mesh.h>
#include <d4est_element_data.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_solver_schwarz_transfer_ghost_data.h>

/* This file was automatically generated.  Do not edit! */
void d4est_solver_schwarz_transfer_ghost_data_and_add_corrections(p4est_t *p4est,d4est_ghost_t *d4est_ghost,d4est_solver_schwarz_metadata_t *schwarz_metadata,d4est_ghost_data_ext_t **d4est_ghost_data_ext,double *u,double *du_over_subdomains);

#endif
