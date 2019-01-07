#ifndef PROBLEM_H
#define PROBLEM_H 

#include <d4est_geometry.h>
#include <d4est_element_data.h>
#include <d4est_operators.h>
#include <d4est_mesh.h>
#include <pXest.h>

void problem_init(p4est_t *p4est,d4est_ghost_t **ghost, d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *geometric_factors,d4est_mesh_initial_extents_t *initial_extents,const char *input_file, sc_MPI_Comm mpicomm);

void
problem_help();


#endif
