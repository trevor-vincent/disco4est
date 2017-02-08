#ifndef PROBLEM_H
#define PROBLEM_H 

#include <d4est_geometry.h>
#include <dgmath.h>
#include <pXest.h>

p4est_connectivity_t*
problem_build_conn();

p4est_t*
problem_build_p4est
(
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t* conn,
 p4est_locidx_t min_quadrants,
 int min_level, 
 int fill_uniform
);

p4est_geometry_t*
problem_build_geom
(
 p4est_connectivity_t* conn
);

void
problem_init
(
 const char* input_file,
 p4est_t* p4est,
 d4est_geometry_t* d4est_geom,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 int proc_size,
 sc_MPI_Comm mpicomm
);


p4est_t*
problem_load_p4est_from_checkpoint
(
 const char* filename,
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t** conn
);


void
problem_help();

#endif
