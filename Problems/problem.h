#ifndef PROBLEM_H
#define PROBLEM_H 

#include "../dGMath/dgmath.h"
#include "../pXest/pXest.h"

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
 int argc,
 char* argv [],
 p4est_t* p4est,
 p4est_geometry_t* p4est_geom,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 int proc_size,
 sc_MPI_Comm mpicomm,
 int load_from_checkpoint
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
