#ifndef MATRIX_SYM_TESTER_H
#define MATRIX_SYM_TESTER_H 

#include <problem_data.h>
#include <problem_weakeqn_ptrs.h>
/* This file was automatically generated.  Do not edit! */
void
matrix_sym_tester_parallel_aux
(
 p4est_t* p4est,
 problem_data_t* vecs, 
 void* fcns,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 int i, /* local node on mpirank_i */
 int mpirank_i,
 int j, /* local node on mpirank_j */
 int mpirank_j, 
 double* Aji,
 double* Aij
);
int matrix_sym_tester_parallel(p4est_t *p4est,problem_data_t *vecs,void *fcns,p4est_ghost_t *ghost,void *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,int print,int num_tests,double sym_eps);
void
serial_matrix_sym_tester
(
 p4est_t* p4est,
 problem_data_t* vecs, /* only needed for # of nodes */
 weakeqn_ptrs_t* fcns,
 double sym_eps,
 d4est_operators_t* d4est_ops,
 int print,
 d4est_geometry_t* geom
);

/* void */
/* parallel_matrix_sym_tester */
/* ( */
/*  p4est_t* p4est, */
/*  problem_data_t* vecs, /\* only needed for # of nodes *\/ */
/*  void* fcns, */
/*  double sym_eps, */
/*  d4est_operators_t* d4est_ops, */
/*  int curved, */
/*  int print */
/* ); */

/* int */
/* matrix_sym_tester_parallel */
/* ( */
/*  p4est_t* p4est, */
/*  problem_data_t* vecs, /\* only needed for # of nodes *\/ */
/*  void* fcns, */
/*  int mpi_rank, */
/*  int num_tests, */
/*  double sym_eps, */
/*  d4est_operators_t* d4est_ops, */
/*  int curved, */
/*  int test_PD, /\* test if positive definite *\/ */
/*  int random, */
/*  int normalize */
/* ); */

int
matrix_spd_tester_parallel
(
 p4est_t* p4est,
 problem_data_t* vecs, 
 void* fcns,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 int print,
 int num_tests
);
#endif
