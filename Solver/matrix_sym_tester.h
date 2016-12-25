#ifndef MATRIX_SYM_TESTER_H
#define MATRIX_SYM_TESTER_H 

#include <problem_data.h>
#include <problem_weakeqn_ptrs.h>

void
serial_matrix_sym_tester
(
 p4est_t* p4est,
 problem_data_t* vecs, /* only needed for # of nodes */
 void* fcns,
 double sym_eps,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 int curved,
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
/*  dgmath_jit_dbase_t* dgmath_jit_dbase, */
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
/*  dgmath_jit_dbase_t* dgmath_jit_dbase, */
/*  int curved, */
/*  int test_PD, /\* test if positive definite *\/ */
/*  int random, */
/*  int normalize */
/* ); */


#endif
