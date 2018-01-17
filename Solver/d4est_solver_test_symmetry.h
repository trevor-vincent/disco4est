#ifndef D4EST_SOLVER_TEST_SYMMETRY_H
#define D4EST_SOLVER_TEST_SYMMETRY_H 

#include <time.h>
#include <stdlib.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <d4est_element_data.h>
#include <d4est_elliptic_eqns.h>
#include <pXest.h>
#include <d4est_solver_test_symmetry.h>
#include <sc_reduce.h>
#include <sc_allgather.h>


typedef enum {SYM_PRINT_UNEQUAL_PAIRS_AND_XYZ,SYM_PRINT_UNEQUAL_PAIRS, SYM_PRINT_MAT_AND_TRANSPOSE_AS_VECS, SYM_PRINT_MAT}  d4est_solver_test_symmetry_print_option_t;

/* This file was automatically generated.  Do not edit! */
void d4est_solver_test_symmetry(p4est_t *p4est,p4est_ghost_t *ghost,d4est_element_data_t *ghost_data,int local_nodes,d4est_elliptic_eqns_t *fcns,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_solver_test_symmetry_print_option_t print,double sym_eps);


#endif
