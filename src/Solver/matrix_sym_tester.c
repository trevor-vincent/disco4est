/**
 * @file   matrix_sym_tester.c
 * @author Trevor Vincent <tvincent@cita.utoronto.ca>
 * @date   Sat Nov 21 18:49:21 2015
 * 
 * @brief  Use only for single process runs.
 * 
 * 
 */

#include <time.h>
#include <stdlib.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <d4est_solver_symmetry_test.h>
#include <sc_reduce.h>
#include <sc_allgather.h>
