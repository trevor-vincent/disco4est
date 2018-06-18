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
#include <d4est_solver_full_matrix.h>
#include <sc_reduce.h>
#include <sc_allgather.h>

void
d4est_solver_full_matrix
(
 p4est_t* p4est,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_elliptic_eqns_t* fcns,
 d4est_elliptic_data_t* vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 double* a_mat
)
{
  int local_nodes = vecs->local_nodes;
  double* u_temp = P4EST_ALLOC(double, local_nodes);
  double* Au_temp = P4EST_ALLOC(double, local_nodes);

  
  d4est_elliptic_data_t vecs_temp;
  d4est_elliptic_data_copy_ptrs(vecs, &vecs_temp);
  vecs_temp.u = u_temp;
  vecs_temp.Au = Au_temp;
  
  d4est_util_fill_array(vecs_temp.u, 0., vecs_temp.local_nodes);
  
  for (int i = 0; i < vecs_temp.local_nodes; i++){
    vecs_temp.u[i] = 1.;

  d4est_elliptic_eqns_apply_lhs
    (
     p4est,
     ghost,
     ghost_data,
     fcns,
     &vecs_temp,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors
    );

    
    
    d4est_linalg_set_column(a_mat, vecs_temp.Au, i, vecs_temp.local_nodes, vecs_temp.local_nodes);
    vecs_temp.u[i] = 0.;
  }

  
  P4EST_FREE(u_temp);
  P4EST_FREE(Au_temp);
}
