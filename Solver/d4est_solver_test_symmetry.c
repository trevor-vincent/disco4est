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
#include <d4est_element_data.h>
#include <d4est_mesh.h>
#include <d4est_elliptic_eqns.h>
#include <pXest.h>
#include <d4est_solver_test_symmetry.h>
#include <sc_reduce.h>
#include <sc_allgather.h>

void
d4est_solver_test_symmetry
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 int local_nodes,
 d4est_elliptic_eqns_t* fcns,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_geometry_storage_t* d4est_factors,
 d4est_solver_test_symmetry_print_option_t print,
 double sym_eps
)
{
  double* a_mat = P4EST_ALLOC(double, local_nodes*local_nodes);
  double* a_mat_trans = P4EST_ALLOC(double, local_nodes*local_nodes);

  double* u_temp = P4EST_ALLOC(double, local_nodes);
  double* Au_temp = P4EST_ALLOC(double, local_nodes);

  d4est_elliptic_data_t vecs;
  vecs.u = u_temp;
  vecs.Au = Au_temp;
  vecs.local_nodes = local_nodes;
  
  d4est_linalg_fill_vec(vecs.u, 0., vecs.local_nodes);
  
  for (int i = 0; i < local_nodes; i++){
    vecs.u[i] = 1.;

    d4est_elliptic_eqns_apply_lhs
      (
     p4est,
     ghost,
     ghost_data,
     fcns,
     &vecs,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors
      );
    
    d4est_linalg_set_column(a_mat, vecs.Au, i, vecs.local_nodes, vecs.local_nodes);
    vecs.u[i] = 0.;
  }

  d4est_linalg_mat_transpose(a_mat, a_mat_trans, local_nodes);
  int compare = d4est_util_compare_vecs(a_mat, a_mat_trans, local_nodes*local_nodes, sym_eps);

  double biggest_err = -1.;
  int biggest_id = -1;

  d4est_util_find_biggest_error(a_mat, a_mat_trans, local_nodes*local_nodes, &biggest_err, &biggest_id);
  
  if (compare == 1){
    printf("A_{DG} is symmetric\n");
    printf("biggest error is at %d with size %.20f\n", biggest_id, biggest_err);
  }
  else{
    printf("A_{DG} is NOT symmetric\n");
    printf("biggest error is at %d with size %.20f\n", biggest_id, biggest_err);
  }
  
  if (print == SYM_PRINT_MAT)
    d4est_util_print_matrix(a_mat, local_nodes, local_nodes, "A_{DG} = ",0);

  if (print == SYM_PRINT_MAT_AND_TRANSPOSE_AS_VECS){
    printf("Printing vectors of %d length\n", local_nodes*local_nodes);
    for (int i = 0; i < local_nodes; i++){
      for (int j = 0; j < local_nodes; j++)
        printf("%d,%d A A^t: %.20f %.20f\n", i,j, a_mat[i*local_nodes + j], a_mat_trans[i*local_nodes + j]);
    }
  }

  if (print == SYM_PRINT_UNEQUAL_PAIRS){
    /* D4EST_ASSERT(curved); */
    printf("(i,j) Pairs that aren't equal\n");
    for (int i = 0; i < local_nodes; i++){
      for (int j = 0; j <= i; j++){
        if (fabs(a_mat[i*local_nodes + j] - a_mat_trans[i*local_nodes + j]) > sym_eps) {
          int e1 = d4est_mesh_debug_find_node(p4est, i);
          int e2 = d4est_mesh_debug_find_node(p4est, j);
          printf("node %d in element %d, node %d in element %d, A/A^t/ERR: %.20f %.20f %.20f\n", i, e1, j, e2, a_mat[i*local_nodes + j], a_mat_trans[i*local_nodes + j], fabs(a_mat[i*local_nodes + j] - a_mat_trans[i*local_nodes + j]));
        }
      }
    }
  }

  if (print == SYM_PRINT_UNEQUAL_PAIRS_AND_XYZ){
    printf("(i,j) Pairs that aren't equal\n");
    for (int i = 0; i < local_nodes; i++){
      for (int j = 0; j <= i; j++){
        if (fabs(a_mat[i*local_nodes + j] - a_mat_trans[i*local_nodes + j]) > sym_eps) {
          int e1 = d4est_mesh_debug_find_node(p4est, i);
          int e2 = d4est_mesh_debug_find_node(p4est, j);
          d4est_element_data_t* e1_data = d4est_mesh_get_element_data(p4est, e1);
          d4est_element_data_t* e2_data = d4est_mesh_get_element_data(p4est, e2);
          printf("node %d in element %d (x,y,z = %f,%f,%f), node %d in element %d (x,y,z = %f,%f,%f), A/A^t/ERR: %.20f %.20f %.20f\n", i, e1, e1_data->xyz[0][0],e1_data->xyz[1][0],  0., j, e2, e2_data->xyz[0][0],e2_data->xyz[1][0], 0., a_mat[i*local_nodes + j], a_mat_trans[i*local_nodes + j], fabs(a_mat[i*local_nodes + j] - a_mat_trans[i*local_nodes + j]));
        }
      }
    }
  }
  
  P4EST_FREE(u_temp);
  P4EST_FREE(Au_temp);
  P4EST_FREE(a_mat);
  P4EST_FREE(a_mat_trans);
}
