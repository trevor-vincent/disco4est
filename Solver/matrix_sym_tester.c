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
#include "../Solver/matrix_sym_tester.h"
#include "../LinearAlgebra/d4est_linalg.h"
#include "../Utilities/util.h"
#include <sc_reduce.h>
#include <sc_allgather.h>

/* #define NON_RANDOM */


/* void */
/* curved_matrix_sym_tester */
/* ( */
/*  p4est_t* p4est, */
/*  d4est_elliptic_data_t* vecs, /\* only needed for # of nodes *\/ */
/*  curved_d4est_elliptic_eqns_t* fcns, */
/*  int print, */
/*  d4est_operators_t* d4est_ops */
/* ) */
/* { */
/*   d4est_element_data_t *ghost_data; */
/*   p4est_ghost_t *ghost; */
/*   int i, local_nodes; */
/*   local_nodes = vecs->local_nodes; */
/*   double* a_mat = P4EST_ALLOC(double, local_nodes*local_nodes); */
/*   double* a_mat_trans = P4EST_ALLOC(double, local_nodes*local_nodes); */
/*   ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE); */
/*   ghost_data = P4EST_ALLOC (d4est_element_data_t, ghost->ghosts.elem_count); */

/*   double* u_temp = P4EST_ALLOC(double, local_nodes); */
/*   double* tmp = vecs->u; */
/*   vecs->u = u_temp; */
  
/*   d4est_linalg_fill_vec(vecs->u, 0., vecs->local_nodes); */
  
/*   for (i = 0; i < vecs->local_nodes; i++){ */
/*     vecs->u[i] = 1.; */
/*     fcns->apply_lhs(p4est, ghost, ghost_data, vecs, d4est_ops); */
/*     d4est_linalg_set_column(a_mat, vecs->Au, i, vecs->local_nodes, vecs->local_nodes); */
/*     vecs->u[i] = 0.; */
/*   } */

/*   d4est_linalg_mat_transpose(a_mat, a_mat_trans, local_nodes); */
/*   int compare = util_compare_vecs(a_mat, a_mat_trans, local_nodes, .000001); */

/*   if (compare == 1) */
/*     printf("A_{DG} is symmetric\n"); */
/*   else */
/*     printf("A_{DG} is NOT symmetric\n"); */

/*   if (print == 1) */
/*     util_print_matrix(a_mat, local_nodes, local_nodes, "A_{DG} = ",0); */

/*   if (print == 2){ */
/*     printf("Printing vectors of %d length\n", local_nodes*local_nodes); */
/*     util_print_matrices(a_mat, a_mat_trans, local_nodes*local_nodes, 1, "A, A^T = "); */
/*   } */
    
/*   vecs->u = tmp; */
/*   P4EST_FREE(u_temp); */
/*   P4EST_FREE(a_mat); */
/*   P4EST_FREE(a_mat_trans); */
/*   P4EST_FREE(ghost_data); */
/*   p4est_ghost_destroy (ghost); */
/* } */

/* /\**  */
/*  * Constructs a few random vectors and tests */
/*  * if v^T A u = u ^ T A v in parallel */
/*  * can be used on large matrices */
/*  *  */
/*  * @param p4est  */
/*  * @param vecs  */
/*  * @param fcns  */
/*  * @param print  */
/*  *\/ */
/* int */
/* matrix_sym_tester_parallel */
/* ( */
/*  p4est_t* p4est, */
/*  d4est_elliptic_data_t* vecs, /\* only needed for # of nodes *\/ */
/*  d4est_elliptic_eqns_t* fcns, */
/*  int mpi_rank, */
/*  int tests, */
/*  int test_PD, /\* test if positive definite *\/ */
/*  double sym_eps, */
/*  int normalize, */
/*  d4est_operators_t* d4est_ops */
/* ) */
/* { */
/*   element_data_t *ghost_data; */
/*   p4est_ghost_t *ghost;   */
/*   ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL); */
/*   ghost_data = P4EST_ALLOC (element_data_t, ghost->ghosts.elem_count); */

/*   double* tmp = vecs->u; */

/*   double* u = P4EST_ALLOC(double, vecs->local_nodes); */
/*   double* v = P4EST_ALLOC(double, vecs->local_nodes); */

/*   int compare_dot = 1; */
/*   double vTAu_local, vTAu_global; */
/*   double uTAv_local, uTAv_global; */

  
/*   /\* srand(time(NULL)); *\/ */
/*   srand(134434); */
/*   int i; */

/*   for (i = 0; i < tests; i++){ */

/*     vTAu_global = 0.; */
/*     uTAv_global = 0.; */
    
/*     /\* generate two random vectors *\/ */

/* #ifndef NON_RANDOM */
/*     for (j = 0; j < vecs->local_nodes; j++){ */
/*       u[j] = (double)rand() / (double)RAND_MAX; */
/*     } */
/*     for (j = 0; j < vecs->local_nodes; j++){ */
/*       v[j] = (double)rand() / (double)RAND_MAX; */
/*     } */
/* #else */
/*     /\* for testing multi-process runs *\/ */
/*     double rand1 = 1.32*(i+1); */
/*     double rand2 = 1.453*(i+1); */
    
/*     d4est_linalg_fill_vec(u, rand1, vecs->local_nodes); */
/*     d4est_linalg_fill_vec(v, rand2, vecs->local_nodes); */
/* #endif    */

/*     if (normalize == 1){ */
/*       d4est_linalg_vec_normalize(u, vecs->local_nodes); */
/*       d4est_linalg_vec_normalize(v, vecs->local_nodes); */
/*     } */
    
/*     vecs->u = u; */
/*     double* Au = vecs->Au; */
/*     fcns->apply_lhs(p4est, ghost, ghost_data, vecs, d4est_ops); */

/*     /\* util_print_matrices(vecs->u, vecs->Au, vecs->local_nodes, 1, "u, Au = "); *\/ */

    
/*     vTAu_local = d4est_linalg_vec_dot(v, Au, vecs->local_nodes); */

/*     sc_reduce */
/*       ( */
/*        &vTAu_local, */
/*        &vTAu_global, */
/*        1, */
/*        sc_MPI_DOUBLE, */
/*        sc_MPI_SUM, */
/*        0, */
/*        sc_MPI_COMM_WORLD */
/*       ); */
    
/*     vecs->u = v; */
/*     double* Av = vecs->Au; */
/*     fcns->apply_lhs(p4est, ghost, ghost_data, vecs, d4est_ops); */
/*     uTAv_local = d4est_linalg_vec_dot(u, Av, vecs->local_nodes); */

/*     sc_reduce */
/*       ( */
/*        &uTAv_local, */
/*        &uTAv_global, */
/*        1, */
/*        sc_MPI_DOUBLE, */
/*        sc_MPI_SUM, */
/*        0, */
/*        sc_MPI_COMM_WORLD */
/*       ); */


/*     if (mpi_rank == 0){ */
/*       compare_dot *= (fabs(uTAv_global - vTAu_global) < sym_eps); */
/*       printf("[MATRIX_SYM_TESTER]: %.20f %.20f\n", uTAv_global, vTAu_global); */
/*       if (test_PD == 1){ */
/*         compare_dot *= (uTAv_global > 0); */
/*         compare_dot *= (vTAu_global > 0); */
/*       } */
/*     } */
/*     /\* sc_MPI_Barrier(sc_MPI_COMM_WORLD); *\/ */
/*   } */

/*   if (mpi_rank == 0){ */
/*     if (test_PD == 1){ */
/*       if (compare_dot == 1) */
/*         printf("[MATRIX_SYM_TESTER]: SPD TEST PASSED\n"); */
/*       else */
/*         printf("[MATRIX_SYM_TESTER]: SPD TEST FAILED\n"); */
/*     } */
/*     else{ */
/*       if (compare_dot == 1) */
/*         printf("[MATRIX_SYM_TESTER]: SYM TEST PASSED\n"); */
/*       else */
/*         printf("[MATRIX_SYM_TESTER]: SYM TEST FAILED\n"); */
/*     } */
/*   } */
  
/*   vecs->u = tmp; */
  
/*   P4EST_FREE(u); */
/*   P4EST_FREE(v); */
/*   /\* P4EST_FREE(Au); *\/ */
/*   P4EST_FREE(ghost_data); */
/*   p4est_ghost_destroy (ghost); */

/*   return compare_dot; */
/* } */

/* int */
/* curved_matrix_sym_tester_parallel */
/* ( */
/*  p4est_t* p4est, */
/*  d4est_elliptic_data_t* vecs, /\* only needed for # of nodes *\/ */
/*  curved_d4est_elliptic_eqns_t* fcns, */
/*  int mpi_rank, */
/*  int tests, */
/*  int test_PD, /\* test if positive definite *\/ */
/*  double sym_eps, */
/*  int random, */
/*  int normalize, */
/*  d4est_operators_t* d4est_ops */
/* ) */
/* { */
/*   d4est_element_data_t *ghost_data; */
/*   p4est_ghost_t *ghost;   */
/*   ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE); */
/*   ghost_data = P4EST_ALLOC (d4est_element_data_t, ghost->ghosts.elem_count); */

/*   double* tmp = vecs->u; */

/*   double* u = P4EST_ALLOC(double, vecs->local_nodes); */
/*   double* v = P4EST_ALLOC(double, vecs->local_nodes); */

/*   int compare_dot = 1; */
/*   double vTAu_local, vTAu_global; */
/*   double uTAv_local, uTAv_global; */

  
/*   /\* srand(time(NULL)); *\/ */
/*   srand(134434); */
/*   int i; */

/*   for (i = 0; i < tests; i++){ */

/*     vTAu_global = 0.; */
/*     uTAv_global = 0.; */
    
/*     /\* generate two random vectors *\/ */

/*     if (random == 1){ */
/*       for (j = 0; j < vecs->local_nodes; j++) */
/*         u[j] = (double)rand() / (double)RAND_MAX; */
/*       for (j = 0; j < vecs->local_nodes; j++) */
/*         v[j] = (double)rand() / (double)RAND_MAX; */
/*     } */
/*     else { */
/*       double rand1 = 1.32*(i+1); */
/*       double rand2 = 1.453*(i+1); */
/*       d4est_linalg_fill_vec(u, rand1, vecs->local_nodes); */
/*       d4est_linalg_fill_vec(v, rand2, vecs->local_nodes); */
/*     } */

/*     if (normalize == 1){ */
/*       d4est_linalg_vec_normalize(u, vecs->local_nodes); */
/*       d4est_linalg_vec_normalize(v, vecs->local_nodes); */
/*     } */
    
    
/*     vecs->u = u; */
/*     /\* util_print_matrix(vecs->u, vecs->local_nodes, 1, "u = ", 0); *\/ */

/*     /\* printf(" u = {\n"); *\/ */
/*     /\* for (i = 0; i < vecs->local_nodes; i++){ *\/ */
/*     /\*   printf("%f, \n", vecs->u[i]); *\/ */
/*     /\* } *\/ */
/*     /\* printf("}\n"); *\/ */

/*     /\* d4est_linalg_fill_vec(vecs->Au, 0., vecs->local_nodes); *\/ */
    
/*     double* Au = vecs->Au; */
/*     fcns->apply_lhs(p4est, ghost, ghost_data, vecs, d4est_ops); */

/*     /\* util_print_matrices(vecs->u, vecs->Au, vecs->local_nodes, 1, "u, Au = "); *\/ */

/*     vTAu_local = d4est_linalg_vec_dot(v, Au, vecs->local_nodes); */

/*     sc_reduce */
/*       ( */
/*        &vTAu_local, */
/*        &vTAu_global, */
/*        1, */
/*        sc_MPI_DOUBLE, */
/*        sc_MPI_SUM, */
/*        0, */
/*        sc_MPI_COMM_WORLD */
/*       ); */
    
/*     vecs->u = v; */
/*     double* Av = vecs->Au; */
/*     fcns->apply_lhs(p4est, ghost, ghost_data, vecs, d4est_ops); */
/*     uTAv_local = d4est_linalg_vec_dot(u, Av, vecs->local_nodes); */

    
/*     /\* printf("\ni = %d\n", i); *\/ */
/*     /\* int f; *\/ */
/*     /\* for (f = 0; f < vecs->local_nodes; f++){ *\/ */
/*     /\*   printf("u, Av, v = %f, %f, %f\n", u[f], Av[f], v[f]); *\/ */
/*     /\* } *\/ */
    
/*     sc_reduce */
/*       ( */
/*        &uTAv_local, */
/*        &uTAv_global, */
/*        1, */
/*        sc_MPI_DOUBLE, */
/*        sc_MPI_SUM, */
/*        0, */
/*        sc_MPI_COMM_WORLD */
/*       ); */

/*    /\* printf(" v = {\n"); *\/ */
/*    /\*  for (f = 0; f < vecs->local_nodes; f++){ *\/ */
/*    /\*    printf("%f, \n", v[f]); *\/ */
/*    /\*  } *\/ */
/*    /\*  printf("}\n"); *\/ */
    

/*     if (mpi_rank == 0){ */
/*       compare_dot *= (fabs(uTAv_global - vTAu_global) < sym_eps); */
/*       printf("[MATRIX_SYM_TESTER]: %.20f %.20f\n", uTAv_global, vTAu_global); */
/*       if (test_PD == 1){ */
/*         compare_dot *= (uTAv_global > 0); */
/*         compare_dot *= (vTAu_global > 0); */
/*       } */
/*     } */
/*     /\* sc_MPI_Barrier(sc_MPI_COMM_WORLD); *\/ */
/*   } */

/*   if (mpi_rank == 0){ */
/*     if (test_PD == 1){ */
/*       if (compare_dot == 1) */
/*         printf("[MATRIX_SYM_TESTER]: SPD TEST PASSED\n"); */
/*       else */
/*         printf("[MATRIX_SYM_TESTER]: SPD TEST FAILED\n"); */
/*     } */
/*     else{ */
/*       if (compare_dot == 1) */
/*         printf("[MATRIX_SYM_TESTER]: SYM TEST PASSED\n"); */
/*       else */
/*         printf("[MATRIX_SYM_TESTER]: SYM TEST FAILED\n"); */
/*     } */
/*   } */
  
/*   vecs->u = tmp; */
  
/*   P4EST_FREE(u); */
/*   P4EST_FREE(v); */
/*   /\* P4EST_FREE(Au); *\/ */
/*   P4EST_FREE(ghost_data); */
/*   p4est_ghost_destroy (ghost); */

/*   return compare_dot; */
/* } */



void
serial_matrix_sym_tester
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs, /* only needed for # of nodes */
 d4est_elliptic_eqns_t* fcns,
 double sym_eps,
 d4est_operators_t* d4est_ops,
 int print,
 d4est_geometry_t* geom
)
{
  int i, local_nodes;
  local_nodes = vecs->local_nodes;
  
  double* a_mat = P4EST_ALLOC(double, local_nodes*local_nodes);
  double* a_mat_trans = P4EST_ALLOC(double, local_nodes*local_nodes);

  void* ghost_data;
  p4est_ghost_t *ghost;  
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);

  if (geom != NULL)
    ghost_data = (void*)P4EST_ALLOC (d4est_element_data_t, ghost->ghosts.elem_count);
  else
    ghost_data = (void*)P4EST_ALLOC (element_data_t, ghost->ghosts.elem_count);

  double* u_temp = P4EST_ALLOC(double, local_nodes);
  double* Au_temp = P4EST_ALLOC(double, local_nodes);
  double* tmp = vecs->u;
  double* tmp1 = vecs->Au;
  vecs->u = u_temp;
  vecs->Au = Au_temp;
  
  d4est_linalg_fill_vec(vecs->u, 0., vecs->local_nodes);
  
  for (i = 0; i < vecs->local_nodes; i++){
    vecs->u[i] = 1.;
    /* if (curved) */
      /* ((curved_d4est_elliptic_eqns_t*)fcns)->apply_lhs(p4est, ghost, (d4est_element_data_t*)ghost_data, vecs, d4est_ops, geom); */
    /* else */
    fcns->apply_lhs(p4est, ghost, ghost_data, vecs, d4est_ops, geom);

    /* to make parallel send stride as well */
    d4est_linalg_set_column(a_mat, vecs->Au, i, vecs->local_nodes, vecs->local_nodes);

    /* util_print_matrix(vecs->Au, vecs->local_nodes, 1, "Au = ", 0); */
    vecs->u[i] = 0.;
  }

  d4est_linalg_mat_transpose(a_mat, a_mat_trans, local_nodes);


  /* printf("sym_eps = %.25f\n", sym_eps); */
  /* for (int i = 0; i < local_nodes*local_nodes; i++){ */
    /* if (fabs(a_mat[i] - a_mat_trans[i]) > sym_eps){ */
    /*   printf("a_mat[%d] = %.25f\n", i, a_mat[i]); */
    /*   printf("a_mat_trans[%d] = %.25f\n", i, a_mat_trans[i]); */
    /*   printf("fabs(a_mat[i] - a_mat_trans[i]) = %.25f\n", fabs(a_mat[i] - a_mat_trans[i])); */
    /* } */
  /* } */
  int compare = util_compare_vecs(a_mat, a_mat_trans, local_nodes*local_nodes, sym_eps);

  double biggest_err = -1.;
  int biggest_id = -1;

  util_find_biggest_error(a_mat, a_mat_trans, local_nodes*local_nodes, &biggest_err, &biggest_id);
  
  if (compare == 1){
    printf("A_{DG} is symmetric\n");
    printf("biggest error is at %d with size %.20f\n", biggest_id, biggest_err);
  }
  else{
    printf("A_{DG} is NOT symmetric\n");
    printf("biggest error is at %d with size %.20f\n", biggest_id, biggest_err);
    /* print = 2; */
  }
  if (print == 1)
    util_print_matrix(a_mat, local_nodes, local_nodes, "A_{DG} = ",0);

  if (print == 2){
    printf("Printing vectors of %d length\n", local_nodes*local_nodes);
    /* util_print_matrices(a_mat, a_mat_trans, local_nodes*local_nodes, 1, "A, A^T = "); */
    for (int i = 0; i < local_nodes; i++){
      for (int j = 0; j < local_nodes; j++)
        printf("%d,%d A A^t: %.20f %.20f\n", i,j, a_mat[i*local_nodes + j], a_mat_trans[i*local_nodes + j]);
    }
  }

  if (print == 3){
    /* D4EST_ASSERT(curved); */
    D4EST_ASSERT(geom != NULL);
    printf("(i,j) Pairs that aren't equal\n");
    for (int i = 0; i < local_nodes; i++){
      for (int j = 0; j <= i; j++){
        if (fabs(a_mat[i*local_nodes + j] - a_mat_trans[i*local_nodes + j]) > sym_eps) {
          int e1 = d4est_element_data_debug_find_node(p4est, i);
          int e2 = d4est_element_data_debug_find_node(p4est, j);
          /* d4est_element_data_t* e1_data = d4est_element_data_get_element_data(p4est, e1); */
          /* d4est_element_data_t* e2_data = curved_elemenet_data_get_element_data(p4est, e2); */
          printf("node %d in element %d, node %d in element %d, A/A^t/ERR: %.20f %.20f %.20f\n", i, e1, j, e2, a_mat[i*local_nodes + j], a_mat_trans[i*local_nodes + j], fabs(a_mat[i*local_nodes + j] - a_mat_trans[i*local_nodes + j]));
        }
      }
    }
  }

  if (print == 4){
    D4EST_ASSERT(geom != NULL);
    printf("(i,j) Pairs that aren't equal\n");
    for (int i = 0; i < local_nodes; i++){
      for (int j = 0; j <= i; j++){
        if (fabs(a_mat[i*local_nodes + j] - a_mat_trans[i*local_nodes + j]) > sym_eps) {
          int e1 = d4est_element_data_debug_find_node(p4est, i);
          int e2 = d4est_element_data_debug_find_node(p4est, j);
          d4est_element_data_t* e1_data = d4est_element_data_get_element_data(p4est, e1);
          d4est_element_data_t* e2_data = d4est_element_data_get_element_data(p4est, e2);
          printf("node %d in element %d (x,y,z = %f,%f,%f), node %d in element %d (x,y,z = %f,%f,%f), A/A^t/ERR: %.20f %.20f %.20f\n", i, e1, e1_data->xyz[0][0],e1_data->xyz[1][0],  0., j, e2, e2_data->xyz[0][0],e2_data->xyz[1][0], 0., a_mat[i*local_nodes + j], a_mat_trans[i*local_nodes + j], fabs(a_mat[i*local_nodes + j] - a_mat_trans[i*local_nodes + j]));
        }
      }
    }
  }
  

  
  
  vecs->u = tmp;
  vecs->Au = tmp1;
  P4EST_FREE(u_temp);
  P4EST_FREE(Au_temp);
  P4EST_FREE(a_mat);
  P4EST_FREE(a_mat_trans);
  P4EST_FREE(ghost_data);
  p4est_ghost_destroy (ghost);
}



/* void */
/* parallel_matrix_sym_tester */
/* ( */
/*  p4est_t* p4est, */
/*  d4est_elliptic_data_t* vecs, /\* only needed for # of nodes *\/ */
/*  void* fcns, */
/*  double sym_eps, */
/*  d4est_operators_t* d4est_ops, */
/*  int curved, */
/*  int print */
/* ) */
/* { */
/*   int local_nodes; */
/*   local_nodes = vecs->local_nodes; */
  
/*   void* ghost_data; */
/*   p4est_ghost_t *ghost;   */
/*   ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE); */

/*   if (curved) */
/*     ghost_data = (void*)P4EST_ALLOC (d4est_element_data_t, ghost->ghosts.elem_count); */
/*   else */
/*     ghost_data = (void*)P4EST_ALLOC (element_data_t, ghost->ghosts.elem_count); */

/*   double* u_temp = P4EST_ALLOC(double, local_nodes); */
/*   double* Au_temp = P4EST_ALLOC(double, local_nodes); */
/*   double* tmp = vecs->u; */
/*   double* tmp1 = vecs->Au; */
/*   vecs->u = u_temp; */
/*   vecs->Au = Au_temp; */

  
/*   /\* seed 1467819770 produces ghosts with 2 elements*\/   */

/*   d4est_linalg_fill_vec(vecs->u, 0., vecs->local_nodes); */

/*   int* local_nodes_array = P4EST_ALLOC_ZERO(int, p4est->mpisize); */
/*   int* local_strides_array = P4EST_ALLOC_ZERO(int, p4est->mpisize + 1); */

/*   sc_allgather */
/*     ( */
/*      &local_nodes, */
/*      1, */
/*      sc_MPI_INT, */
/*      local_nodes_array, */
/*      1, */
/*      sc_MPI_INT, */
/*      sc_MPI_COMM_WORLD */
/*     ); */

/*   int stride = 0; */
/*   int global_nodes = 0; */
/*   for(int i = 0; i < p4est->mpisize; i++){ */
/*     global_nodes += local_nodes_array[i]; */
/*     printf("rank %d: local_nodes_array[%d] = %d\n", p4est->mpirank, i, local_nodes_array[i]); */
/*   } */
/*   printf("global_nodes = %d\n", global_nodes); */
/*   for (int i = 0; i < p4est->mpisize + 1; i++){ */
/*     local_strides_array[i] = stride; */
/*     printf("stride = %d\n", stride); */
/*     if (i < p4est->mpisize) */
/*       stride += local_nodes_array[i]; */
/*   } */

/*   double* a_mat = P4EST_ALLOC(double, local_nodes*global_nodes); */
/*   int  rstart = local_strides_array[p4est->mpirank]; */
/*   int rend = local_strides_array[p4est->mpirank + 1]; */

/*   printf("rank %d: [rstart,rend] = [%d,%d]\n", p4est->mpirank, rstart, rend);  */

/*   for (int c = 0; c < global_nodes; c++){ */

/*     if (c >= rstart && c < rend) */
/*       vecs->u[c - rstart] = 1.; */
    
/*     if (curved) */
/*       ((curved_d4est_elliptic_eqns_t*)fcns)->apply_lhs(p4est, ghost, (d4est_element_data_t*)ghost_data, vecs, d4est_ops, geom); */
/*     else */
/*       ((d4est_elliptic_eqns_t*)fcns)->apply_lhs(p4est, ghost, (element_data_t*)ghost_data, vecs, d4est_ops); */

/*     /\* to make parallel send stride as well *\/ */
/*     d4est_linalg_set_column_opt(a_mat, vecs->Au, c, vecs->local_nodes, global_nodes); */
    
/*     if (c >= rstart && c < rend) */
/*       vecs->u[c - rstart] = 0.; */
/*   } */

/*   double* amat_global = NULL; */
/*   double* amat_trans_global = NULL; */
/*   if (p4est->mpirank == 0){ */
/*     amat_global = P4EST_ALLOC(double, global_nodes*global_nodes); */
/*     amat_trans_global = P4EST_ALLOC(double, global_nodes*global_nodes); */
/*   } */
/*   MPI_Gather( */
/*     a_mat, */
/*     local_nodes*global_nodes, */
/*     sc_MPI_DOUBLE, */
/*     amat_global, */
/*     local_nodes*global_nodes, */
/*     sc_MPI_DOUBLE, */
/*     0, */
/*     sc_MPI_COMM_WORLD); */


/*   if (p4est->mpirank == 0){ */
/*     d4est_linalg_mat_transpose(amat_global, amat_trans_global, global_nodes); */
/*     int compare = util_compare_vecs(amat_global, amat_trans_global, global_nodes*global_nodes, sym_eps); */

/*     double biggest_err; */
/*     int biggest_id; */
/*   util_find_biggest_error(amat_global, amat_trans_global, global_nodes*global_nodes, &biggest_err, &biggest_id); */
  
/*   if (compare == 1){ */
/*     printf("A_{DG} is symmetric\n"); */
/*     printf("biggest error is at %d with size %.20f\n", biggest_id, biggest_err); */
/*   } */
/*   else{ */
/*     printf("A_{DG} is NOT symmetric\n"); */
/*     printf("biggest error is at %d with size %.20f\n", biggest_id, biggest_err); */
/*     /\* print = 2; *\/ */
/*   } */
  
/*     if (print == 1) */
/*       util_print_matrix(amat_global, global_nodes, global_nodes, "A_{DG} = ",0); */

/*     if (print == 2){ */
/*       printf("Printing vectors of %d length\n", global_nodes*global_nodes); */
/*       /\* util_print_matrices(a_mat, a_mat_trans, global_nodes*global_nodes, 1, "A, A^T = "); *\/ */
/*       for (int i = 0; i < global_nodes*global_nodes; i++){ */
/*         printf("A A^t [%d]: %.20f %.20f\n", i, amat_global[i], amat_trans_global[i]); */
/*       } */
/*     } */
/*   } */
  
/*   P4EST_FREE(amat_global); */
/*   P4EST_FREE(amat_trans_global); */
/*   vecs->u = tmp; */
/*   vecs->Au = tmp1; */
/*   P4EST_FREE(local_nodes_array);    */
/*   P4EST_FREE(local_strides_array);    */
/*   P4EST_FREE(u_temp); */
/*   P4EST_FREE(Au_temp); */
/*   P4EST_FREE(a_mat); */
/*   P4EST_FREE(ghost_data); */
/*   p4est_ghost_destroy (ghost); */
/* } */

int
matrix_sym_tester_parallel
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs, 
 void* fcns,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 int print,
 int num_tests,
 double sym_eps
)
{
  int tests = 0;
  int tests_passed = 0;
  int tests_failed = 0;
  int local_nodes = vecs->local_nodes;
  
  for (int i = 0; i < local_nodes; i++){
    if (tests > num_tests)
      break;
    for (int j = i; j < local_nodes; j++){
      tests++;
      double Aji, Aij;
      matrix_sym_tester_parallel_aux
        (
         p4est,
         vecs, 
         fcns,
         ghost,
         ghost_data,
         d4est_ops,
         d4est_geom,
         i,
         0,
         j,
         0,
         &Aji,
         &Aij
        );

      if (fabs(Aji - Aij) < sym_eps){
        if (print && p4est->mpirank == 0){
          printf("i = %d, j = %d, Aji = %.25f, Aij = %.25f = SYMMETRIC\n", i, j, Aji, Aij);
        }
        tests_passed++;
      }
      else {
        if (print && p4est->mpirank == 0){
          printf("i = %d, j = %d, Aji = %.25f, Aij = %.25f = NOT SYMMETRIC\n", i, j, Aji, Aij);
        }
        tests_failed++;
      }
      if (tests > num_tests)
        break;
    }

  }

  if(print && p4est->mpirank == 0){
    printf("Passed %d and Failed %d out of %d Tests\n", tests_passed, tests_failed, num_tests);
  }

  return (tests_failed == 0);
}

/* ghetto test of spd, just tests if diagonal is >0, which only means something if its <0 */
int
matrix_spd_tester_parallel
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs, 
 void* fcns,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 int print,
 int num_tests
)
{
  int tests = 0;
  int tests_passed = 0;
  int tests_failed = 0;
  int local_nodes = vecs->local_nodes;
  
  for (int i = 0; i < local_nodes; i++){
    if (tests > num_tests)
      break;
    int j = i;
    tests++;
      double Aji, Aij;
      matrix_sym_tester_parallel_aux
        (
         p4est,
         vecs, 
         fcns,
         ghost,
         ghost_data,
         d4est_ops,
         d4est_geom,
         i,
         0,
         j,
         0,
         &Aji,
         &Aij
        );

      if (Aij > 0){
        if (print && p4est->mpirank == 0){
          printf("i = %d, j = %d, Aji = %.25f, Aij = %.25f = PD\n", i, j, Aji, Aij);
        }
        tests_passed++;
      }
      else {
        if (print && p4est->mpirank == 0){
          printf("i = %d, j = %d, Aji = %.25f, Aij = %.25f = NOT PD\n", i, j, Aji, Aij);
        }
        tests_failed++;
      }
  }

  if(print && p4est->mpirank == 0){
    printf("Passed %d and Failed %d out of %d Tests\n", tests_passed, tests_failed, num_tests);
  }

  return (tests_failed == 0);
}

void
matrix_sym_tester_parallel_aux
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs, 
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
)
{
  d4est_elliptic_data_t vecs_for_sym_test;
  d4est_elliptic_data_copy_ptrs(vecs, &vecs_for_sym_test);
  
  double uj_A_ui_local; //Aji
  double ui_A_uj_local; //Aij

  double* ui = P4EST_ALLOC_ZERO(double, vecs->local_nodes);
  double* ui0 = P4EST_ALLOC_ZERO(double, vecs->local_nodes);
  double* uj = P4EST_ALLOC_ZERO(double, vecs->local_nodes);
  double* uj0 = P4EST_ALLOC_ZERO(double, vecs->local_nodes);
  double* A_ui = P4EST_ALLOC(double, vecs->local_nodes);
  double* A_uj = P4EST_ALLOC(double, vecs->local_nodes);

  
  if (p4est->mpirank == mpirank_i){
    D4EST_ASSERT(i < vecs->local_nodes);
    ui[i] = 1.;
  }
  if (p4est->mpirank == mpirank_j){
    D4EST_ASSERT(j < vecs->local_nodes);
    uj[j] = 1.;
  }
  
  vecs_for_sym_test.u = ui;
  vecs_for_sym_test.u0 = ui0;
  vecs_for_sym_test.Au = A_ui;

  /* if (d4est_geom != NULL){ */
    /* ((curved_d4est_elliptic_eqns_t*)fcns)->apply_lhs( */
                                              /* p4est, */
                                              /* ghost, */
                                              /* (d4est_element_data_t*)ghost_data, */
                                              /* &vecs_for_sym_test, */
                                              /* d4est_ops, */
                                              /* d4est_geom */
                                             /* ); */
  /* }   */
  /* else{ */
    ((d4est_elliptic_eqns_t*)fcns)->apply_lhs(
                                       p4est,
                                       ghost,
                                       ghost_data,
                                       &vecs_for_sym_test,
                                       d4est_ops,
                                       d4est_geom
                                      );
  /* } */
  uj_A_ui_local = d4est_linalg_vec_dot(uj, A_ui, vecs->local_nodes);

 
    
  vecs_for_sym_test.u = uj;
  vecs_for_sym_test.u0 = uj0;
  vecs_for_sym_test.Au = A_uj;
    
  /* if (d4est_geom != NULL) */
    /* ((curved_d4est_elliptic_eqns_t*)fcns)->apply_lhs(p4est, ghost, (d4est_element_data_t*)ghost_data, &vecs_for_sym_test, d4est_ops, d4est_geom); */
  /* else */
    ((d4est_elliptic_eqns_t*)fcns)->apply_lhs(p4est, ghost, ghost_data, &vecs_for_sym_test, d4est_ops, d4est_geom);
    
  ui_A_uj_local = d4est_linalg_vec_dot(ui, A_uj, vecs->local_nodes);

  double vAu_locals [2];
  double vAu_globals [2] = {0,0};
  vAu_locals[0] = uj_A_ui_local;
  vAu_locals[1] = ui_A_uj_local;

  sc_reduce
    (
     &vAu_locals,
     &vAu_globals,
     2,
     sc_MPI_DOUBLE,
     sc_MPI_SUM,
     0,
     sc_MPI_COMM_WORLD
    );

  *Aji = vAu_globals[0];
  *Aij = vAu_globals[1];

  P4EST_FREE(ui);
  P4EST_FREE(ui0);
  P4EST_FREE(uj);
  P4EST_FREE(uj0);
  P4EST_FREE(A_ui);
  P4EST_FREE(A_uj);
}
