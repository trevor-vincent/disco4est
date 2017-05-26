#include "d4est_operators.h"
#include "GL_and_GLL_nodes_and_weights.h"
#include <assert.h>
#include <math.h>
#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../pXest/pXest.h"

static const int d4est_operators_p8est_FToF_code[6][6] = {
  {0, 1, 1, 0, 0, 1}, {2, 0, 0, 1, 1, 0}, {2, 0, 0, 1, 1, 0},
  {0, 2, 2, 0, 0, 1}, {0, 2, 2, 0, 0, 1}, {2, 0, 0, 2, 2, 0}};

static const int d4est_operators_p8est_code_to_perm[3][4] = {
  {1, 2, 5, 6}, {0, 3, 4, 7}, {0, 4, 3, 7}};

static const int d4est_operators_p8est_perm_to_order[8][4] = {
  {0, 1, 2, 3}, {0, 2, 1, 3}, {1, 0, 3, 2}, {1, 3, 0, 2},
  {2, 0, 3, 1}, {2, 3, 0, 1}, {3, 1, 2, 0}, {3, 2, 1, 0}};
  

#define ORDERING 3

d4est_operators_t* d4est_ops_init() {
  d4est_operators_t* d4est_ops = P4EST_ALLOC(d4est_operators_t, 1);
  int deg_max = d4est_ops->d4est_operators_max_storage = (MAX_DEGREE);
  
  d4est_ops->d4est_operators_ref_Mij_1d_table = P4EST_ALLOC(double*, deg_max);
  d4est_ops->d4est_operators_ref_invMij_1d_table = P4EST_ALLOC(double*, deg_max);
  d4est_ops->d4est_operators_ref_invVij_1d_table = P4EST_ALLOC(double*, deg_max);
  d4est_ops->d4est_operators_ref_invGaussVij_1d_table = P4EST_ALLOC(double*, deg_max);

  d4est_ops->d4est_operators_GLL_nodes_1d_table = P4EST_ALLOC(double*, deg_max);
  d4est_ops->d4est_operators_GLL_weights_1d_table = P4EST_ALLOC(double*, deg_max);
  d4est_ops->d4est_operators_GL_nodes_1d_table = P4EST_ALLOC(double*, deg_max);
  d4est_ops->d4est_operators_GL_weights_1d_table = P4EST_ALLOC(double*, deg_max);
  
  d4est_ops->d4est_operators_ref_Dij_1d_table = P4EST_ALLOC(double*, deg_max);
  d4est_ops->d4est_operators_LIFT_1d_table = P4EST_ALLOC(double*, deg_max);
  d4est_ops->d4est_operators_ref_xyz_nd_table = P4EST_ALLOC(double*, deg_max*(2)); /* 2-> DIM = 2 or 3  */

  d4est_ops->d4est_operators_vtk_rst_2d_table = P4EST_ALLOC(double*, deg_max); /* 2-> DIM = 2 or 3  */
  d4est_ops->d4est_operators_vtk_rst_3d_table = P4EST_ALLOC(double*, deg_max); /* 2-> DIM = 2 or 3  */

  d4est_ops->d4est_operators_ref_Gauss_xyz_nd_table = P4EST_ALLOC(double*, deg_max*(2)); /* 2-> DIM = 2 or 3  */
  d4est_ops->d4est_operators_FLIP_1d_table = P4EST_ALLOC(double*, deg_max);

  d4est_ops->d4est_operators_hp_prolong_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  d4est_ops->d4est_operators_p_prolong_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  d4est_ops->d4est_operators_p_prolong_1d_inverse_table = P4EST_ALLOC(double*, deg_max * deg_max);
  d4est_ops->d4est_operators_hp_restrict_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  d4est_ops->d4est_operators_p_restrict_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  d4est_ops->d4est_operators_hp_prolong_transpose_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  d4est_ops->d4est_operators_p_prolong_transpose_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);  
  d4est_ops->d4est_operators_p_prolong_transpose_1d_inverse_table = P4EST_ALLOC(double*, deg_max * deg_max);  
  d4est_ops->d4est_operators_ref_GaussVij_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  d4est_ops->d4est_operators_ref_GaussDij_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_inverse_table = P4EST_ALLOC(double*, deg_max * deg_max);
  d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_inverse_table = P4EST_ALLOC(double*, deg_max * deg_max);
  d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  
  d4est_ops->d4est_operators_hp_restrict_interp_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  d4est_ops->d4est_operators_p_restrict_interp_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  d4est_ops->d4est_operators_vtk_interp_1d_table = P4EST_ALLOC(double*, deg_max);
  
  int i;
  for (i = 0; i < d4est_ops->d4est_operators_max_storage * d4est_ops->d4est_operators_max_storage; i++){
    d4est_ops->d4est_operators_hp_prolong_1d_table[i] = NULL;
    d4est_ops->d4est_operators_p_prolong_1d_table[i] = NULL;
    d4est_ops->d4est_operators_p_prolong_1d_inverse_table[i] = NULL;
    d4est_ops->d4est_operators_hp_restrict_1d_table[i] = NULL;
    d4est_ops->d4est_operators_hp_prolong_transpose_1d_table[i] = NULL;
    d4est_ops->d4est_operators_p_restrict_1d_table[i] = NULL;
    d4est_ops->d4est_operators_p_prolong_transpose_1d_table[i] = NULL;
    d4est_ops->d4est_operators_p_prolong_transpose_1d_inverse_table[i] = NULL;
    d4est_ops->d4est_operators_hp_restrict_interp_1d_table[i] = NULL;
    d4est_ops->d4est_operators_p_restrict_interp_1d_table[i] = NULL;
    d4est_ops->d4est_operators_ref_GaussVij_1d_table[i] = NULL;
    d4est_ops->d4est_operators_ref_GaussDij_1d_table[i] = NULL;
    d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_table[i] = NULL;
    d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_inverse_table[i] = NULL;
    d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_table[i] = NULL;
    d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_inverse_table[i] = NULL;
    if (i < d4est_ops->d4est_operators_max_storage){
      d4est_ops->d4est_operators_vtk_interp_1d_table[i] = NULL;
      d4est_ops->d4est_operators_ref_Mij_1d_table[i] = NULL;
      d4est_ops->d4est_operators_GLL_nodes_1d_table[i] = NULL;
      d4est_ops->d4est_operators_GLL_weights_1d_table[i] = NULL;
      d4est_ops->d4est_operators_GL_nodes_1d_table[i] = NULL;
      d4est_ops->d4est_operators_GL_weights_1d_table[i] = NULL;
      d4est_ops->d4est_operators_ref_invMij_1d_table[i] = NULL;
      d4est_ops->d4est_operators_vtk_rst_3d_table[i] = NULL;
      d4est_ops->d4est_operators_vtk_rst_2d_table[i] = NULL;
      d4est_ops->d4est_operators_ref_invVij_1d_table[i] = NULL;
      d4est_ops->d4est_operators_ref_invGaussVij_1d_table[i] = NULL;
      d4est_ops->d4est_operators_ref_Dij_1d_table[i] = NULL;
      d4est_ops->d4est_operators_LIFT_1d_table[i] = NULL;
      d4est_ops->d4est_operators_FLIP_1d_table[i] = NULL;
    }
    if (i < d4est_ops->d4est_operators_max_storage * (2)){
      d4est_ops->d4est_operators_ref_xyz_nd_table[i] = NULL;    
      d4est_ops->d4est_operators_ref_Gauss_xyz_nd_table[i] = NULL;
    }
  }

  return d4est_ops;
}

void d4est_ops_destroy(d4est_operators_t* d4est_ops) {

  int i;
  for (i = 0; i < d4est_ops->d4est_operators_max_storage; i++) {
    P4EST_FREE(d4est_ops->d4est_operators_ref_Mij_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_vtk_rst_3d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_vtk_rst_2d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_GLL_nodes_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_GLL_weights_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_GL_nodes_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_GL_weights_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_ref_invMij_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_ref_invVij_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_ref_invGaussVij_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_ref_Dij_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_LIFT_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_FLIP_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_vtk_interp_1d_table[i]);
  }
  P4EST_FREE(d4est_ops->d4est_operators_ref_Mij_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_vtk_rst_3d_table);
  P4EST_FREE(d4est_ops->d4est_operators_vtk_rst_2d_table);
  P4EST_FREE(d4est_ops->d4est_operators_vtk_interp_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_GLL_nodes_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_GLL_weights_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_GL_nodes_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_GL_weights_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_ref_invMij_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_ref_invVij_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_ref_invGaussVij_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_ref_Dij_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_LIFT_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_FLIP_1d_table);

  
  for (i = 0; i < d4est_ops->d4est_operators_max_storage * d4est_ops->d4est_operators_max_storage; i++) {
    P4EST_FREE(d4est_ops->d4est_operators_hp_prolong_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_hp_prolong_transpose_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_p_prolong_transpose_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_p_prolong_transpose_1d_inverse_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_p_prolong_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_p_prolong_1d_inverse_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_hp_restrict_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_hp_restrict_interp_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_p_restrict_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_p_restrict_interp_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_ref_GaussVij_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_ref_GaussDij_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_inverse_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_inverse_table[i]);
  }
  P4EST_FREE(d4est_ops->d4est_operators_hp_prolong_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_hp_prolong_transpose_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_p_prolong_transpose_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_p_prolong_transpose_1d_inverse_table);
  P4EST_FREE(d4est_ops->d4est_operators_hp_restrict_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_hp_restrict_interp_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_p_restrict_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_p_restrict_interp_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_ref_GaussVij_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_ref_GaussDij_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_inverse_table);
  P4EST_FREE(d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_inverse_table);
  P4EST_FREE(d4est_ops->d4est_operators_p_prolong_1d_table);
  P4EST_FREE(d4est_ops->d4est_operators_p_prolong_1d_inverse_table);
  
  for (i = 0; i < d4est_ops->d4est_operators_max_storage * (2); i++) {
    P4EST_FREE(d4est_ops->d4est_operators_ref_xyz_nd_table[i]);
    P4EST_FREE(d4est_ops->d4est_operators_ref_Gauss_xyz_nd_table[i]);
  }
  P4EST_FREE(d4est_ops->d4est_operators_ref_xyz_nd_table);
  P4EST_FREE(d4est_ops->d4est_operators_ref_Gauss_xyz_nd_table);
  P4EST_FREE(d4est_ops);
}

int d4est_operators_get_nodes(int dim, int deg) {
  mpi_assert(dim > 0 && dim < 4 && deg > 0);
  int i;
  int volume_nodes = 1;
  for (i = 0; i < dim; i++) volume_nodes *= deg + 1;
  return volume_nodes;
}

int d4est_operators_is_child_left_or_right(int c, int dir) {
  if (dir == 0)
    if (c == 0 || c == 2 || c == 4 || c == 6)
      return 0;
    else
      return 1;
  else if (dir == 1)
    if (c == 0 || c == 1 || c == 4 || c == 5)
      return 0;
    else
      return 1;
  else if (dir == 2)
    if (c == 0 || c == 1 || c == 2 || c == 3)
      return 0;
    else
      return 1;
  else {
    mpi_abort("Error: d4est_operators_is_child_left_or_right\n");
    return -1;
  }
  /* printf("Error: d4est_operators_is_child_left_or_right\n"); */
}

void d4est_operators_child_to_parent_ref_coords(int c, double* r) {
  *r *= .5;
  if (c == 0) {
    *r -= .5;
  } else if (c == 1) {
    *r += .5;
  } else {
    mpi_abort("Error: d4est_operators_child_to_parent_ref_coords\n");
    /* printf("ERROR: Boundary translation\n"); */
  }
}


void d4est_operators_build_Vij_1d(d4est_operators_t* d4est_ops, double* Vij_1d,
                                int deg) {
  int i, j, rows, cols;
  const double* lobatto_nodes = d4est_operators_fetch_GLL_nodes_1d(d4est_ops, deg);
  rows = cols = deg + 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      Vij_1d[i * cols + j] = d4est_operators_jacobi(lobatto_nodes[i], 0., 0., j);
}



static void d4est_operators_build_invVij_1d(d4est_operators_t* d4est_ops,
                                   double* invVij_1d, int deg)
{
  d4est_operators_build_Vij_1d(d4est_ops, invVij_1d, deg);
  linalg_invert(invVij_1d, deg + 1);
}


double* d4est_operators_fetch_invVij_1d(d4est_operators_t* d4est_ops, int deg) {
  mpi_assert(deg < d4est_ops->d4est_operators_max_storage);
  if (d4est_ops->d4est_operators_ref_invVij_1d_table[deg] != NULL) {
    return d4est_ops->d4est_operators_ref_invVij_1d_table[deg];
  }
  else {
    int size = (deg + 1) * (deg + 1);
    d4est_ops->d4est_operators_ref_invVij_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_ref_invVij_1d_table[deg];
    d4est_operators_build_invVij_1d(d4est_ops, op, deg);
    return op;
  }
}

void d4est_operators_hp_apply_nd_prolong_with_ptr(double* Uh, int degh, double* UH,
                                         int degH, int dim, int c,
                                         double* hp_prolong_matrix_1d) {
  mpi_assert((degH <= degh) && (c < (1 << dim)));
  int nodesh = degh + 1;
  int nodesH = degH + 1;

  if (dim == 1) {
    linalg_matvec_plus_vec(1.0, &hp_prolong_matrix_1d[c * nodesh * nodesH], UH,
                           0., Uh, nodesh, nodesH);
  } else if (dim == 2) {
    int cx = d4est_operators_is_child_left_or_right(c, 0);
    int cy = d4est_operators_is_child_left_or_right(c, 1);
#if (ORDERING==1)
       linalg_kron_A1A2x_NONSQR(Uh, &hp_prolong_matrix_1d[cx * nodesh * nodesH],
                          &hp_prolong_matrix_1d[cy * nodesh * nodesH], UH,
                          nodesh, nodesH, nodesh, nodesH);
    /* mpi_abort("Ordering = 1"); */
#endif

#if (ORDERING==3)
       linalg_kron_A1A2x_NONSQR(Uh, &hp_prolong_matrix_1d[cy * nodesh * nodesH],
                           &hp_prolong_matrix_1d[cx * nodesh * nodesH], UH,
                           nodesh, nodesH, nodesh, nodesH);
#endif
    
  } else if (dim == 3) {
    int cx = d4est_operators_is_child_left_or_right(c, 0);
    int cy = d4est_operators_is_child_left_or_right(c, 1);
    int cz = d4est_operators_is_child_left_or_right(c, 2);
#if (ORDERING==1)
    linalg_kron_A1A2A3x_NONSQR(Uh, &hp_prolong_matrix_1d[cx * nodesh * nodesH],
                               &hp_prolong_matrix_1d[cy * nodesh * nodesH],
                               &hp_prolong_matrix_1d[cz * nodesh * nodesH], UH,
                               nodesh, nodesH, nodesh, nodesH, nodesh, nodesH);
    /* mpi_abort("Ordering = 1"); */
#endif
#if (ORDERING==3)
    linalg_kron_A1A2A3x_NONSQR(Uh,
                               &hp_prolong_matrix_1d[cz * nodesh * nodesH],
                               &hp_prolong_matrix_1d[cy * nodesh * nodesH],
                               &hp_prolong_matrix_1d[cx * nodesh * nodesH], UH,
                               nodesh, nodesH, nodesh, nodesH, nodesh, nodesH);
#endif
  } else {
    mpi_abort("ERROR: d4est_operators_prolong_nd_U\n");
  }
}


void d4est_operators_build_ref_GaussVij_1d
(
 d4est_operators_t* d4est_ops,
 double* ref_Vij_Gauss_1d,
 int Lobatto_degree,
 int Gauss_degree
)
{
  double* Gauss_nodes = d4est_operators_fetch_GL_nodes_1d(d4est_ops, Gauss_degree);
  
  int rows = Gauss_degree + 1;
  int cols = Lobatto_degree + 1;
  
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      ref_Vij_Gauss_1d[i * cols + j] = d4est_operators_jacobi(Gauss_nodes[i], 0., 0., j);

  // util_print_matrix(Gauss_nodes, Gauss_degree + 1, 1, "Gauss_nodes = ", 0);
  
}

double* d4est_operators_fetch_ref_GaussVij_1d(d4est_operators_t* d4est_ops, int deg_Lobatto,
                                         int deg_Gauss) {
  mpi_assert(deg_Gauss <= d4est_ops->d4est_operators_max_storage);
  mpi_assert(deg_Lobatto <= d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_ref_GaussVij_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] != NULL) {
    return d4est_ops->d4est_operators_ref_GaussVij_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss];
  }
  else {
    int size = (deg_Gauss + 1) * (deg_Lobatto + 1);
    d4est_ops->d4est_operators_ref_GaussVij_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_ref_GaussVij_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss];
    d4est_operators_build_ref_GaussVij_1d(d4est_ops, op, deg_Lobatto, deg_Gauss);
    return op;
  }
}

void d4est_operators_build_ref_invGaussVij_1d
(
 d4est_operators_t* d4est_ops,
 double* ref_invGaussVij_1d,
 int Gauss_degree
)
{
  double* GaussVij = d4est_operators_fetch_ref_GaussVij_1d(d4est_ops, Gauss_degree, Gauss_degree);
  linalg_invert_and_copy(GaussVij, ref_invGaussVij_1d, Gauss_degree + 1);
}


double* d4est_operators_fetch_invGaussVij_1d(d4est_operators_t* d4est_ops, int deg) {
  mpi_assert(deg < d4est_ops->d4est_operators_max_storage);
  if (d4est_ops->d4est_operators_ref_invGaussVij_1d_table[deg] != NULL) {
    return d4est_ops->d4est_operators_ref_invGaussVij_1d_table[deg];
  }
  else {
    int size = (deg + 1) * (deg + 1);
    d4est_ops->d4est_operators_ref_invGaussVij_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_ref_invGaussVij_1d_table[deg];
    d4est_operators_build_ref_invGaussVij_1d(d4est_ops, op, deg);
    return op;
  }
}

void d4est_operators_build_ref_GLL_to_GL_interp_1d
(
 d4est_operators_t* d4est_ops,
 double* ref_GLL_to_GL_interp_1d,
 int Lobatto_degree,
 int Gauss_degree
)
{
  double* Gauss_nodes = d4est_operators_fetch_GL_nodes_1d(d4est_ops, Gauss_degree);
  double* ref_GaussVij = P4EST_ALLOC(double, (Gauss_degree + 1)*(Lobatto_degree+1));
  
  int rows = Gauss_degree + 1;
  int cols = Lobatto_degree + 1;
  
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      ref_GaussVij[i * cols + j] = d4est_operators_jacobi(Gauss_nodes[i], 0., 0., j);

  double* invVij = d4est_operators_fetch_invVij_1d(d4est_ops, Lobatto_degree);
  linalg_mat_multiply(ref_GaussVij, invVij, ref_GLL_to_GL_interp_1d, Gauss_degree + 1, Lobatto_degree + 1, Lobatto_degree + 1);

  P4EST_FREE(ref_GaussVij);
}

void d4est_operators_build_custom_GL_interp_1d
(
 d4est_operators_t* d4est_ops,
 double* custom_GL_interp_1d,
 int Lobatto_degree,
 int custom_degree,
 double* custom_points /* size = custom_degree + 1 */
)
{
  double* customVij = P4EST_ALLOC(double, (custom_degree + 1)*(Lobatto_degree+1));
  
  int rows = custom_degree + 1;
  int cols = Lobatto_degree + 1;
  
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      customVij[i * cols + j] = d4est_operators_jacobi(custom_points[i], 0., 0., j);

  double* invVij = d4est_operators_fetch_invVij_1d(d4est_ops, Lobatto_degree);
  linalg_mat_multiply(customVij, invVij, custom_GL_interp_1d, custom_degree + 1, Lobatto_degree + 1, Lobatto_degree + 1);

  P4EST_FREE(customVij);
}


/* void d4est_operators_build_ref_GL_to_GLL_interp_1d */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* ref_GL_to_GLL_interp_1d, */
/*  int Lobatto_degree, */
/*  int Gauss_degree */
/* ) */
/* { */
/*   double* Lobatto_nodes = d4est_operators_fetch_GLL_nodes_1d(d4est_ops, Gauss_degree); */
/*   double* ref_LobattoVij = P4EST_ALLOC(double, (Gauss_degree + 1)*(Lobatto_degree+1)); */
  
/*   int cols = Gauss_degree + 1; */
/*   int rows = Lobatto_degree + 1; */
  
/*   for (int i = 0; i < rows; i++) */
/*     for (int j = 0; j < cols; j++) */
/*       ref_LobattoVij[i * cols + j] = d4est_operators_jacobi(Lobatto_nodes[i], 0., 0., j); */

/*   double* invGaussVij = d4est_operators_fetch_invGaussVij_1d(d4est_ops, Gauss_degree); */
/*   linalg_mat_multiply(ref_LobattoVij, invGaussVij, ref_GL_to_GLL_interp_1d, Lobatto_degree + 1, Gauss_degree + 1, Gauss_degree + 1); */

/*   P4EST_FREE(ref_LobattoVij); */
/* } */



double*
d4est_operators_fetch_ref_GLL_to_GL_interp_1d(d4est_operators_t* d4est_ops, int deg_Lobatto,
                                         int deg_Gauss) {
  mpi_assert(deg_Gauss <= d4est_ops->d4est_operators_max_storage);
  mpi_assert(deg_Lobatto <= d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] != NULL) {
    return d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss];
  }
  else {
    int size = (deg_Gauss + 1) * (deg_Lobatto + 1);
    d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss];
    d4est_operators_build_ref_GLL_to_GL_interp_1d(d4est_ops, op, deg_Lobatto, deg_Gauss);
    return op;
  }
}

double* d4est_operators_fetch_ref_GLL_to_GL_interp_1d_inverse(d4est_operators_t* d4est_ops, int deg_Lobatto,
                                         int deg_Gauss) {
  mpi_assert(deg_Gauss <= d4est_ops->d4est_operators_max_storage);
  mpi_assert(deg_Lobatto <= d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_inverse_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] != NULL) {
    return d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_inverse_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss];
  }
  else {
    int size = (deg_Gauss + 1) * (deg_Lobatto + 1);
    d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_inverse_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_ref_GLL_to_GL_interp_1d_inverse_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss];
    d4est_operators_build_ref_GLL_to_GL_interp_1d_inverse(d4est_ops, op, deg_Lobatto, deg_Gauss);
    return op;
  }
}


void d4est_operators_build_ref_GLL_to_GL_interp_1d_inverse
(
 d4est_operators_t* d4est_ops,
 double* ref_GLL_to_GL_interp_1d_inverse,
 int Lobatto_degree,
 int Gauss_degree
)
{
  double* ref_GLL_to_GL_interp_1d = d4est_operators_fetch_ref_GLL_to_GL_interp_1d(d4est_ops, Lobatto_degree, Gauss_degree);
  linalg_leftinverse(ref_GLL_to_GL_interp_1d, ref_GLL_to_GL_interp_1d_inverse, Gauss_degree + 1, Lobatto_degree + 1);  
}

/* double* d4est_operators_fetch_ref_GL_to_GLL_interp_1d(d4est_operators_t* d4est_ops, int deg_Lobatto, */
/*                                          int deg_Gauss) { */
/*   mpi_assert(deg_Gauss <= d4est_ops->d4est_operators_max_storage); */
/*   mpi_assert(deg_Lobatto <= d4est_ops->d4est_operators_max_storage); */

/*   if (d4est_ops->d4est_operators_ref_GL_to_GLL_interp_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] != NULL) { */
/*     return d4est_ops->d4est_operators_ref_GL_to_GLL_interp_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss]; */
/*   } */
/*   else { */
/*     int size = (deg_Gauss + 1) * (deg_Lobatto + 1); */
/*     d4est_ops->d4est_operators_ref_GL_to_GLL_interp_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size); */
/*     double* op = d4est_ops->d4est_operators_ref_GL_to_GLL_interp_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss]; */
/*     d4est_operators_build_ref_GL_to_GLL_interp_1d(d4est_ops, op, deg_Lobatto, deg_Gauss); */
/*     return op; */
/*   } */
/* } */




/* void d4est_operators_build_ref_GL_to_GLL_interp_trans_1d */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* ref_GL_to_GLL_interp_trans_1d, */
/*  int Lobatto_degree, */
/*  int Gauss_degree */
/* ) */
/* { */
/*   double* ref_GL_to_GLL_interp_1d = d4est_operators_fetch_ref_GL_to_GLL_interp_1d(d4est_ops, Lobatto_degree, Gauss_degree); */
/*   linalg_mat_transpose_nonsqr(ref_GL_to_GLL_interp_1d, ref_GL_to_GLL_interp_trans_1d, Lobatto_degree + 1, Gauss_degree + 1); */
/* } */


void d4est_operators_build_ref_GLL_to_GL_interp_trans_1d
(
 d4est_operators_t* d4est_ops,
 double* ref_GLL_to_GL_interp_trans_1d,
 int Lobatto_degree,
 int Gauss_degree
)
{
  double* ref_GLL_to_GL_interp_1d = d4est_operators_fetch_ref_GLL_to_GL_interp_1d(d4est_ops, Lobatto_degree, Gauss_degree);
  linalg_mat_transpose_nonsqr(ref_GLL_to_GL_interp_1d, ref_GLL_to_GL_interp_trans_1d, Gauss_degree + 1, Lobatto_degree + 1);
}

double* d4est_operators_fetch_ref_GLL_to_GL_interp_trans_1d(d4est_operators_t* d4est_ops, int deg_Lobatto,
                                         int deg_Gauss) {
  mpi_assert(deg_Gauss <= d4est_ops->d4est_operators_max_storage);
  mpi_assert(deg_Lobatto <= d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] != NULL) {
    return d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss];
  }
  else {
    int size = (deg_Gauss + 1) * (deg_Lobatto + 1);
    d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss];
    d4est_operators_build_ref_GLL_to_GL_interp_trans_1d(d4est_ops, op, deg_Lobatto, deg_Gauss);
    return op;
  }
}


/* double* d4est_operators_fetch_ref_GLL_to_GL_interp_trans_1d_inverse(d4est_operators_t* d4est_ops, int deg_Lobatto, */
/*                                          int deg_Gauss) { */
/*   mpi_assert(deg_Gauss <= d4est_ops->d4est_operators_max_storage); */
/*   mpi_assert(deg_Lobatto <= d4est_ops->d4est_operators_max_storage); */

/*   if (d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_inverse_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] != NULL) { */
/*     return d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_inverse_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss]; */
/*   } */
/*   else { */
/*     int size = (deg_Gauss + 1) * (deg_Lobatto + 1); */
/*     d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_inverse_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size); */
/*     double* op = d4est_ops->d4est_operators_ref_GLL_to_GL_interp_trans_1d_inverse_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss]; */
/*     d4est_operators_build_ref_GLL_to_GL_interp_trans_1d_inverse(d4est_ops, op, deg_Lobatto, deg_Gauss); */
/*     return op; */
/*   } */
/* } */




/* void d4est_operators_build_ref_GLL_to_GL_interp_trans_1d_inverse */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* ref_GLL_to_GL_interp_trans_1d_inverse, */
/*  int Lobatto_degree, */
/*  int Gauss_degree */
/* ) */
/* { */
/*   double* ref_GLL_to_GL_interp_trans_1d = d4est_operators_fetch_ref_GLL_to_GL_interp_trans_1d(d4est_ops, Lobatto_degree, Gauss_degree); */
/*   linalg_leftinverse(ref_GLL_to_GL_interp_trans_1d, ref_GLL_to_GL_interp_trans_1d_inverse, Lobatto_degree + 1, Gauss_degree + 1); */
/* } */


/* /\* double* d4est_operators_fetch_ref_GL_to_GLL_interp_trans_1d(d4est_operators_t* d4est_ops, int deg_Lobatto, *\/ */
/* /\*                                          int deg_Gauss) { *\/ */
/* /\*   mpi_assert(deg_Gauss <= d4est_ops->d4est_operators_max_storage); *\/ */
/* /\*   mpi_assert(deg_Lobatto <= d4est_ops->d4est_operators_max_storage); *\/ */

/* /\*   if (d4est_ops->d4est_operators_ref_GL_to_GLL_interp_trans_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] != NULL) { *\/ */
/* /\*     return d4est_ops->d4est_operators_ref_GL_to_GLL_interp_trans_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss]; *\/ */
/* /\*   } *\/ */
/* /\*   else { *\/ */
/* /\*     int size = (deg_Gauss + 1) * (deg_Lobatto + 1); *\/ */
/* /\*     d4est_ops->d4est_operators_ref_GL_to_GLL_interp_trans_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size); *\/ */
/* /\*     double* op = d4est_ops->d4est_operators_ref_GL_to_GLL_interp_trans_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss]; *\/ */
/* /\*     d4est_operators_build_ref_GL_to_GLL_interp_trans_1d(d4est_ops, op, deg_Lobatto, deg_Gauss); *\/ */
/* /\*     return op; *\/ */
/* /\*   } *\/ */
/* /\* } *\/ */

/* void d4est_operators_interp_old */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* u_in, */
/*  quadrature_type_t quad_type_in, */
/*  int deg_in, */
/*  double* u_out, */
/*  quadrature_type_t quad_type_out, */
/*  int deg_out, */
/*  int dim */
/* ) */
/* { */
/*   mpi_assert(deg_in <= deg_out); */
/*   if (quad_type_in == QUAD_LOBATTO && quad_type_out == QUAD_LOBATTO){ */
/*     d4est_operators_apply_p_prolong(d4est_ops, u_in, deg_in, dim, deg_out, u_out); */
/*   } */
/*   else if (quad_type_in == QUAD_LOBATTO && quad_type_out == QUAD_GAUSS){ */
/*     d4est_operators_interp_Lobatto_to_Gauss(d4est_ops, u_in, deg_in, deg_out, u_out, dim); */
/*   } */
/*   else { */
/*     mpi_abort("[D4EST_ERROR]: Currently we only support interpolation from Lobatto to Lobatto or Lobatto to Gauss"); */
/*   } */
/* } */



/* void d4est_operators_interp_Lobatto_to_Gauss */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* u_Lobatto_in, */
/*  int deg_Lobatto, */
/*  int deg_Gauss, */
/*  double* u_Gauss_out, */
/*  int dim */
/* ) */
/* { */
/*   mpi_assert(dim == 1 || dim == 2 || dim == 3); */
/*   double* GLL_to_GL_interp = d4est_operators_fetch_ref_GLL_to_GL_interp_1d(d4est_ops, deg_Lobatto, deg_Gauss); */

/*   int nodes_Lobatto = deg_Lobatto + 1; */
/*   int nodes_Gauss = deg_Gauss + 1; */
  
/*   if (dim == 1){ */
/*     linalg_matvec_plus_vec(1.0, GLL_to_GL_interp, u_Lobatto_in, 0., u_Gauss_out, nodes_Gauss, nodes_Lobatto); */
/*   } */
/*   else if (dim == 2) { */
/*     linalg_kron_A1A2x_NONSQR(u_Gauss_out, GLL_to_GL_interp, GLL_to_GL_interp, u_Lobatto_in, nodes_Gauss, nodes_Lobatto, */
/*                              nodes_Gauss, nodes_Lobatto); */
/*   } */
/*   else if (dim == 3){ */
/*     linalg_kron_A1A2A3x_NONSQR(u_Gauss_out, GLL_to_GL_interp, GLL_to_GL_interp, GLL_to_GL_interp, u_Lobatto_in, */
/*                                nodes_Gauss, nodes_Lobatto, nodes_Gauss, nodes_Lobatto, nodes_Gauss, nodes_Lobatto); */
/*   } */
/*   else { */
/*     mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension."); */
/*   } */
/* } */



/* /\* void d4est_operators_interp_GL_to_GLL *\/ */
/* /\* ( *\/ */
/* /\*  d4est_operators_t* d4est_ops, *\/ */
/* /\*  double* u_Gauss_in, *\/ */
/* /\*  int deg_Gauss, *\/ */
/* /\*  int deg_Lobatto, *\/ */
/* /\*  double* u_Lobatto_out, *\/ */
/* /\*  int dim *\/ */
/* /\* ) *\/ */
/* /\* { *\/ */
/* /\*   mpi_assert(dim == 1 || dim == 2 || dim == 3); *\/ */
/* /\*   double* GL_to_GLL_interp = d4est_operators_fetch_ref_GL_to_GLL_interp_1d(d4est_ops, deg_Lobatto, deg_Gauss); *\/ */

/* /\*   int nodes_Gauss = deg_Gauss + 1; *\/ */
/* /\*   int nodes_Lobatto = deg_Lobatto + 1; *\/ */
  
/* /\*   if (dim == 1){ *\/ */
/* /\*     linalg_matvec_plus_vec(1.0, GL_to_GLL_interp, u_Gauss_in, 0., u_Lobatto_out, nodes_Lobatto, nodes_Gauss); *\/ */
/* /\*   } *\/ */
/* /\*   else if (dim == 2) { *\/ */
/* /\*     linalg_kron_A1A2x_NONSQR(u_Lobatto_out, GL_to_GLL_interp, GL_to_GLL_interp, u_Gauss_in, nodes_Lobatto, nodes_Gauss, *\/ */
/* /\*                              nodes_Lobatto, nodes_Gauss); *\/ */
/* /\*   } *\/ */
/* /\*   else if (dim == 3){ *\/ */
/* /\*     linalg_kron_A1A2A3x_NONSQR(u_Lobatto_out, GL_to_GLL_interp, GL_to_GLL_interp, GL_to_GLL_interp, u_Gauss_in, *\/ */
/* /\*                                nodes_Lobatto, nodes_Gauss, nodes_Lobatto, nodes_Gauss, nodes_Lobatto, nodes_Gauss); *\/ */
/* /\*   } *\/ */
/* /\*   else { *\/ */
/* /\*     mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension."); *\/ */
/* /\*   } *\/ */
/* /\* } *\/ */


/* void d4est_operators_apply_curved_mass_matrix_old */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* in, */
/*  int deg_lobatto, */
/*  double* jac_quad, */
/*  int deg_quad, */
/*  quadrature_type_t quad_type, */
/*  int dim, */
/*  double* out */
/* ){ */
/*   mpi_assert(dim == 1 || dim == 2 || dim == 3); */

/*   /\* for now assert this, can get rid of by p-prolonging then p-restricting *\/ */
/*   /\* mpi_assert(deg_lobatto == deg_quad); *\/ */
/*   int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad); */

/*   double* in_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*   double* w_j_in_quad = P4EST_ALLOC(double, volume_nodes_quad); */

/*   double* interp_trans = NULL; */
/*   double* quad_weights = NULL; */

/*   if (quad_type == QUAD_GAUSS){ */
/*     quad_weights = d4est_operators_fetch_GL_weights_1d(d4est_ops, deg_quad); */
/*     interp_trans = d4est_operators_fetch_ref_GLL_to_GL_interp_trans_1d */
/*                                    (d4est_ops, deg_lobatto, deg_quad);   */
/*   } */
/*   else if (quad_type == QUAD_LOBATTO){ */
/*     quad_weights = d4est_operators_fetch_GLL_weights_1d(d4est_ops, deg_quad); */
/*     interp_trans = d4est_operators_fetch_p_prolong_transpose_1d(d4est_ops, deg_lobatto, deg_quad); */
/*   } */
/*   else { */
/*     mpi_abort("[D4EST_ERROR]: Not a support quad type"); */
/*   } */
  
/*   d4est_operators_interp(d4est_ops, in, QUAD_LOBATTO, deg_lobatto, in_quad, quad_type, deg_quad, dim); */
  
/*   int nodes_lobatto = deg_lobatto + 1; */
/*   int nodes_quad = deg_quad + 1; */
  
/*   if (dim == 1){ */
/*     linalg_kron_vec_dot_xy(quad_weights, jac_quad, in_quad, nodes_quad, w_j_in_quad); */
/*     linalg_matvec_plus_vec(1.0, interp_trans, w_j_in_quad, 0., out, nodes_lobatto, nodes_quad); */
/*   } */
/*   else if (dim == 2) { */
/*     linalg_kron_vec_o_vec_dot_xy(quad_weights, jac_quad, in_quad, nodes_quad, w_j_in_quad); */
/*     linalg_kron_A1A2x_NONSQR(out, interp_trans, interp_trans, w_j_in_quad, nodes_lobatto, nodes_quad, */
/*                              nodes_lobatto, nodes_quad); */
/*   } else if (dim == 3) { */
/*     linalg_kron_vec_o_vec_o_vec_dot_xy(quad_weights, jac_quad, in_quad, nodes_quad, w_j_in_quad); */
/*     linalg_kron_A1A2A3x_NONSQR(out, interp_trans, interp_trans, interp_trans, w_j_in_quad, */
/*                                nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad); */
/*   } else {     */
/*     P4EST_FREE(w_j_in_quad); */
/*     P4EST_FREE(in_quad); */
/*     mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension."); */
/*   }   */
/*   P4EST_FREE(w_j_in_quad); */
/*   P4EST_FREE(in_quad); */
/* } */




/* void d4est_operators_compute_curved_inverse_stiffness_matrix */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  int deg_lobatto, */
/*  double* jac_quad, */
/*  double* rst_xyz [(P4EST_DIM)][(P4EST_DIM)], /\* must be padded with NULL if you want lower dim *\/ */
/*  int deg_quad, */
/*  quadrature_type_t quad_type, */
/*  int dim, */
/*  double* invS */
/* ) */
/* { */
/*   int volume_nodes_lobatto = d4est_operators_get_nodes(dim,deg_lobatto); */
  
/*   double* S = P4EST_ALLOC(double, volume_nodes_lobatto*volume_nodes_lobatto); */
/*   double* u = P4EST_ALLOC_ZERO(double, volume_nodes_lobatto); */
/*   double* Su = P4EST_ALLOC_ZERO(double, volume_nodes_lobatto); */

/*   for (int i = 0; i < volume_nodes_lobatto; i++){ */
/*     u[i] = 1.; */
/*     d4est_operators_apply_curved_stiffness_matrix(d4est_ops, u, deg_lobatto, jac_quad, rst_xyz, deg_quad, quad_type, dim, Su); */
/*     linalg_set_column(S, Su, i, volume_nodes_lobatto, volume_nodes_lobatto); */
/*     u[i] = 0.; */
/*   } */
/*   linalg_invert_and_copy(S, invS, volume_nodes_lobatto); */

/*   P4EST_FREE(Su); */
/*   P4EST_FREE(u); */
/*   P4EST_FREE(S); */
/* } */


/* void d4est_operators_apply_curved_stiffness_matrix_old */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* in, */
/*  int deg_lobatto, */
/*  double* jac_quad, */
/*  double* rst_xyz [(P4EST_DIM)][(P4EST_DIM)], /\* must be padded with NULL if you want lower dim *\/ */
/*  int deg_quad, */
/*  quadrature_type_t quad_type, */
/*  int dim, */
/*  double* out */
/* ){ */
/*   mpi_assert(dim == 1 || dim == 2 || dim == 3); */

/*   /\* for now assert this, can get rid of by p-prolonging then p-restricting *\/ */
/*   /\* mpi_assert(deg_lobatto == deg_quad); *\/ */
/*   int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad); */
/*   int volume_nodes_lobatto = d4est_operators_get_nodes(dim, deg_lobatto); */

/*   double* Dl_in = P4EST_ALLOC(double, volume_nodes_lobatto); */
/*   double* V_Dl_in = P4EST_ALLOC(double, volume_nodes_quad); */
/*   double* W_V_Dl_in = P4EST_ALLOC(double, volume_nodes_quad); */
/*   double* VT_W_V_Dl_in = P4EST_ALLOC(double, volume_nodes_lobatto); */
/*   double* DTlp_VT_W_V_Dl_in = P4EST_ALLOC(double, volume_nodes_lobatto); */
  
/*   double* interp_trans = NULL; */
/*   double* quad_weights = NULL; */

/*   if (quad_type == QUAD_GAUSS){ */
/*     quad_weights = d4est_operators_fetch_GL_weights_1d(d4est_ops, deg_quad); */
/*     interp_trans = d4est_operators_fetch_ref_GLL_to_GL_interp_trans_1d */
/*                                    (d4est_ops, deg_lobatto, deg_quad);   */
/*   } */
/*   else if (quad_type == QUAD_LOBATTO){ */
/*     quad_weights = d4est_operators_fetch_GLL_weights_1d(d4est_ops, deg_quad); */
/*     interp_trans = d4est_operators_fetch_p_prolong_transpose_1d(d4est_ops, deg_lobatto, deg_quad); */
/*   } */
/*   else { */
/*     mpi_abort("[D4EST_ERROR]: Not a support quad type"); */
/*   } */
  
  
/*   int nodes_lobatto = deg_lobatto + 1; */
/*   int nodes_quad = deg_quad + 1; */
/*   linalg_fill_vec(out, 0., volume_nodes_lobatto); */
  
/*   for (int k = 0; k < dim; k++){ */
/*     for (int lp = 0; lp < dim; lp++){ */
/*       for (int l = 0; l < dim; l++){ */
/*         /\* Apply Dbar in lth direction *\/ */
/*         d4est_operators_apply_Dij(d4est_ops, in, dim, deg_lobatto, l, Dl_in); */

/*         /\* Interpolate to quadrature points *\/ */
/*         d4est_operators_interp(d4est_ops, Dl_in, QUAD_LOBATTO, deg_lobatto, V_Dl_in, quad_type, deg_quad, dim); */
        
/*         if (dim == 3) { */
/*           linalg_kron_vec_o_vec_o_vec_dot_wxyz(quad_weights, jac_quad, rst_xyz[l][k], rst_xyz[lp][k], V_Dl_in, nodes_quad, W_V_Dl_in); */
/*           linalg_kron_A1A2A3x_NONSQR(VT_W_V_Dl_in, interp_trans, interp_trans, interp_trans, W_V_Dl_in, */
/*                                      nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad); */
/*         } */
/*         else if (dim == 2) { */
/*           linalg_kron_vec_o_vec_dot_wxyz(quad_weights, jac_quad, rst_xyz[l][k], rst_xyz[lp][k], V_Dl_in, nodes_quad, W_V_Dl_in); */
/*           linalg_kron_A1A2x_NONSQR(VT_W_V_Dl_in, interp_trans, interp_trans, W_V_Dl_in, nodes_lobatto, nodes_quad, */
/*                                    nodes_lobatto, nodes_quad); */
/*         } */
/*         else if (dim == 1){ */
/*           linalg_kron_vec_dot_wxyz(quad_weights, jac_quad, V_Dl_in, rst_xyz[l][k], rst_xyz[lp][k], nodes_quad, W_V_Dl_in); */
/*           linalg_matvec_plus_vec(1.0, interp_trans, W_V_Dl_in, 0., VT_W_V_Dl_in, nodes_lobatto, nodes_quad); */
/*         } */
          
/*         else {     */
/*           P4EST_FREE(DTlp_VT_W_V_Dl_in); */
/*           P4EST_FREE(VT_W_V_Dl_in); */
/*           P4EST_FREE(W_V_Dl_in); */
/*           P4EST_FREE(V_Dl_in); */
/*           P4EST_FREE(Dl_in); */
/*           mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension."); */
/*         } */

/*         d4est_operators_apply_Dij_transpose(d4est_ops, VT_W_V_Dl_in, dim, deg_lobatto, lp, DTlp_VT_W_V_Dl_in); */
/*         linalg_vec_axpy(1., DTlp_VT_W_V_Dl_in, out, volume_nodes_lobatto); */
/*       } */
/*     } */
/*   } */
  
/*   P4EST_FREE(DTlp_VT_W_V_Dl_in); */
/*   P4EST_FREE(VT_W_V_Dl_in); */
/*   P4EST_FREE(W_V_Dl_in); */
/*   P4EST_FREE(V_Dl_in); */
/*   P4EST_FREE(Dl_in); */
/* } */

/* void d4est_operators_apply_curved_stiffness_matrix */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* in, */
/*  int deg_lobatto, */
/*  double* jac_quad, */
/*  double* rst_xyz [(P4EST_DIM)][(P4EST_DIM)], /\* must be padded with NULL if you want lower dim *\/ */
/*  int deg_quad, */
/*  double** quad_weights, */
/*  double** interp_lobatto_to_quad, */
/*  double** interp_lobatto_to_quad_trans, */
/*  int dim, */
/*  double* out */
/* ){ */
/*   mpi_assert(dim == 1 || dim == 2 || dim == 3); */

/*   /\* for now assert this, can get rid of by p-prolonging then p-restricting *\/ */
/*   /\* mpi_assert(deg_lobatto == deg_quad); *\/ */
/*   int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad); */
/*   int volume_nodes_lobatto = d4est_operators_get_nodes(dim, deg_lobatto); */

/*   double* Dl_in = P4EST_ALLOC(double, volume_nodes_lobatto); */
/*   double* V_Dl_in = P4EST_ALLOC(double, volume_nodes_quad); */
/*   double* W_V_Dl_in = P4EST_ALLOC(double, volume_nodes_quad); */
/*   double* VT_W_V_Dl_in = P4EST_ALLOC(double, volume_nodes_lobatto); */
/*   double* DTlp_VT_W_V_Dl_in = P4EST_ALLOC(double, volume_nodes_lobatto); */
  
  
/*   int nodes_lobatto = deg_lobatto + 1; */
/*   int nodes_quad = deg_quad + 1; */
/*   linalg_fill_vec(out, 0., volume_nodes_lobatto); */
  
/*   for (int k = 0; k < dim; k++){ */
/*     for (int lp = 0; lp < dim; lp++){ */
/*       for (int l = 0; l < dim; l++){ */
/*         /\* Apply Dbar in lth direction *\/ */
/*         d4est_operators_apply_Dij(d4est_ops, in, dim, deg_lobatto, l, Dl_in); */

/*         /\* Interpolate to quadrature points *\/ */
/*         /\* d4est_operators_interp(d4est_ops, Dl_in, QUAD_LOBATTO, deg_lobatto, V_Dl_in, quad_type, deg_quad, dim); *\/ */
        
/*         if (dim == 3) { */

/*               linalg_kron_A1A2A3x_NONSQR(V_Dl_in, interp_lobatto_to_quad[2], interp_lobatto_to_quad[1], interp_lobatto_to_quad[0], Dl_in, */
/*                                nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto); */
              
/*               linalg_kron_vec1_o_vec2_o_vec3_dot_wxyz(quad_weights[2],quad_weights[1],quad_weights[0], jac_quad, rst_xyz[l][k], rst_xyz[lp][k], V_Dl_in, nodes_quad, W_V_Dl_in); */
/*           linalg_kron_A1A2A3x_NONSQR(VT_W_V_Dl_in, interp_lobatto_to_quad_trans[2], interp_lobatto_to_quad_trans[1], interp_lobatto_to_quad_trans[0], W_V_Dl_in, */
/*                                      nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad); */
/*         } */
/*         else if (dim == 2) { */
/*           linalg_kron_A1A2x_NONSQR(V_Dl_in, interp_lobatto_to_quad[1], interp_lobatto_to_quad[0], Dl_in, nodes_quad, nodes_lobatto, */
/*                             nodes_quad, nodes_lobatto); */
          
/*           linalg_kron_vec1_o_vec2_dot_wxyz(quad_weights[1], quad_weights[0], jac_quad, rst_xyz[l][k], rst_xyz[lp][k], V_Dl_in, nodes_quad, W_V_Dl_in); */
/*           linalg_kron_A1A2x_NONSQR(VT_W_V_Dl_in, interp_lobatto_to_quad_trans[1], interp_lobatto_to_quad_trans[0], W_V_Dl_in, nodes_lobatto, nodes_quad, */
/*                                    nodes_lobatto, nodes_quad); */
/*         } */
/*         else if (dim == 1){ */
/*           linalg_matvec_plus_vec(1.0, interp_lobatto_to_quad[0], Dl_in, 0.,  V_Dl_in, nodes_quad, nodes_lobatto); */
/*           linalg_kron_vec_dot_wxyz(quad_weights[0], jac_quad, V_Dl_in, rst_xyz[l][k], rst_xyz[lp][k], nodes_quad, W_V_Dl_in); */
/*           linalg_matvec_plus_vec(1.0, interp_lobatto_to_quad_trans[0], W_V_Dl_in, 0., VT_W_V_Dl_in, nodes_lobatto, nodes_quad); */
/*         } */
          
/*         else {     */
/*           P4EST_FREE(DTlp_VT_W_V_Dl_in); */
/*           P4EST_FREE(VT_W_V_Dl_in); */
/*           P4EST_FREE(W_V_Dl_in); */
/*           P4EST_FREE(V_Dl_in); */
/*           P4EST_FREE(Dl_in); */
/*           mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension."); */
/*         } */

/*         d4est_operators_apply_Dij_transpose(d4est_ops, VT_W_V_Dl_in, dim, deg_lobatto, lp, DTlp_VT_W_V_Dl_in); */
/*         linalg_vec_axpy(1., DTlp_VT_W_V_Dl_in, out, volume_nodes_lobatto); */
/*       } */
/*     } */
/*   } */
  
/*   P4EST_FREE(DTlp_VT_W_V_Dl_in); */
/*   P4EST_FREE(VT_W_V_Dl_in); */
/*   P4EST_FREE(W_V_Dl_in); */
/*   P4EST_FREE(V_Dl_in); */
/*   P4EST_FREE(Dl_in); */
/* } */


/* /\* void d4est_operators_apply_curvedGaussStiff_withJdrdxdrdx *\/ */
/* /\* ( *\/ */
/* /\*  d4est_operators_t* d4est_ops, *\/ */
/* /\*  double* in, *\/ */
/* /\*  int deg_Lobatto, *\/ */
/* /\*  double* Jdrdxdrdx [(P4EST_DIM)][(P4EST_DIM)], *\/ */
/* /\*  int deg_Gauss, *\/ */
/* /\*  int dim, *\/ */
/* /\*  double* out *\/ */
/* /\* ){ *\/ */
/* /\*   mpi_assert(dim == 1 || dim == 2 || dim == 3); *\/ */

/* /\*   /\\* for now assert this, can get rid of by p-prolonging then p-restricting *\\/ *\/ */
/* /\*   /\\* mpi_assert(deg_Lobatto == deg_Gauss); *\\/ *\/ */
/* /\*   int volume_nodes_Gauss = d4est_operators_get_nodes(dim, deg_Gauss); *\/ */
/* /\*   int volume_nodes_Lobatto = d4est_operators_get_nodes(dim, deg_Lobatto); *\/ */

/* /\*   double* Dl_in = P4EST_ALLOC(double, volume_nodes_Lobatto); *\/ */
/* /\*   double* V_Dl_in = P4EST_ALLOC(double, volume_nodes_Gauss); *\/ */
/* /\*   double* W_V_Dl_in = P4EST_ALLOC(double, volume_nodes_Gauss); *\/ */
/* /\*   double* VT_W_V_Dl_in = P4EST_ALLOC(double, volume_nodes_Lobatto); *\/ */
/* /\*   double* DTlp_VT_W_V_Dl_in = P4EST_ALLOC(double, volume_nodes_Lobatto); *\/ */
  
/* /\*   double* Gauss_weights = d4est_operators_fetch_GL_weights_1d(d4est_ops, deg_Gauss); *\/ */
/* /\*   double* GLL_to_GL_interp_trans = d4est_operators_fetch_ref_GLL_to_GL_interp_trans_1d *\/ */
/* /\*                                    (d4est_ops, deg_Lobatto, deg_Gauss); *\/ */

/* /\*   int nodes_Lobatto = deg_Lobatto + 1; *\/ */
/* /\*   int nodes_Gauss = deg_Gauss + 1; *\/ */
/* /\*   linalg_fill_vec(out, 0., volume_nodes_Lobatto); *\/ */
  
/* /\*   for (int k = 0; k < dim; k++){ *\/ */
/* /\*     for (int lp = 0; lp < dim; lp++){ *\/ */
/* /\*       for (int l = 0; l < dim; l++){ *\/ */
/* /\*         /\\* Apply Dbar in lth direction *\\/ *\/ */
/* /\*         d4est_operators_apply_Dij(d4est_ops, in, dim, deg_Lobatto, l, Dl_in); *\/ */

/* /\*           /\\* Apply Gauss quadrature vandermonde*\\/ *\/ */
/* /\*         d4est_operators_interp_GLL_to_GL(d4est_ops, Dl_in, deg_Lobatto, deg_Gauss, V_Dl_in, dim); *\/ */

/* /\*         if (dim == 3) { *\/ */
/* /\*           linalg_kron_vec_o_vec_o_vec_dot_xy(Gauss_weights, Jdrdxdrdx[l][lp], V_Dl_in, nodes_Gauss, W_V_Dl_in); *\/ */
/* /\*           linalg_kron_A1A2A3x_NONSQR(VT_W_V_Dl_in, GLL_to_GL_interp_trans, GLL_to_GL_interp_trans, GLL_to_GL_interp_trans, W_V_Dl_in, *\/ */
/* /\*                                      nodes_Lobatto, nodes_Gauss, nodes_Lobatto, nodes_Gauss, nodes_Lobatto, nodes_Gauss); *\/ */
/* /\*         } *\/ */
/* /\*         else if (dim == 2) { *\/ */
/* /\*           linalg_kron_vec_o_vec_dot_xy(Gauss_weights, Jdrdxdrdx[l][lp], V_Dl_in, nodes_Gauss, W_V_Dl_in); *\/ */
/* /\*           linalg_kron_A1A2x_NONSQR(VT_W_V_Dl_in, GLL_to_GL_interp_trans, GLL_to_GL_interp_trans, W_V_Dl_in, nodes_Lobatto, nodes_Gauss, *\/ */
/* /\*                                    nodes_Lobatto, nodes_Gauss); *\/ */
/* /\*         } *\/ */
/* /\*         else if (dim == 1){ *\/ */
/* /\*           linalg_kron_vec_dot_xy(Gauss_weights, Jdrdxdrdx[l][lp], V_Dl_in, nodes_Gauss, W_V_Dl_in); *\/ */
/* /\*           linalg_matvec_plus_vec(1.0, GLL_to_GL_interp_trans, W_V_Dl_in, 0., VT_W_V_Dl_in, nodes_Lobatto, nodes_Gauss); *\/ */
/* /\*         } *\/ */
          
/* /\*         else {     *\/ */
/* /\*           P4EST_FREE(DTlp_VT_W_V_Dl_in); *\/ */
/* /\*           P4EST_FREE(VT_W_V_Dl_in); *\/ */
/* /\*           P4EST_FREE(W_V_Dl_in); *\/ */
/* /\*           P4EST_FREE(V_Dl_in); *\/ */
/* /\*           P4EST_FREE(Dl_in); *\/ */
/* /\*           mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension."); *\/ */
/* /\*         } *\/ */

/* /\*         d4est_operators_apply_Dij_transpose(d4est_ops, VT_W_V_Dl_in, dim, deg_Lobatto, lp, DTlp_VT_W_V_Dl_in); *\/ */
/* /\*         linalg_vec_axpy(1., DTlp_VT_W_V_Dl_in, out, volume_nodes_Lobatto); *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */
/* /\*   } *\/ */
  
/* /\*   P4EST_FREE(DTlp_VT_W_V_Dl_in); *\/ */
/* /\*   P4EST_FREE(VT_W_V_Dl_in); *\/ */
/* /\*   P4EST_FREE(W_V_Dl_in); *\/ */
/* /\*   P4EST_FREE(V_Dl_in); *\/ */
/* /\*   P4EST_FREE(Dl_in); *\/ */
/* /\* } *\/ */


/* /\* integral(f(x,y,z) lj) where f is some function to integrate on the quadrature nodes and lj is the j-th lagrange function*\/ */
/* void d4est_operators_apply_curved_galerkin_quadral_old */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* in_quad, */
/*  int deg_lobatto, */
/*  double* jac_quad, */
/*  int deg_quad, */
/*  quadrature_type_t quad_type, */
/*  int dim, */
/*  double* out */
/* ){ */
/*   mpi_assert(dim == 1 || dim == 2 || dim == 3); */
/*   int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad); */

/*   /\* double* in_quad = P4EST_ALLOC(double, volume_nodes_quad); *\/ */
/*   double* w_j_in_quad = P4EST_ALLOC(double, volume_nodes_quad); */

/*   double* interp_trans = NULL; */
/*   double* quad_weights = NULL; */

/*   if (quad_type == QUAD_GAUSS){ */
/*     quad_weights = d4est_operators_fetch_GL_weights_1d(d4est_ops, deg_quad); */
/*     interp_trans = d4est_operators_fetch_ref_GLL_to_GL_interp_trans_1d */
/*                                    (d4est_ops, deg_lobatto, deg_quad);   */
/*   } */
/*   else if (quad_type == QUAD_LOBATTO){ */
/*     quad_weights = d4est_operators_fetch_GLL_weights_1d(d4est_ops, deg_quad); */
/*     interp_trans = d4est_operators_fetch_p_prolong_transpose_1d(d4est_ops, deg_lobatto, deg_quad); */
/*   } */
/*   else { */
/*     mpi_abort("[D4EST_ERROR]: Not a support quad type"); */
/*   } */
  
/*   int nodes_lobatto = deg_lobatto + 1; */
/*   int nodes_quad = deg_quad + 1; */
  
/*   if (dim == 1){ */
/*     linalg_kron_vec_dot_xy(quad_weights, jac_quad, in_quad, nodes_quad, w_j_in_quad); */
/*     linalg_matvec_plus_vec(1.0, interp_trans, w_j_in_quad, 0., out, nodes_lobatto, nodes_quad); */
/*   } */
/*   else if (dim == 2) { */
/*     linalg_kron_vec_o_vec_dot_xy(quad_weights, jac_quad, in_quad, nodes_quad, w_j_in_quad); */
/*     linalg_kron_A1A2x_NONSQR(out, interp_trans, interp_trans, w_j_in_quad, nodes_lobatto, nodes_quad, */
/*                              nodes_lobatto, nodes_quad); */
/*   } */
/*   else if (dim == 3) { */
/*     linalg_kron_vec_o_vec_o_vec_dot_xy(quad_weights, jac_quad, in_quad, nodes_quad, w_j_in_quad); */
/*     linalg_kron_A1A2A3x_NONSQR(out, interp_trans, interp_trans, interp_trans, w_j_in_quad, */
/*                                nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad); */
/*   } */
/*   else {     */
/*     P4EST_FREE(w_j_in_quad); */
/*     P4EST_FREE(in_quad); */
/*     mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension."); */
/*   }   */
/*   P4EST_FREE(w_j_in_quad); */
/* } */


/* void d4est_operators_compute_curved_inverse_mass_matrix */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  int deg_lobatto, */
/*  int deg_quad, */
/*  quadrature_type_t quad_type, */
/*  int dim, */
/*  double* jac_quad, */
/*  double* invM */
/* ) */
/* { */
/*   int volume_nodes_lobatto = d4est_operators_get_nodes(dim,deg_lobatto); */
  
/*   double* M = P4EST_ALLOC(double, volume_nodes_lobatto*volume_nodes_lobatto); */
/*   double* u = P4EST_ALLOC_ZERO(double, volume_nodes_lobatto); */
/*   double* Mu = P4EST_ALLOC_ZERO(double, volume_nodes_lobatto); */

/*   for (int i = 0; i < volume_nodes_lobatto; i++){ */
/*     u[i] = 1.; */
/*     d4est_operators_apply_curved_mass_matrix(d4est_ops, u, deg_lobatto, jac_quad, deg_quad, quad_type, dim, Mu); */
/*     linalg_set_column(M, Mu, i, volume_nodes_lobatto, volume_nodes_lobatto); */
/*     u[i] = 0.; */
/*   } */
/*   linalg_invert_and_copy(M, invM, volume_nodes_lobatto); */

/*   P4EST_FREE(Mu); */
/*   P4EST_FREE(u); */
/*   P4EST_FREE(M); */
/* } */



void d4est_operators_compute_prolong_matrix
(
 d4est_operators_t* d4est_ops,
 int degH,
 int dim,
 int* degh,
 int children,
 double* prolong_mat
)
{
  int volume_nodes_h = 0;
  for (int i = 0; i < children; i++){
    volume_nodes_h += d4est_operators_get_nodes(dim, degh[i]);
  }
  int volume_nodes_H = d4est_operators_get_nodes(dim, degH);

  double* u = P4EST_ALLOC_ZERO(double, volume_nodes_H);
  double* Mu = P4EST_ALLOC_ZERO(double, volume_nodes_h);

  for (int i = 0; i < volume_nodes_H; i++){
    u[i] = 1.;
    if (children == (P4EST_CHILDREN)){
      d4est_operators_apply_hp_prolong(d4est_ops, u, degH, dim, degh, Mu);
    }
    else {
      d4est_operators_apply_p_prolong(d4est_ops, u, degH, dim, degh[0], Mu);
    }
    linalg_set_column(prolong_mat, Mu, i, volume_nodes_h, volume_nodes_H);
    u[i] = 0.;
  }

  P4EST_FREE(Mu);
  P4EST_FREE(u);
}


void d4est_operators_compute_PT_mat_P
(
 d4est_operators_t* d4est_ops,
 double* mat, /* dimensions of the degh space \sumi(deghi+1)^DIM x \sumi(deghi+1)^DIM */
 int degH,
 int dim,
 int* degh,
 int children,
 double* PT_mat_P /* dimensions of the degH space: (degH+1)^DIM x (degH+1)^DIM */
)
{
  int volume_nodes_h [P4EST_CHILDREN];
  int total_volume_nodes_h = 0;
  int max_volume_nodes_h = -1;
  
  for (int i = 0; i < children; i++){
    volume_nodes_h[i] = d4est_operators_get_nodes(dim, degh[i]);
    total_volume_nodes_h += volume_nodes_h[i];
    max_volume_nodes_h = (volume_nodes_h[i] > max_volume_nodes_h) ? volume_nodes_h[i] : max_volume_nodes_h;
  }

  int volume_nodes_H = d4est_operators_get_nodes(dim, degH);

  linalg_fill_vec(PT_mat_P, 0., volume_nodes_H*volume_nodes_H);
  /* DEBUG_PRINT_ARR_DBL(PT_mat_P, volume_nodes_H*volume_nodes_H); */

  /* printf("total_volume_nodes_h*volume_nodes_H = %d\n",total_volume_nodes_h*volume_nodes_H); */
  /* printf("volume_nodes_H*volume_nodes_H = %d\n",volume_nodes_H*volume_nodes_H); */
  
  double* P = P4EST_ALLOC(double, total_volume_nodes_h*volume_nodes_H);
  double* PT = P4EST_ALLOC(double, total_volume_nodes_h*volume_nodes_H);


  
  double* mat_P_i = P4EST_ALLOC(double, max_volume_nodes_h*volume_nodes_H);
  double* PT_mat_P_i = P4EST_ALLOC(double, max_volume_nodes_h*volume_nodes_H);


  
  d4est_operators_compute_prolong_matrix(d4est_ops, degH, dim, degh, children, P);
  linalg_mat_transpose_nonsqr(P, PT, total_volume_nodes_h, volume_nodes_H);

  int stride_mat = 0;
  int stride_P = 0;
  for (int i = 0; i < children; i++){
    linalg_mat_multiply(&mat[stride_mat],
                        &P[stride_P],
                        mat_P_i,
                        volume_nodes_h[i],
                        volume_nodes_h[i],
                        volume_nodes_H);
    
    linalg_mat_multiply(&PT[stride_P], mat_P_i, PT_mat_P_i, volume_nodes_H, volume_nodes_h[i], volume_nodes_H);
    linalg_vec_axpy(1., PT_mat_P_i, PT_mat_P, volume_nodes_H*volume_nodes_H);
    
    stride_mat += volume_nodes_h[i]*volume_nodes_h[i];
    stride_P += volume_nodes_h[i]*volume_nodes_H;
  }

  P4EST_FREE(PT_mat_P_i);
  P4EST_FREE(mat_P_i);
  P4EST_FREE(PT);
  P4EST_FREE(P);
}

/* void d4est_operators_form_fofufofvlilj_matrix */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* u, */
/*  double* v, */
/*  int deg_lobatto, */
/*  double* xyz_quad [(P4EST_DIM)], */
/*  double* jac_quad, */
/*  int deg_quad, */
/*  quadrature_type_t quad_type, */
/*  int dim, */
/*  double* out, */
/*  grid_fcn_ext_t fofu_fcn, */
/*  void* fofu_ctx, */
/*  grid_fcn_ext_t fofv_fcn, */
/*  void* fofv_ctx */
/* ) */
/* { */
/*   if (fofu_fcn == NULL) */
/*     fofu_fcn = identity_fcn; */
/*   if (fofv_fcn == NULL) */
/*     fofv_fcn = identity_fcn; */
  
/*   double* u_quad = NULL; */
/*   double* v_quad = NULL; */

/*   int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad); */
  
/*   if (u != NULL){ */
/*     u_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*     d4est_operators_interp(d4est_ops, u, QUAD_LOBATTO, deg_lobatto, u_quad, quad_type, deg_quad, dim); */
/*   } */
/*   if (v != NULL){ */
/*     v_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*     d4est_operators_interp(d4est_ops, v, QUAD_LOBATTO, deg_lobatto, v_quad, quad_type, deg_quad, dim); */
/*   } */
  
/*   double* fofu_fofv_jac = P4EST_ALLOC(double, volume_nodes_quad); */
/*   for (int i = 0; i < volume_nodes_quad; i++){ */
/*     fofu_fofv_jac[i] = jac_quad[i]; */
/*     if (u != NULL || fofu_fcn != identity_fcn){ */
/*       fofu_fofv_jac[i] *= fofu_fcn(xyz_quad[0][i], */
/*                                    xyz_quad[1][i], */
/* #if (P4EST_DIM)==3 */
/*                                    xyz_quad[2][i], */
/* #endif */
/*                                    (u != NULL) ? u_quad[i] : 0, */
/*                                    fofu_ctx); */
/*     } */
/*     if (v != NULL || fofv_fcn != identity_fcn){ */
/*       fofu_fofv_jac[i] *= fofv_fcn(xyz_quad[0][i], */
/*                                    xyz_quad[1][i], */
/* #if (P4EST_DIM)==3 */
/*                                    xyz_quad[2][i], */
/* #endif */
/*                                    (v != NULL) ? v_quad[i] : 0, */
/*                                    fofv_ctx); */
/*     } */
/*   } */
  
/*   d4est_operators_compute_curved_mass_matrix */
/*     ( */
/*      d4est_ops, */
/*      deg_lobatto, */
/*      deg_quad, */
/*      quad_type, */
/*      dim, */
/*      fofu_fofv_jac, */
/*      out */
/*     ); */
  
/*   if (u != NULL){ */
/*     P4EST_FREE(u_quad); */
/*   } */
/*   if (v != NULL){ */
/*     P4EST_FREE(v_quad); */
/*   } */
  
/* } */


/* void d4est_operators_apply_fofufofvlilj */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* vec, */
/*  double* u, */
/*  double* v, */
/*  int deg_lobatto, */
/*  double* jac_quad, */
/*  double* xyz_quad [(P4EST_DIM)], */
/*  int deg_quad, */
/*  quadrature_type_t quad_type, */
/*  int dim, */
/*  double* Mvec, */
/*  grid_fcn_ext_t fofu_fcn, */
/*  void* fofu_ctx, */
/*  grid_fcn_ext_t fofv_fcn, */
/*  void* fofv_ctx */
/* ) */
/* { */
/*   if (fofu_fcn == NULL){ */
/*     fofu_fcn = identity_fcn; */
/*   } */
/*   if (fofv_fcn == NULL){ */
/*     fofv_fcn = identity_fcn; */
/*   } */
  
/*   double* u_quad = NULL; */
/*   double* v_quad = NULL; */

/*   int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad); */
  
/*   if (u != NULL){ */
/*     u_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*     d4est_operators_interp(d4est_ops, u, QUAD_LOBATTO, deg_lobatto, u_quad, quad_type, deg_quad, dim); */
/*   } */
/*   if (v != NULL){ */
/*     v_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*     d4est_operators_interp(d4est_ops, v, QUAD_LOBATTO, deg_lobatto, v_quad, quad_type, deg_quad, dim); */
/*   } */
  
/*   double* fofu_fofv_jac = P4EST_ALLOC(double, volume_nodes_quad); */
/*   for (int i = 0; i < volume_nodes_quad; i++){ */
/*     fofu_fofv_jac[i] = jac_quad[i]; */
/*     if (u != NULL || fofu_fcn != identity_fcn){ */
/*       fofu_fofv_jac[i] *= fofu_fcn(xyz_quad[0][i], */
/*                                    xyz_quad[1][i], */
/* #if (P4EST_DIM)==3 */
/*                                    xyz_quad[2][i], */
/* #endif */
/*                                    (u != NULL) ? u_quad[i] : 0, */
/*                                    fofu_ctx); */
/*     } */
/*     if (v != NULL || fofv_fcn != identity_fcn){ */
/*       fofu_fofv_jac[i] *= fofv_fcn(xyz_quad[0][i], */
/*                                    xyz_quad[1][i], */
/* #if (P4EST_DIM)==3 */
/*                                    xyz_quad[2][i], */
/* #endif */
/*                                    (v != NULL) ? v_quad[i] : 0, */
/*                                    fofv_ctx); */
/*     } */
/*   } */

/*   d4est_operators_apply_curved_mass_matrix */
/*     ( */
/*      d4est_ops, */
/*      vec, */
/*      deg_lobatto, */
/*      fofu_fofv_jac, */
/*      deg_quad, */
/*      quad_type, */
/*      dim, */
/*      Mvec */
/*     ); */
  
/*   P4EST_FREE(u_quad); */
/*   P4EST_FREE(v_quad); */
/*   P4EST_FREE(fofu_fofv_jac); */
  
/* } */

/* void d4est_operators_apply_fofufofvlj */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* u, */
/*  double* v, */
/*  int deg_lobatto, */
/*  double* jac_quad, */
/*  double* xyz_quad [(P4EST_DIM)], */
/*  int deg_quad, */
/*  quadrature_type_t quad_type, */
/*  int dim, */
/*  double* out, */
/*  grid_fcn_ext_t fofu_fcn, */
/*  void* fofu_ctx, */
/*  grid_fcn_ext_t fofv_fcn, */
/*  void* fofv_ctx */
/* ) */
/* { */
/*   if (fofu_fcn == NULL){ */
/*     fofu_fcn = identity_fcn; */
/*   } */
/*   if (fofv_fcn == NULL){ */
/*     fofv_fcn = identity_fcn; */
/*   } */
  
/*   double* u_quad = NULL; */
/*   double* v_quad = NULL; */

/*   int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad); */
  
/*   if (u != NULL){ */
/*     u_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*     d4est_operators_interp(d4est_ops, u, QUAD_LOBATTO, deg_lobatto, u_quad, quad_type, deg_quad, dim); */
/*   } */
/*   if (v != NULL){ */
/*     v_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*     d4est_operators_interp(d4est_ops, v, QUAD_LOBATTO, deg_lobatto, v_quad, quad_type, deg_quad, dim); */
/*   } */
  
/*   double* fofu_fofv = P4EST_ALLOC(double, volume_nodes_quad); */
/*   for (int i = 0; i < volume_nodes_quad; i++){ */
/*     fofu_fofv[i] = 1.0; */
/*     if (u != NULL || fofu_fcn != identity_fcn){ */
/*       fofu_fofv[i] *= fofu_fcn(xyz_quad[0][i], */
/*                                    xyz_quad[1][i], */
/* #if (P4EST_DIM)==3 */
/*                                    xyz_quad[2][i], */
/* #endif */
/*                                    (u != NULL) ? u_quad[i] : 0, */
/*                                    fofu_ctx); */
/*       /\* printf("xyz_quad[0][i], xyz_quad[1][i], xyz_quad[2][i], fofu_fofv[i] = %f,%f,%f,%f\n",xyz_quad[0][i], xyz_quad[1][i], xyz_quad[2][i], fofu_fofv[i]); *\/ */
/*     } */
/*     if (v != NULL || fofv_fcn != identity_fcn){ */
/*       fofu_fofv[i] *= fofv_fcn(xyz_quad[0][i], */
/*                                    xyz_quad[1][i], */
/* #if (P4EST_DIM)==3 */
/*                                    xyz_quad[2][i], */
/* #endif */
/*                                    (v != NULL) ? v_quad[i] : 0, */
/*                                    fofv_ctx); */
/*     } */

/*   } */
 
/*   d4est_operators_apply_curved_galerkin_quadral */
/*     ( */
/*      d4est_ops, */
/*      fofu_fofv, */
/*      deg_lobatto, */
/*      jac_quad, */
/*      deg_quad, */
/*      quad_type, */
/*      dim, */
/*      out */
/*     ); */

/*   /\* DEBUG_PRINT_2ARR_DBL(fofu_fofv,jac_quad, d4est_operators_get_nodes((P4EST_DIM), deg_quad)); *\/ */
/*   /\* DEBUG_PRINT_ARR_DBL(out, d4est_operators_get_nodes((P4EST_DIM), deg_lobatto)); *\/ */

  
/*   P4EST_FREE(u_quad); */
/*   P4EST_FREE(v_quad); */
/*   P4EST_FREE(fofu_fofv); */
  
/* } */


/* void d4est_operators_apply_curvedInverseGaussMass */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* in, */
/*  int deg_Lobatto, */
/*  double* jac_Gauss, */
/*  int deg_Gauss, */
/*  int dim, */
/*  double* out */
/* ){ */
/*   mpi_assert(dim == 1 || dim == 2 || dim == 3); */
/*   int volume_nodes_Gauss = d4est_operators_get_nodes(dim, deg_Gauss); */

/*   double* in_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss); */
/*   double* one_over_w_j_in_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss); */
  
/*   double* Gauss_weights = d4est_operators_fetch_GL_weights_1d(d4est_ops, deg_Gauss); */
/*   double* GLL_to_GL_interp_trans_inverse = d4est_operators_fetch_ref_GLL_to_GL_interp_trans_1d_inverse(d4est_ops, deg_Lobatto, deg_Gauss); */
/*   double* GLL_to_GL_interp_inverse = d4est_operators_fetch_ref_GLL_to_GL_interp_1d_inverse(d4est_ops, deg_Lobatto, deg_Gauss); */
  
/*   int nodes_Lobatto = deg_Lobatto + 1; */
/*   int nodes_Gauss = deg_Gauss + 1; */
  
/*   if (dim == 1){ */
/*     linalg_matvec_plus_vec(1.0, GLL_to_GL_interp_trans_inverse, in, 0., in_Gauss, nodes_Gauss, nodes_Lobatto); */
/*     linalg_kron_oneover_vec_dot_oneover_x_dot_y(Gauss_weights, jac_Gauss, in_Gauss, nodes_Gauss, one_over_w_j_in_Gauss); */
/*     linalg_matvec_plus_vec(1.0, GLL_to_GL_interp_inverse, one_over_w_j_in_Gauss, 0., out, nodes_Lobatto, nodes_Gauss); */
/*   } */
  
/*   else if (dim == 2){ */
/*     linalg_kron_A1A2x_NONSQR */
/*       ( */
/*        in_Gauss, */
/*        GLL_to_GL_interp_trans_inverse, */
/*        GLL_to_GL_interp_trans_inverse, */
/*        in, */
/*        nodes_Gauss, */
/*        nodes_Lobatto, */
/*        nodes_Gauss, */
/*        nodes_Lobatto */
/*       ); */
    
/*     linalg_kron_oneover_vec_o_vec_dot_oneover_x_dot_y */
/*       ( */
/*        Gauss_weights, */
/*        jac_Gauss, */
/*        in_Gauss, */
/*        nodes_Gauss, */
/*        one_over_w_j_in_Gauss */
/*       ); */
    
/*     linalg_kron_A1A2x_NONSQR */
/*       ( */
/*        out, */
/*        GLL_to_GL_interp_inverse, */
/*        GLL_to_GL_interp_inverse, */
/*        one_over_w_j_in_Gauss, */
/*        nodes_Lobatto, */
/*        nodes_Gauss, */
/*        nodes_Lobatto, */
/*        nodes_Gauss */
/*       ); */
/*   } */
  
/*   else if (dim == 3) { */
/*     linalg_kron_A1A2A3x_NONSQR( */
/*                                in_Gauss, */
/*                                GLL_to_GL_interp_trans_inverse, */
/*                                GLL_to_GL_interp_trans_inverse, */
/*                                GLL_to_GL_interp_trans_inverse, */
/*                                in, */
/*                                nodes_Gauss, */
/*                                nodes_Lobatto, */
/*                                nodes_Gauss, */
/*                                nodes_Lobatto, */
/*                                nodes_Gauss, */
/*                                nodes_Lobatto */
/*                               ); */
    
/*     linalg_kron_oneover_vec_o_vec_o_vec_dot_oneover_x_dot_y */
/*       ( */
/*        Gauss_weights, */
/*        jac_Gauss, */
/*        in_Gauss, */
/*        nodes_Gauss, */
/*        one_over_w_j_in_Gauss */
/*       ); */
    
/*     linalg_kron_A1A2A3x_NONSQR( */
/*                                out, */
/*                                GLL_to_GL_interp_inverse, */
/*                                GLL_to_GL_interp_inverse, */
/*                                GLL_to_GL_interp_inverse, */
/*                                one_over_w_j_in_Gauss, */
/*                                nodes_Lobatto, */
/*                                nodes_Gauss, */
/*                                nodes_Lobatto, */
/*                                nodes_Gauss, */
/*                                nodes_Lobatto, */
/*                                nodes_Gauss */
/*                               ); */
/*   } else { */
/*     P4EST_FREE(one_over_w_j_in_Gauss); */
/*     P4EST_FREE(in_Gauss); */
/*     mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension."); */
/*   }   */
/*   P4EST_FREE(one_over_w_j_in_Gauss); */
/*   P4EST_FREE(in_Gauss); */
/* } */



/* void d4est_operators_apply_curvedInverseLobattoMass */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* in, */
/*  int deg_Lobatto, */
/*  double* jac_Lobatto_quad, */
/*  int deg_Lobatto_quad, */
/*  int dim, */
/*  double* out */
/* ){ */
/*   mpi_assert(dim == 1 || dim == 2 || dim == 3); */
/*   int volume_nodes_Lobatto_quad = d4est_operators_get_nodes(dim, deg_Lobatto_quad); */

/*   double* in_Lobatto_quad = P4EST_ALLOC(double, volume_nodes_Lobatto_quad); */
/*   double* one_over_w_j_in_Lobatto_quad = P4EST_ALLOC(double, volume_nodes_Lobatto_quad); */
  
/*   double* Lobatto_quad_weights = d4est_operators_fetch_GLL_weights_1d(d4est_ops, deg_Lobatto_quad); */
/*   double* p_prolong_trans_inverse = d4est_operators_fetch_p_prolong_transpose_1d_inverse(d4est_ops, deg_Lobatto, deg_Lobatto_quad); */
/*   double* p_prolong_inverse = d4est_operators_fetch_p_prolong_1d_inverse(d4est_ops, deg_Lobatto, deg_Lobatto_quad); */
  
/*   int nodes_Lobatto = deg_Lobatto + 1; */
/*   int nodes_Lobatto_quad = deg_Lobatto_quad + 1; */
  
/*   if (dim == 1){ */
/*     linalg_matvec_plus_vec(1.0, p_prolong_trans_inverse, in, 0., in_Lobatto_quad, nodes_Lobatto_quad, nodes_Lobatto); */
/*     linalg_kron_oneover_vec_dot_oneover_x_dot_y(Lobatto_quad_weights, jac_Lobatto_quad, in_Lobatto_quad, nodes_Lobatto_quad, one_over_w_j_in_Lobatto_quad); */
/*     linalg_matvec_plus_vec(1.0, p_prolong_inverse, one_over_w_j_in_Lobatto_quad, 0., out, nodes_Lobatto, nodes_Lobatto_quad); */
/*   } */
  
/*   else if (dim == 2){ */
/*     linalg_kron_A1A2x_NONSQR */
/*       ( */
/*        in_Lobatto_quad, */
/*        p_prolong_trans_inverse, */
/*        p_prolong_trans_inverse, */
/*        in, */
/*        nodes_Lobatto_quad, */
/*        nodes_Lobatto, */
/*        nodes_Lobatto_quad, */
/*        nodes_Lobatto */
/*       ); */
    
/*     linalg_kron_oneover_vec_o_vec_dot_oneover_x_dot_y */
/*       ( */
/*        Lobatto_quad_weights, */
/*        jac_Lobatto_quad, */
/*        in_Lobatto_quad, */
/*        nodes_Lobatto_quad, */
/*        one_over_w_j_in_Lobatto_quad */
/*       ); */
    
/*     linalg_kron_A1A2x_NONSQR */
/*       ( */
/*        out, */
/*        p_prolong_inverse, */
/*        p_prolong_inverse, */
/*        one_over_w_j_in_Lobatto_quad, */
/*        nodes_Lobatto, */
/*        nodes_Lobatto_quad, */
/*        nodes_Lobatto, */
/*        nodes_Lobatto_quad */
/*       ); */
/*   } */
  
/*   else if (dim == 3) { */
/*     linalg_kron_A1A2A3x_NONSQR( */
/*                                in_Lobatto_quad, */
/*                                p_prolong_trans_inverse, */
/*                                p_prolong_trans_inverse, */
/*                                p_prolong_trans_inverse, */
/*                                in, */
/*                                nodes_Lobatto_quad, */
/*                                nodes_Lobatto, */
/*                                nodes_Lobatto_quad, */
/*                                nodes_Lobatto, */
/*                                nodes_Lobatto_quad, */
/*                                nodes_Lobatto */
/*                               ); */
    
/*     linalg_kron_oneover_vec_o_vec_o_vec_dot_oneover_x_dot_y */
/*       ( */
/*        Lobatto_quad_weights, */
/*        jac_Lobatto_quad, */
/*        in_Lobatto_quad, */
/*        nodes_Lobatto_quad, */
/*        one_over_w_j_in_Lobatto_quad */
/*       ); */
    
/*     linalg_kron_A1A2A3x_NONSQR( */
/*                                out, */
/*                                p_prolong_inverse, */
/*                                p_prolong_inverse, */
/*                                p_prolong_inverse, */
/*                                one_over_w_j_in_Lobatto_quad, */
/*                                nodes_Lobatto, */
/*                                nodes_Lobatto_quad, */
/*                                nodes_Lobatto, */
/*                                nodes_Lobatto_quad, */
/*                                nodes_Lobatto, */
/*                                nodes_Lobatto_quad */
/*                               ); */
/*   } else { */
/*     P4EST_FREE(one_over_w_j_in_Lobatto_quad); */
/*     P4EST_FREE(in_Lobatto_quad); */
/*     mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension."); */
/*   }   */
/*   P4EST_FREE(one_over_w_j_in_Lobatto_quad); */
/*   P4EST_FREE(in_Lobatto_quad); */
/* } */


/* void d4est_operators_interp_to_GL_nodes */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  double* u_Lobatto_in, */
/*  int deg_Lobatto, */
/*  int deg_Gauss, */
/*  double* u_Gauss_out, */
/*  int dim */
/* ) */
/* { */
/*   mpi_assert(dim == 1 || dim == 2 || dim == 3); */
/*   double* invVij = d4est_operators_fetch_invVij_1d(d4est_ops, deg_Lobatto); */
/*   double* GaussVij = d4est_operators_fetch_ref_GaussVij_1d(d4est_ops, deg_Lobatto, deg_Gauss); */
/*   double* temp = P4EST_ALLOC(double, d4est_operators_get_nodes((dim), deg_Lobatto)); */

/*   int nodes_Lobatto = deg_Lobatto + 1; */
/*   int nodes_Gauss = deg_Gauss + 1; */
  
/*   if (dim == 1){ */
/*     linalg_matvec_plus_vec(1.0, invVij, u_Lobatto_in, 0., temp, nodes_Lobatto, nodes_Lobatto); */
/*     linalg_matvec_plus_vec(1.0, GaussVij, temp, 0., u_Gauss_out, nodes_Gauss, nodes_Lobatto); */
/*   } */
/*   else if (dim == 2) { */
/*     linalg_kron_A1A2x_NONSQR(temp, invVij, invVij, u_Lobatto_in, nodes_Lobatto, nodes_Lobatto, */
/*                              nodes_Lobatto, nodes_Lobatto);     */
/*     linalg_kron_A1A2x_NONSQR(u_Gauss_out, GaussVij, GaussVij, temp, nodes_Gauss, nodes_Lobatto, */
/*                              nodes_Gauss, nodes_Lobatto); */
/*   } */
/*   else if (dim == 3){ */
/*     linalg_kron_A1A2A3x_NONSQR(temp, invVij, invVij, invVij, u_Lobatto_in, */
/*                                nodes_Lobatto, nodes_Lobatto, nodes_Lobatto, nodes_Lobatto, nodes_Lobatto, nodes_Lobatto); */
/*     linalg_kron_A1A2A3x_NONSQR(u_Gauss_out, GaussVij, GaussVij, GaussVij, temp, */
/*                                nodes_Gauss, nodes_Lobatto, nodes_Gauss, nodes_Lobatto, nodes_Gauss, nodes_Lobatto); */
/*   } */
/*   else { */
/*     P4EST_FREE(temp); */
/*     mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension."); */
/*   } */
/*   P4EST_FREE(temp); */
/* } */




void d4est_operators_hp_apply_nd_restrict_with_ptr(double* uH, int degH, double* uh,
                                          int degh, int dim, int c,
                                          double* hp_restrict_matrix_1d) {
  /* printf("c = %d, (1 << dim) = %d, degh = %d, degH = %d \n",c, (1 << dim),
   * degh, degH); */
  mpi_assert((degH <= degh) && (c < (1 << dim)));

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  if (dim == 1)
    linalg_matvec_plus_vec(1.0, &hp_restrict_matrix_1d[c * nodesh * nodesH], uh,
                           0., uH, nodesH, nodesh);
  else if (dim == 2) {
    int cx = d4est_operators_is_child_left_or_right(c, 0);
    int cy = d4est_operators_is_child_left_or_right(c, 1);

#if (ORDERING==1)
    linalg_kron_A1A2x_NONSQR(uH, &hp_restrict_matrix_1d[cx * nodesh * nodesH],
                             &hp_restrict_matrix_1d[cy * nodesh * nodesH], uh,
                             nodesH, nodesh, nodesH, nodesh);
    /* mpi_abort("Ordering = 1"); */
#endif
#if (ORDERING==3)
    linalg_kron_A1A2x_NONSQR(uH, &hp_restrict_matrix_1d[cy * nodesh * nodesH],
                             &hp_restrict_matrix_1d[cx * nodesh * nodesH], uh,
                             nodesH, nodesh, nodesH, nodesh);
#endif
  } else if (dim == 3) {
    int cx = d4est_operators_is_child_left_or_right(c, 0);
    int cy = d4est_operators_is_child_left_or_right(c, 1);
    int cz = d4est_operators_is_child_left_or_right(c, 2);
#if (ORDERING==1)
    linalg_kron_A1A2A3x_NONSQR(uH, &hp_restrict_matrix_1d[cx * nodesh * nodesH],
                               &hp_restrict_matrix_1d[cy * nodesh * nodesH],
                               &hp_restrict_matrix_1d[cz * nodesh * nodesH], uh,
                               nodesH, nodesh, nodesH, nodesh, nodesH, nodesh);
    /* mpi_abort("Ordering = 1"); */
#endif
#if (ORDERING==3)
    linalg_kron_A1A2A3x_NONSQR(uH, &hp_restrict_matrix_1d[cz * nodesh * nodesH],
                               &hp_restrict_matrix_1d[cy * nodesh * nodesH],
                               &hp_restrict_matrix_1d[cx * nodesh * nodesH], uh,
                               nodesH, nodesh, nodesH, nodesh, nodesH, nodesh);
#endif
    
  } else {
    mpi_abort("ERROR: d4est_operators_restrict_nd_U\n");
  }
}

/* double* d4est_operators_fetch_GLL_1d(d4est_operators_t* d4est_ops, int deg) { */
/*   /\* printf("deg,d4est_ops->degmath_max_storage = %d,%d" ,deg,d4est_ops->d4est_operators_max_storage); *\/ */
/*   mpi_assert(deg < d4est_ops->d4est_operators_max_storage); */
/*   if (deg > d4est_ops->d4est_operators_max_degree_used) */
/*     d4est_ops->d4est_operators_max_degree_used = deg; */
/*   return &(d4est_operators_GLLnodes[deg - 1][0]); */
/* } */

/* TODO: put this into the database */


/* TODO: put this into the database */
static void d4est_operators_build_dVij_1d(d4est_operators_t* d4est_ops, double* dVij_1d,
                                 int deg) {
  int i, j, rows, cols;
  double* lobatto_nodes = d4est_operators_fetch_GLL_nodes_1d(d4est_ops, deg);
  rows = cols = deg + 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      dVij_1d[i * cols + j] = d4est_operators_gradjacobi(lobatto_nodes[i], 0., 0., j);
}

/* TODO: put this into the database */
static void d4est_operators_build_dGaussVij_1d(d4est_operators_t* d4est_ops, double* dVij_1d,
                                      int Lobatto_degree, int Gauss_degree) {

  double* Gauss_nodes = d4est_operators_fetch_GL_nodes_1d(d4est_ops, Gauss_degree);
  
  int rows = Gauss_degree + 1;
  int cols = Lobatto_degree + 1;
  
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      dVij_1d[i * cols + j] = d4est_operators_gradjacobi(Gauss_nodes[i], 0., 0., j);
}

void d4est_operators_build_Mij_1d(d4est_operators_t* d4est_ops, double* Mij_1d,
                                int deg) {
  int edge_nodes = deg + 1;
  double* v1d = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);
  double* v1d_trans = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);
  d4est_operators_build_Vij_1d(d4est_ops, v1d, deg);
  linalg_mat_transpose(v1d, v1d_trans, edge_nodes);
  linalg_mat_multiply(v1d, v1d_trans, Mij_1d, edge_nodes, edge_nodes,
                      edge_nodes);
  linalg_invert(Mij_1d, edge_nodes);
  P4EST_FREE(v1d_trans);
  P4EST_FREE(v1d);
}

void
d4est_operators_build_GLL_nodes_and_weights_1d(d4est_operators_t* d4est_ops, double* GLL_nodes, double* GLL_weights, int deg)
{
  d4est_operators_GLL_nodes_and_weights(deg+1, GLL_nodes, GLL_weights);
}


double* d4est_operators_fetch_GLL_nodes_1d(d4est_operators_t* d4est_ops, int deg) {
  mpi_assert(deg < d4est_ops->d4est_operators_max_storage);
  if (d4est_ops->d4est_operators_GLL_nodes_1d_table[deg] != NULL) {
    return d4est_ops->d4est_operators_GLL_nodes_1d_table[deg];
  }
  else {
    int size = (deg + 1);
    d4est_ops->d4est_operators_GLL_nodes_1d_table[deg] = P4EST_ALLOC(double, size);
    d4est_ops->d4est_operators_GLL_weights_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_GLL_nodes_1d_table[deg];
    double* op1 = d4est_ops->d4est_operators_GLL_weights_1d_table[deg];
    d4est_operators_build_GLL_nodes_and_weights_1d(d4est_ops, op, op1, deg);
    return op;
  }
}

double* d4est_operators_fetch_GLL_weights_1d(d4est_operators_t* d4est_ops, int deg) {
  mpi_assert(deg < d4est_ops->d4est_operators_max_storage);
  if (d4est_ops->d4est_operators_GLL_weights_1d_table[deg] != NULL) {
    return d4est_ops->d4est_operators_GLL_weights_1d_table[deg];
  }
  else {
    int size = (deg + 1);
    d4est_ops->d4est_operators_GLL_nodes_1d_table[deg] = P4EST_ALLOC(double, size);
    d4est_ops->d4est_operators_GLL_weights_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_GLL_nodes_1d_table[deg];
    double* op1 = d4est_ops->d4est_operators_GLL_weights_1d_table[deg];
    d4est_operators_build_GLL_nodes_and_weights_1d(d4est_ops, op, op1, deg);
    return op1;
  }
}

void
d4est_operators_build_GL_nodes_and_weights_1d(d4est_operators_t* d4est_ops, double* GL_nodes, double* GL_weights, int deg)
{
  d4est_operators_GL_nodes_and_weights(deg+1, GL_nodes, GL_weights);
}

double* d4est_operators_fetch_GL_nodes_1d(d4est_operators_t* d4est_ops, int deg) {
  mpi_assert(deg < d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_GL_nodes_1d_table[deg] != NULL) {
    return d4est_ops->d4est_operators_GL_nodes_1d_table[deg];
  }
  else {
    int size = (deg + 1);
    d4est_ops->d4est_operators_GL_nodes_1d_table[deg] = P4EST_ALLOC(double, size);
    d4est_ops->d4est_operators_GL_weights_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_GL_nodes_1d_table[deg];
    double* op1 = d4est_ops->d4est_operators_GL_weights_1d_table[deg];
    d4est_operators_build_GL_nodes_and_weights_1d(d4est_ops, op, op1, deg);
    return op;
  }
}

double* d4est_operators_fetch_GL_weights_1d(d4est_operators_t* d4est_ops, int deg) {
  mpi_assert(deg < d4est_ops->d4est_operators_max_storage);
  if (d4est_ops->d4est_operators_GL_weights_1d_table[deg] != NULL) {
    return d4est_ops->d4est_operators_GL_weights_1d_table[deg];
  }
  else {
    int size = (deg + 1);
    d4est_ops->d4est_operators_GL_nodes_1d_table[deg] = P4EST_ALLOC(double, size);
    d4est_ops->d4est_operators_GL_weights_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_GL_nodes_1d_table[deg];
    double* op1 = d4est_ops->d4est_operators_GL_weights_1d_table[deg];
    d4est_operators_build_GL_nodes_and_weights_1d(d4est_ops, op, op1, deg);
    return op1;
  }
}



/* double* d4est_operators_fetch_Mij_1d(d4est_operators_t* d4est_ops, int deg) { */
/*   mpi_assert(deg < d4est_ops->d4est_operators_max_storage); */
/*   int* last_stride = &(d4est_ops->d4est_operators_ref_Mij_1d_table_last_stride); */
/*   if (d4est_ops->d4est_operators_ref_Mij_1d_table[deg] != -1) { */
/*     return &(d4est_ops->d4est_operators_ref_Mij_1d[d4est_ops->d4est_operators_ref_Mij_1d_table[deg]]); */
/*   } */

/*   else { */
/*     int size = (deg + 1) * (deg + 1); */
/*     d4est_ops->d4est_operators_ref_Mij_1d = */
/*         P4EST_REALLOC(d4est_ops->d4est_operators_ref_Mij_1d, double, *last_stride + size); */
/*     d4est_operators_build_Mij_1d(d4est_ops, &(d4est_ops->d4est_operators_ref_Mij_1d)[*last_stride], */
/*                         deg); */
/*     d4est_ops->d4est_operators_ref_Mij_1d_table[deg] = *last_stride; */
/*     *last_stride += size; */
/*     return &(d4est_ops->d4est_operators_ref_Mij_1d)[*last_stride - size]; */
/*   } */
/* } */

double* d4est_operators_fetch_Mij_1d(d4est_operators_t* d4est_ops, int deg) {
  mpi_assert(deg < d4est_ops->d4est_operators_max_storage);
  if (d4est_ops->d4est_operators_ref_Mij_1d_table[deg] != NULL) {
    return d4est_ops->d4est_operators_ref_Mij_1d_table[deg];
  }
  else {
    int size = (deg + 1) * (deg + 1);
    d4est_ops->d4est_operators_ref_Mij_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_ref_Mij_1d_table[deg];
    d4est_operators_build_Mij_1d(d4est_ops, op, deg);
    return op;
  }
}



/* double* d4est_operators_fetch_Mij_1d(d4est_operators_t* d4est_ops, int deg) { */
/*   mpi_assert(deg < d4est_ops->d4est_operators_max_storage); */
/*   if (d4est_ops->d4est_operators_ref_Mij_1d_table[deg] != NULL) { */
/*     return d4est_ops->d4est_operators_ref_Mij_1d_table[deg]; */
/*   } */
/*   else { */
/*     int size = (deg + 1) * (deg + 1); */
/*     d4est_ops->d4est_operators_ref_Mij_1d_table[deg] = P4EST_ALLOC(double, size); */
/*     double* op = d4est_ops->d4est_operators_ref_Mij_1d_table[deg]; */
/*     d4est_operators_build_Mij_1d(d4est_ops, op, deg); */
/*     return op; */
/*   } */
/* } */




static void d4est_operators_build_invMij_1d(d4est_operators_t* d4est_ops,
                                   double* invMij_1d, int deg) {
  double* Mij_1d = d4est_operators_fetch_Mij_1d(d4est_ops, deg);
  linalg_invert_and_copy(Mij_1d, invMij_1d, deg + 1);
}

static void d4est_operators_build_Dij_1d(d4est_operators_t* d4est_ops, double* Dij_1d,
                                int deg) {
  int edge_nodes = deg + 1;
  double* v1d = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);
  double* inv_v1d = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);
  double* v1dr = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);

  d4est_operators_build_Vij_1d(d4est_ops, v1d, deg);
  d4est_operators_build_dVij_1d(d4est_ops, v1dr, deg);
  linalg_invert_and_copy(v1d, inv_v1d, edge_nodes);

  linalg_mat_multiply(v1dr, inv_v1d, Dij_1d, edge_nodes, edge_nodes,
                      edge_nodes);

  P4EST_FREE(v1dr);
  P4EST_FREE(inv_v1d);
  P4EST_FREE(v1d);
}


static void d4est_operators_build_GaussDij_1d(d4est_operators_t* d4est_ops, double* Dij_1d,
                                     int deg_Lobatto, int deg_Gauss) {
  int nodes_Lobatto = deg_Lobatto + 1;
  int nodes_Gauss = deg_Gauss + 1;
  
  double* inv_v1d = (double*)P4EST_ALLOC(double, nodes_Lobatto*nodes_Lobatto);
  double* v1dr = (double*)P4EST_ALLOC(double, nodes_Gauss*nodes_Lobatto);

  d4est_operators_build_Vij_1d(d4est_ops, inv_v1d, deg_Lobatto);
  linalg_invert(inv_v1d, nodes_Lobatto);
  
  d4est_operators_build_dGaussVij_1d(d4est_ops, v1dr, deg_Lobatto, deg_Gauss);

  linalg_mat_multiply(v1dr, inv_v1d, Dij_1d, nodes_Gauss, nodes_Lobatto,
                      nodes_Lobatto);

  P4EST_FREE(v1dr);
  P4EST_FREE(inv_v1d);
}

static double* d4est_operators_fetch_invMij_1d(d4est_operators_t* d4est_ops, int deg){
  mpi_assert(deg < d4est_ops->d4est_operators_max_storage);
  if (d4est_ops->d4est_operators_ref_invMij_1d_table[deg] != NULL) {
    return d4est_ops->d4est_operators_ref_invMij_1d_table[deg];
  }
  else {
    int size = (deg + 1) * (deg + 1);
    d4est_ops->d4est_operators_ref_invMij_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_ref_invMij_1d_table[deg];
    d4est_operators_build_invMij_1d(d4est_ops, op, deg);
    return op;
  }
}

void d4est_operators_apply_Mij(d4est_operators_t* d4est_ops, double* in, int dim, int deg,
                      double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  double* mass_ref_1d = d4est_operators_fetch_Mij_1d(d4est_ops, deg);

  int nodes = deg + 1;
  if (dim == 1)
    linalg_matvec_plus_vec(1.0, mass_ref_1d, in, 0., out, nodes, nodes);
  else if (dim == 2) {
    /* linalg_kron_A1A2x_SQR(out, mass_ref_1d, mass_ref_1d, in, nodes, nodes);
     */
    linalg_kron_A1A2x_NONSQR(out, mass_ref_1d, mass_ref_1d, in, nodes, nodes,
                             nodes, nodes);
  } else if (dim == 3) {
    /* linalg_kron_A1A2A3x_SQR(out, mass_ref_1d, mass_ref_1d, mass_ref_1d, in,
     * nodes, nodes, nodes); */
    linalg_kron_A1A2A3x_NONSQR(out, mass_ref_1d, mass_ref_1d, mass_ref_1d, in,
                               nodes, nodes, nodes, nodes, nodes, nodes);
  } else {
    mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension.");
  }
}

void d4est_operators_apply_invMij(d4est_operators_t* d4est_ops, double* in, int dim,
                         int deg, double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  double* inv_mass_ref_1d = d4est_operators_fetch_invMij_1d(d4est_ops, deg);

  int nodes = deg + 1;
  if (dim == 1)
    linalg_matvec_plus_vec(1.0, inv_mass_ref_1d, in, 0., out, nodes, nodes);
  else if (dim == 2) {
    linalg_kron_A1A2x_NONSQR(out, inv_mass_ref_1d, inv_mass_ref_1d, in, nodes,
                             nodes, nodes, nodes);
  } else if (dim == 3) {
    linalg_kron_A1A2A3x_NONSQR(out, inv_mass_ref_1d, inv_mass_ref_1d,
                               inv_mass_ref_1d, in, nodes, nodes, nodes, nodes,
                               nodes, nodes);
  } else {
    mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension.");
  }
}

/* static double* d4est_operators_fetch_Dij_1d(d4est_operators_t* d4est_ops, int deg) { */
/*   mpi_assert(deg < d4est_ops->d4est_operators_max_storage); */
/*   int* last_stride = &(d4est_ops->d4est_operators_ref_Dij_1d_table_last_stride); */
/*   if (d4est_ops->d4est_operators_ref_Dij_1d_table[deg] != -1) { */
/*     return &(d4est_ops->d4est_operators_ref_Dij_1d[d4est_ops->d4est_operators_ref_Dij_1d_table[deg]]); */
/*   } else { */
/*     int size = (deg + 1) * (deg + 1); */
/*     d4est_ops->d4est_operators_ref_Dij_1d = */
/*         P4EST_REALLOC(d4est_ops->d4est_operators_ref_Dij_1d, double, *last_stride + size); */

/*     d4est_operators_build_Dij_1d(d4est_ops, &(d4est_ops->d4est_operators_ref_Dij_1d)[*last_stride], */
/*                         deg); */
/*     d4est_ops->d4est_operators_ref_Dij_1d_table[deg] = *last_stride; */
/*     *last_stride += size; */
/*     return &(d4est_ops->d4est_operators_ref_Dij_1d)[*last_stride - size]; */
/*   } */
/* } */


static double* d4est_operators_fetch_Dij_1d(d4est_operators_t* d4est_ops, int deg) {
  mpi_assert(deg < d4est_ops->d4est_operators_max_storage);
  if (d4est_ops->d4est_operators_ref_Dij_1d_table[deg] != NULL) {
    return d4est_ops->d4est_operators_ref_Dij_1d_table[deg];
  }
  else {
    int size = (deg + 1) * (deg + 1);
    d4est_ops->d4est_operators_ref_Dij_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_ref_Dij_1d_table[deg];
    d4est_operators_build_Dij_1d(d4est_ops, op, deg);
    return op;
  }
}

double* d4est_operators_fetch_GaussDij_1d(d4est_operators_t* d4est_ops, int deg_Lobatto,
                                         int deg_Gauss) {
  mpi_assert(deg_Gauss <= d4est_ops->d4est_operators_max_storage);
  mpi_assert(deg_Lobatto <= d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_ref_GaussDij_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] != NULL) {
    return d4est_ops->d4est_operators_ref_GaussDij_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss];
  }
  else {
    int size = (deg_Gauss + 1) * (deg_Lobatto + 1);
    d4est_ops->d4est_operators_ref_GaussDij_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_ref_GaussDij_1d_table[d4est_ops->d4est_operators_max_storage * deg_Lobatto + deg_Gauss];
    d4est_operators_build_GaussDij_1d(d4est_ops, op, deg_Lobatto, deg_Gauss);
    return op;
  }
}

  void
    d4est_operators_apply_GaussDij(d4est_operators_t* d4est_ops, double* in, int dim, int deg_Lobatto, int deg_Gauss, int dir, double* out)
{
  /* printf("dir = %d, dim = %d\n", dir, dim); */

    mpi_assert(dir < dim && dir >= 0 && (dim == 1 || dim == 3 || dim == 2));

    /* TODO get rid of this */
  /* mpi_assert(deg_Lobatto == deg_Gauss); */
  
  double* Dr_1d = d4est_operators_fetch_GaussDij_1d(d4est_ops, deg_Lobatto, deg_Gauss);
  int nodes = deg_Lobatto + 1;
  
  if (dim == 1)
    linalg_matvec_plus_vec(1.0, Dr_1d, in, 0., out, nodes, nodes);
  else if (dim == 2) {
#if ORDERING == 1
    if (dir == 0) linalg_kron_MAToIx_SQR(out, Dr_1d, in, nodes);
    if (dir == 1) linalg_kron_IoMATx_SQR(out, Dr_1d, in, nodes);
    /* mpi_abort("ORDERING == 1"); */
#endif
#if ORDERING == 3
    if (dir == 0) linalg_kron_IoMATx_SQR(out, Dr_1d, in, nodes);
    if (dir == 1) linalg_kron_MAToIx_SQR(out, Dr_1d, in, nodes);
#endif
    
  } else if (dim == 3) {
#if ORDERING == 1
    if (dir == 0) linalg_kron_MAToIoIx_SQR(out, Dr_1d, in, nodes);
    if (dir == 1) linalg_kron_IoMAToIx_SQR(out, Dr_1d, in, nodes);
    if (dir == 2) linalg_kron_IoIoMATx_SQR(out, Dr_1d, in, nodes);
    /* mpi_abort("Ordering == 1"); */
#endif
#if ORDERING == 2
    /* if (dir == 0) linalg_kron_IoMAToIx_SQR(out, Dr_1d, in, nodes); */
    /* if (dir == 1) linalg_kron_IoIoMATx_SQR(out, Dr_1d, in, nodes); */
    /* if (dir == 2) linalg_kron_MAToIoIx_SQR(out, Dr_1d, in, nodes); */
    mpi_abort("Ordering == 2");
#endif
#if ORDERING == 3
    if (dir == 0) linalg_kron_IoIoMATx_SQR(out, Dr_1d, in, nodes);
    if (dir == 1) linalg_kron_IoMAToIx_SQR(out, Dr_1d, in, nodes);
    if (dir == 2) linalg_kron_MAToIoIx_SQR(out, Dr_1d, in, nodes);
#endif
  } else {
    mpi_abort("ERROR: d4est_operators_apply_GaussDij");
  }
}

void d4est_operators_apply_Dij(d4est_operators_t* d4est_ops, double* in, int dim, int deg,
                      int dir, double* out) {
  mpi_assert(dir < dim && dir >= 0 && (dim == 1 || dim == 3 || dim == 2));
  double* Dr_1d = d4est_operators_fetch_Dij_1d(d4est_ops, deg);
  int nodes = deg + 1;
  if (dim == 1)
    linalg_matvec_plus_vec(1.0, Dr_1d, in, 0., out, nodes, nodes);
  else if (dim == 2) {
#if ORDERING == 1
    if (dir == 0) linalg_kron_MAToIx_SQR(out, Dr_1d, in, nodes);
    if (dir == 1) linalg_kron_IoMATx_SQR(out, Dr_1d, in, nodes);
    /* mpi_abort("ORDERING == 1"); */
#endif
#if ORDERING == 3
    if (dir == 0) linalg_kron_IoMATx_SQR(out, Dr_1d, in, nodes);
    if (dir == 1) linalg_kron_MAToIx_SQR(out, Dr_1d, in, nodes);
#endif
    
  } else if (dim == 3) {
#if ORDERING == 1
    if (dir == 0) linalg_kron_MAToIoIx_SQR(out, Dr_1d, in, nodes);
    if (dir == 1) linalg_kron_IoMAToIx_SQR(out, Dr_1d, in, nodes);
    if (dir == 2) linalg_kron_IoIoMATx_SQR(out, Dr_1d, in, nodes);
    /* mpi_abort("Ordering == 1"); */
#endif
#if ORDERING == 2
    /* if (dir == 0) linalg_kron_IoMAToIx_SQR(out, Dr_1d, in, nodes); */
    /* if (dir == 1) linalg_kron_IoIoMATx_SQR(out, Dr_1d, in, nodes); */
    /* if (dir == 2) linalg_kron_MAToIoIx_SQR(out, Dr_1d, in, nodes); */
    mpi_abort("Ordering == 2");
#endif
#if ORDERING == 3
    if (dir == 0) linalg_kron_IoIoMATx_SQR(out, Dr_1d, in, nodes);
    if (dir == 1) linalg_kron_IoMAToIx_SQR(out, Dr_1d, in, nodes);
    if (dir == 2) linalg_kron_MAToIoIx_SQR(out, Dr_1d, in, nodes);
#endif
  } else {
    mpi_abort("ERROR: d4est_operators_apply_Dij");
  }
}


void d4est_operators_apply_Dij_transpose(d4est_operators_t* d4est_ops, double* in, int dim, int deg,
                      int dir, double* out) {
  mpi_assert(dir < dim && dir >= 0 && (dim == 1 || dim == 3 || dim == 2));
  double* Dr_1d = d4est_operators_fetch_Dij_1d(d4est_ops, deg);

  /* bad change this later */
  double* Dr_1d_transpose = P4EST_ALLOC(double, (deg+1)*(deg+1));

  linalg_mat_transpose(Dr_1d, Dr_1d_transpose, (deg+1));
  int nodes = deg + 1;
  if (dim == 1)
    linalg_matvec_plus_vec(1.0, Dr_1d_transpose, in, 0., out, nodes, nodes);
  else if (dim == 2) {
#if ORDERING == 1
    if (dir == 0) linalg_kron_MAToIx_SQR(out, Dr_1d_transpose, in, nodes);
    if (dir == 1) linalg_kron_IoMATx_SQR(out, Dr_1d_transpose, in, nodes);
    /* mpi_abort("ORDERING == 1"); */
#endif
#if ORDERING == 3
    if (dir == 0) linalg_kron_IoMATx_SQR(out, Dr_1d_transpose, in, nodes);
    if (dir == 1) linalg_kron_MAToIx_SQR(out, Dr_1d_transpose, in, nodes);
#endif
    
  } else if (dim == 3) {
#if ORDERING == 1
    if (dir == 0) linalg_kron_MAToIoIx_SQR(out, Dr_1d_transpose, in, nodes);
    if (dir == 1) linalg_kron_IoMAToIx_SQR(out, Dr_1d_transpose, in, nodes);
    if (dir == 2) linalg_kron_IoIoMATx_SQR(out, Dr_1d_transpose, in, nodes);
    /* mpi_abort("Ordering == 1"); */
#endif
#if ORDERING == 2
    /* if (dir == 0) linalg_kron_IoMAToIx_SQR(out, Dr_1d_transpose, in, nodes); */
    /* if (dir == 1) linalg_kron_IoIoMATx_SQR(out, Dr_1d_transpose, in, nodes); */
    /* if (dir == 2) linalg_kron_MAToIoIx_SQR(out, Dr_1d_transpose, in, nodes); */
    mpi_abort("Ordering == 2");
#endif
#if ORDERING == 3
    if (dir == 0) linalg_kron_IoIoMATx_SQR(out, Dr_1d_transpose, in, nodes);
    if (dir == 1) linalg_kron_IoMAToIx_SQR(out, Dr_1d_transpose, in, nodes);
    if (dir == 2) linalg_kron_MAToIoIx_SQR(out, Dr_1d_transpose, in, nodes);
#endif
  } else {
    P4EST_FREE(Dr_1d_transpose);
    mpi_abort("ERROR: d4est_operators_apply_Dij");
  }

  P4EST_FREE(Dr_1d_transpose);
}

static void d4est_operators_build_hp_prolong_1d_aux(int degH, int degh, int c,
                                           double* inv_v1d_trans_degH,
                                           double* lobatto_nodes_degh,
                                           double* hp_prolong_matrix_1d) {
  mpi_assert(degH <= degh);
  int n, i;
  int nodesH = degH + 1;
  int nodesh = degh + 1;
  double* phi_modal = (double*)P4EST_ALLOC(double, nodesH);
  double r;
  for (n = 0; n < nodesh; n++) {
    r = lobatto_nodes_degh[n];
    if (c == 0 || c == 1) d4est_operators_child_to_parent_ref_coords(c, &r);
    for (i = 0; i < nodesH; i++) {
      phi_modal[i] = d4est_operators_jacobi(r, 0, 0, i);
    }

    linalg_matvec_plus_vec(1.0, inv_v1d_trans_degH, phi_modal, 0.0,
                           &hp_prolong_matrix_1d[n * nodesH], nodesH, nodesH);
  }
  P4EST_FREE(phi_modal);
}

static void d4est_operators_build_hp_prolong_1d(d4est_operators_t* d4est_ops,
                                       double* hp_prolong_1d, int degH,
                                       int degh) {
  int edge_nodes_degH = degH + 1;
  int edge_nodes_degh = degh + 1;

  double* v1d = (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH);
  double* inv_v1d =
      (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH);
  double* inv_v1d_trans_degH =
      (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH);

  d4est_operators_build_Vij_1d(d4est_ops, v1d, degH);
  linalg_invert_and_copy(v1d, inv_v1d, edge_nodes_degH);
  linalg_mat_transpose(inv_v1d, inv_v1d_trans_degH, edge_nodes_degH);
  double* lobatto_nodes_degh = d4est_operators_fetch_GLL_nodes_1d(d4est_ops, degh);

  d4est_operators_build_hp_prolong_1d_aux(degH, degh, 0, inv_v1d_trans_degH,
                                 lobatto_nodes_degh, &hp_prolong_1d[0]);
  d4est_operators_build_hp_prolong_1d_aux(
      degH, degh, 1, inv_v1d_trans_degH, lobatto_nodes_degh,
      &hp_prolong_1d[edge_nodes_degH * edge_nodes_degh]);

  P4EST_FREE(inv_v1d);
  P4EST_FREE(inv_v1d_trans_degH);
  P4EST_FREE(v1d);
}

void d4est_operators_build_p_prolong_1d(d4est_operators_t* d4est_ops,
                                      double* p_prolong_1d, int degH,
                                      int degh) {
  /* int edge_nodes_degH = degH + 1; */
  /* /\* int edge_nodes_degh = degh + 1; *\/ */

  /* double* v1d = (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH); */
  /* double* inv_v1d = */
  /*     (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH); */
  /* double* inv_v1d_trans_degH = */
  /*     (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH); */

  /* d4est_operators_build_Vij_1d(d4est_ops, v1d, degH); */
  /* linalg_invert_and_copy(v1d, inv_v1d, edge_nodes_degH); */
  /* linalg_mat_transpose(inv_v1d, inv_v1d_trans_degH, edge_nodes_degH); */
  /* double* lobatto_nodes_degh = d4est_operators_fetch_GLL_nodes_1d(d4est_ops, degh); */

  /* d4est_operators_build_hp_prolong_1d_aux(degH, degh, -1, inv_v1d_trans_degH, */
  /*                                lobatto_nodes_degh, &p_prolong_1d[0]); */

  /* P4EST_FREE(inv_v1d); */
  /* P4EST_FREE(inv_v1d_trans_degH); */
  /* P4EST_FREE(v1d); */
  double* degh_nodes = d4est_operators_fetch_GLL_nodes_1d(d4est_ops, degh);
  double* invVij = d4est_operators_fetch_invVij_1d(d4est_ops, degH);
  double* Vij_degh = P4EST_ALLOC(double, (degH+1)*(degh+1));
  
  int rows = degh + 1;
  int cols = degH + 1;
  
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      Vij_degh[i * cols + j] = d4est_operators_jacobi(degh_nodes[i], 0., 0., j);

  linalg_mat_multiply(Vij_degh, invVij, p_prolong_1d, degh+1, degH+1, degH+1);
  P4EST_FREE(Vij_degh);
}




void d4est_operators_hp_apply_nd_prolong_transpose_with_ptr(
    double* uH, int degH, double* uh, int degh, int dim, int c,
    double* hp_prolong_transpose_matrix_1d) {
  /* sanity check */
  mpi_assert((degH <= degh) && (c < (1 << dim)));

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  if (dim == 1)
    linalg_matvec_plus_vec(1.0,
                           &hp_prolong_transpose_matrix_1d[c * nodesh * nodesH],
                           uh, 0., uH, nodesH, nodesh);
  else if (dim == 2) {
    int cx = d4est_operators_is_child_left_or_right(c, 0);
    int cy = d4est_operators_is_child_left_or_right(c, 1);
#if (ORDERING) == 1
    linalg_kron_A1A2x_NONSQR(
        uH, &hp_prolong_transpose_matrix_1d[cx * nodesh * nodesH],
        &hp_prolong_transpose_matrix_1d[cy * nodesh * nodesH], uh, nodesH,
        nodesh, nodesH, nodesh);
    /* mpi_abort("ORDERING == 1"); */
#endif
#if (ORDERING) == 3
    linalg_kron_A1A2x_NONSQR(
        uH, &hp_prolong_transpose_matrix_1d[cy * nodesh * nodesH],
        &hp_prolong_transpose_matrix_1d[cx * nodesh * nodesH], uh, nodesH,
        nodesh, nodesH, nodesh);
#endif    
  } else if (dim == 3) {
    int cx = d4est_operators_is_child_left_or_right(c, 0);
    int cy = d4est_operators_is_child_left_or_right(c, 1);
    int cz = d4est_operators_is_child_left_or_right(c, 2);
#if (ORDERING) == 1
    linalg_kron_A1A2A3x_NONSQR(
        uH, &hp_prolong_transpose_matrix_1d[cx * nodesh * nodesH],
        &hp_prolong_transpose_matrix_1d[cy * nodesh * nodesH],
        &hp_prolong_transpose_matrix_1d[cz * nodesh * nodesH], uh, nodesH,
        nodesh, nodesH, nodesh, nodesH, nodesh);
    /* mpi_abort("ORDERING==1"); */
#endif
#if (ORDERING) == 3
    linalg_kron_A1A2A3x_NONSQR(
        uH, &hp_prolong_transpose_matrix_1d[cz * nodesh * nodesH],
        &hp_prolong_transpose_matrix_1d[cy * nodesh * nodesH],
        &hp_prolong_transpose_matrix_1d[cx * nodesh * nodesH], uh, nodesH,
        nodesh, nodesH, nodesh, nodesH, nodesh);
#endif    
  } else {
    mpi_abort("ERROR: d4est_operators_prolong_transpose_nd_U\n");
  }
}

static double* d4est_operators_fetch_hp_prolong_1d(d4est_operators_t* d4est_ops, int degH,
                                          int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < d4est_ops->d4est_operators_max_storage);
  
  if (d4est_ops->d4est_operators_hp_prolong_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh] != NULL) {
    return d4est_ops->d4est_operators_hp_prolong_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh];
  }
  else {
    int size = 2 * (degh + 1) * (degH + 1);
    d4est_ops->d4est_operators_hp_prolong_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_hp_prolong_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh];
    d4est_operators_build_hp_prolong_1d(d4est_ops, op, degH, degh);
    return op;
  }
}

double* d4est_operators_fetch_p_prolong_1d(d4est_operators_t* d4est_ops, int degH,
                                         int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh <= d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_p_prolong_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh] != NULL) {
    return d4est_ops->d4est_operators_p_prolong_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh];
  }
  else {
    int size = (degh + 1) * (degH + 1);
    d4est_ops->d4est_operators_p_prolong_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_p_prolong_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh];
    d4est_operators_build_p_prolong_1d(d4est_ops, op, degH, degh);
    return op;
  }
}

double* d4est_operators_fetch_p_prolong_1d_inverse(d4est_operators_t* d4est_ops, int degH,
                                         int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh <= d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_p_prolong_1d_inverse_table[d4est_ops->d4est_operators_max_storage * degH + degh] != NULL) {
    return d4est_ops->d4est_operators_p_prolong_1d_inverse_table[d4est_ops->d4est_operators_max_storage * degH + degh];
  }
  else {
    int size = (degh + 1) * (degH + 1);
    d4est_ops->d4est_operators_p_prolong_1d_inverse_table[d4est_ops->d4est_operators_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_p_prolong_1d_inverse_table[d4est_ops->d4est_operators_max_storage * degH + degh];
    d4est_operators_build_p_prolong_1d_inverse(d4est_ops, op, degH, degh);
    return op;
  }
}


void d4est_operators_apply_hp_prolong(d4est_operators_t* d4est_ops, double* in, int degH,
                             int dim, int* degh, double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int children = (1 << dim);
  int c;
  double* hp_prolong_1d;

  int stride = 0;
  for (c = 0; c < children; c++) {
    hp_prolong_1d = d4est_operators_fetch_hp_prolong_1d(d4est_ops, degH, degh[c]);
    d4est_operators_hp_apply_nd_prolong_with_ptr(&out[stride], degh[c], in, degH, dim, c,
                                        hp_prolong_1d);
    stride += d4est_operators_get_nodes(dim, degh[c]);
  }
}

void d4est_operators_apply_p_prolong(d4est_operators_t* d4est_ops, double* in, int degH,
                            int dim, int degh, double* out) {
  double* p_prolong_1d = d4est_operators_fetch_p_prolong_1d(d4est_ops, degH, degh);

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  /* if degrees are the same, p_prolong matrix should be the identity */
  if (degh == degH) {
    linalg_copy_1st_to_2nd(in, out, d4est_operators_get_nodes(dim, degh));
    return;
  }

  if (dim == 1) {
    linalg_matvec_plus_vec(1.0, p_prolong_1d, in, 0., out, nodesh, nodesH);
  } else if (dim == 2) {
    linalg_kron_A1A2x_NONSQR(out, p_prolong_1d, p_prolong_1d, in, nodesh,
                             nodesH, nodesh, nodesH);
  } else if (dim == 3) {
    linalg_kron_A1A2A3x_NONSQR(out, p_prolong_1d, p_prolong_1d, p_prolong_1d,
                               in, nodesh, nodesH, nodesh, nodesH, nodesh,
                               nodesH);
  } else {
    mpi_abort("ERROR: d4est_operators_prolong_nd_U\n");
  }
}

void d4est_operators_build_hp_restrict_1d_aux(int degh, int degH,
                                            double* hp_prolong_matrix_1d,
                                            double* mass_matrix_rs_degh,
                                            double* inv_mass_matrix_rs_degH,
                                            double* hp_restrict_matrix_1d) {
  int nodesH = degH + 1;
  int nodesh = degh + 1;
  int s;

  double* tmp = P4EST_ALLOC_ZERO(double, nodesH);
  double* c_s_x_Mh = P4EST_ALLOC_ZERO(double, nodesh* nodesh);
  double* c_s = P4EST_ALLOC_ZERO(double, nodesh* nodesH);

  for (s = 0; s < nodesH; s++) {
    tmp[s] = 1.0;
    linalg_matvec_plus_vec(.5, hp_prolong_matrix_1d, tmp, 0.0, &c_s[s * nodesh],
                           nodesh, nodesH);
    tmp[s] = 0.0;
  }

  /* printf("CALLING DG_INTERP_HP_RESTRICT_1d\n"); */
  /* util_print_matrix(hp_prolong_matrix_1d, nodesh, nodesH,
   * "hp_prolong_matrix_1d = ", 0); */

  linalg_mat_multiply(c_s, mass_matrix_rs_degh, c_s_x_Mh, nodesH, nodesh,
                      nodesh);

  linalg_mat_multiply(inv_mass_matrix_rs_degH, c_s_x_Mh,
                      &hp_restrict_matrix_1d[0], nodesH, nodesH, nodesh);

  /* util_print_matrix(mass_matrix_rs_degh, nodesh, nodesh, "mass_matrix_rs_degh
   * = ", 0); */
  /* util_print_matrix(inv_mass_matrix_rs_degH, nodesH, nodesH,
   * "inv_mass_matrix_rs_degH = ", 0); */
  /* util_print_matrix(c_s, nodesh, nodesH, "c_s = ", 0); */
  /* util_print_matrix(c_s_x_Mh, nodesh, nodesH, "c_s_x_Mh = ", 0); */
  /* util_print_matrix(&hp_restrict_matrix_1d[0], nodesH, nodesh,
   * "&hp_restrict_matrix_1d[0] = ", 0); */

  P4EST_FREE(tmp);
  P4EST_FREE(c_s);
  P4EST_FREE(c_s_x_Mh);
}

static void d4est_operators_build_p_restrict_1d(d4est_operators_t* d4est_ops,
                                       double* p_restrict_1d, int degH,
                                       int degh) {
  int edge_nodes_degH = degH + 1;
  int edge_nodes_degh = degh + 1;

  /* printf("degH = %d\n",degH); */
  /* printf("degh = %d\n",degh); */

  double* p_prolong_1d = d4est_operators_fetch_p_prolong_1d(d4est_ops, degH, degh);
  double* inv_Mij_degH =
      (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH);

  /* util_print_matrix(Mij_degH, edge_nodes_degH, edge_nodes_degH, "Mij_degH =
   * ", 0); */
  /* util_print_matrix(inv_Mij_degH, edge_nodes_degH, edge_nodes_degH,
   * "Mij_degH_inv = ", 0); */
  double* Mij_degH = d4est_operators_fetch_Mij_1d(d4est_ops, degH);
  linalg_invert_and_copy(Mij_degH, inv_Mij_degH, edge_nodes_degH);

  double* Mij_degh = d4est_operators_fetch_Mij_1d(d4est_ops, degh);
  d4est_operators_build_hp_restrict_1d_aux(degh, degH, &p_prolong_1d[0], Mij_degh,
                                  inv_Mij_degH, &p_restrict_1d[0]);

  linalg_vec_scale(2., p_restrict_1d, edge_nodes_degH * edge_nodes_degh);

  /* util_print_matrix(p_restrict_1d, edge_nodes_degH, edge_nodes_degh,
   * "p_restrict_1d[0] = ", 0); */
  /* util_print_matrix(&p_restrict_1d[edge_nodes_degH*edge_nodes_degh],
   * edge_nodes_degH, edge_nodes_degh, "p_restrict_1d[1] = ", 0); */

  P4EST_FREE(inv_Mij_degH);
}

 double* d4est_operators_fetch_p_restrict_1d(d4est_operators_t* d4est_ops, int degH,
                                          int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_p_restrict_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh] != NULL) {
    return d4est_ops->d4est_operators_p_restrict_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh];
  }
  else {
    int size = (degh + 1) * (degH + 1);
    d4est_ops->d4est_operators_p_restrict_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_p_restrict_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh];
    d4est_operators_build_p_restrict_1d(d4est_ops, op, degH, degh);
    return op;
  }
}

void d4est_operators_apply_p_restrict(d4est_operators_t* d4est_ops, double* in, int degh,
                             int dim, int degH, double* out) {
  double* p_restrict_1d = d4est_operators_fetch_p_restrict_1d(d4est_ops, degH, degh);

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  /* restrict matrix should be the identity if degrees are the same */
  if (degh == degH) {
    linalg_copy_1st_to_2nd(in, out, d4est_operators_get_nodes(dim, degh));
    return;
  }

  if (dim == 1)
    linalg_matvec_plus_vec(1.0, p_restrict_1d, in, 0., out, nodesH, nodesh);
  else if (dim == 2) {
    linalg_kron_A1A2x_NONSQR(out, p_restrict_1d, p_restrict_1d, in, nodesH,
                             nodesh, nodesH, nodesh);
  } else if (dim == 3) {
    linalg_kron_A1A2A3x_NONSQR(out, p_restrict_1d, p_restrict_1d, p_restrict_1d,
                               in, nodesH, nodesh, nodesH, nodesh, nodesH,
                               nodesh);
  } else {
    mpi_abort("ERROR: d4est_operators_restrict_nd_U\n");
  }
}

static void d4est_operators_build_hp_restrict_1d(d4est_operators_t* d4est_ops,
                                        double* hp_restrict_1d, int degH,
                                        int degh) {
  int edge_nodes_degH = degH + 1;
  int edge_nodes_degh = degh + 1;

  /* printf("degH = %d\n",degH); */
  /* printf("degh = %d\n",degh); */

  double* hp_prolong_1d = d4est_operators_fetch_hp_prolong_1d(d4est_ops, degH, degh);

  double* inv_Mij_degH =
      (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH);

  /* TODO: replace this with a call to the database for inverse Mij */
  double* Mij_degH = d4est_operators_fetch_Mij_1d(d4est_ops, degH);
  linalg_invert_and_copy(Mij_degH, inv_Mij_degH, edge_nodes_degH);

  double* Mij_degh = d4est_operators_fetch_Mij_1d(d4est_ops, degh);
  d4est_operators_build_hp_restrict_1d_aux(degh, degH, &hp_prolong_1d[0], Mij_degh,
                                  inv_Mij_degH, &hp_restrict_1d[0]);

  d4est_operators_build_hp_restrict_1d_aux(
      degh, degH, &hp_prolong_1d[edge_nodes_degH * edge_nodes_degh], Mij_degh,
      inv_Mij_degH, &hp_restrict_1d[edge_nodes_degH * edge_nodes_degh]);

  P4EST_FREE(inv_Mij_degH);
}

static double* d4est_operators_fetch_hp_restrict_1d(d4est_operators_t* d4est_ops, int degH,
                                           int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_hp_restrict_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh] != NULL) {
    return d4est_ops->d4est_operators_hp_restrict_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh];
  }
  else {
    int size = 2 * (degh + 1) * (degH + 1);
    d4est_ops->d4est_operators_hp_restrict_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_hp_restrict_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh];
    d4est_operators_build_hp_restrict_1d(d4est_ops, op, degH, degh);
    return op;
  }
}

void d4est_operators_apply_hp_restrict(d4est_operators_t* d4est_ops, double* in, int* degh,
                              int dim, int degH, double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int children = (1 << dim);
  int c;
  int nodesH = d4est_operators_get_nodes(dim, degH);
  double* hp_restrict_1d;

  int stride = 0;
  double* tmp = P4EST_ALLOC(double, nodesH);
  linalg_fill_vec(out, 0., nodesH);

  for (c = 0; c < children; c++) {
    hp_restrict_1d = d4est_operators_fetch_hp_restrict_1d(d4est_ops, degH, degh[c]);
    /* int nodesh = d4est_operators_get_nodes(dim,degh[c]); */

    d4est_operators_hp_apply_nd_restrict_with_ptr(tmp, degH, &in[stride], degh[c], dim,
                                         c, hp_restrict_1d);
    /* printf("child = %d\n",c); */
    /* util_print_matrix(tmp, nodesH, 1, "tmp = ", 0); */
    /* util_print_matrix(&uh[stride], nodesh, 1, "uh[stride] = ", 0); */
    linalg_vec_axpy(1.0, tmp, out, nodesH);
    /* util_print_matrix(uH, nodesH, 1, "uH = ", 0); */
    stride += d4est_operators_get_nodes(dim, degh[c]);
  }

  P4EST_FREE(tmp);
}

int d4est_operators_fetch_max_degree_used(d4est_operators_t* d4est_ops) {
  return d4est_ops->d4est_operators_max_degree_used;
}

static void d4est_operators_build_xyz_nd(d4est_operators_t* d4est_ops, double* ref_xyz_nd,
                                int dim, int deg) {
  int nodes = deg + 1;
  int vol_nodes = d4est_operators_get_nodes(dim, deg);
  double* lgl = d4est_operators_fetch_GLL_nodes_1d(d4est_ops, deg);
  int i;

  double* eye = P4EST_ALLOC_ZERO(double, nodes);
  for (i = 0; i < nodes; i++) eye[i] = 1.;

  if (dim == 2) {
#if ORDERING == 1
    /* mpi_abort("ORDERING == 1"); */
    linalg_kron_AoB(lgl, eye, &ref_xyz_nd[0], (nodes), 1, (nodes), 1);
    linalg_kron_AoB(eye, lgl, &ref_xyz_nd[nodes * nodes], (nodes), 1, (nodes), 1);
#endif
#if ORDERING == 3
    linalg_kron_AoB(eye, lgl, &ref_xyz_nd[0], (nodes), 1, (nodes), 1);
    linalg_kron_AoB(lgl, eye, &ref_xyz_nd[nodes * nodes], (nodes), 1, (nodes), 1);
#endif
  } else if (dim == 3) {
#if ORDERING == 1
    linalg_kron_AoBoC(lgl, eye, eye, &ref_xyz_nd[0], nodes, 1, nodes, 1, nodes,
                      1);
    linalg_kron_AoBoC(eye, lgl, eye, &ref_xyz_nd[vol_nodes], nodes, 1, nodes, 1,
                      nodes, 1);
    linalg_kron_AoBoC(eye, eye, lgl, &ref_xyz_nd[2 * vol_nodes], nodes, 1,
                      nodes, 1, nodes, 1);
    /* mpi_abort("Ordering == 1"); */
#endif
#if ORDERING == 2
    /* linalg_kron_AoBoC(eye, lgl, eye, &ref_xyz_nd[0], nodes, 1, nodes, 1, nodes, */
    /*                   1); */
    /* linalg_kron_AoBoC(eye, eye, lgl, &ref_xyz_nd[vol_nodes], nodes, 1, nodes, 1, */
    /*                   nodes, 1); */
    /* linalg_kron_AoBoC(lgl, eye, eye, &ref_xyz_nd[2 * vol_nodes], nodes, 1, */
    /*                   nodes, 1, nodes, 1); */
    mpi_abort("Ordering == 2");
#endif
#if ORDERING == 3
    linalg_kron_AoBoC(eye, eye, lgl, &ref_xyz_nd[0], nodes, 1, nodes, 1, nodes,
                      1);
    linalg_kron_AoBoC(eye, lgl, eye, &ref_xyz_nd[vol_nodes], nodes, 1, nodes, 1,
                      nodes, 1);
    linalg_kron_AoBoC(lgl, eye, eye, &ref_xyz_nd[2 * vol_nodes], nodes, 1,
                      nodes, 1, nodes, 1);
#endif
  } else
    mpi_abort("ERROR 1: d4est_operators_build_xyz_nd");

  P4EST_FREE(eye);
}


static void d4est_operators_build_Gauss_xyz_nd(d4est_operators_t* d4est_ops, double* ref_xyz_nd,
                                int dim, int deg) {
  int nodes = deg + 1;
  int vol_nodes = d4est_operators_get_nodes(dim, deg);
  double* gl = d4est_operators_fetch_GL_nodes_1d(d4est_ops, deg);
  int i;

  double* eye = P4EST_ALLOC_ZERO(double, nodes);
  for (i = 0; i < nodes; i++) eye[i] = 1.;

  if (dim == 2) {
#if ORDERING == 1
    /* mpi_abort("ORDERING == 1"); */
    linalg_kron_AoB(gl, eye, &ref_xyz_nd[0], (nodes), 1, (nodes), 1);
    linalg_kron_AoB(eye, gl, &ref_xyz_nd[nodes * nodes], (nodes), 1, (nodes), 1);
#endif
#if ORDERING == 3
    linalg_kron_AoB(eye, gl, &ref_xyz_nd[0], (nodes), 1, (nodes), 1);
    linalg_kron_AoB(gl, eye, &ref_xyz_nd[nodes * nodes], (nodes), 1, (nodes), 1);
#endif
  } else if (dim == 3) {
#if ORDERING == 1
    linalg_kron_AoBoC(gl, eye, eye, &ref_xyz_nd[0], nodes, 1, nodes, 1, nodes,
                      1);
    linalg_kron_AoBoC(eye, gl, eye, &ref_xyz_nd[vol_nodes], nodes, 1, nodes, 1,
                      nodes, 1);
    linalg_kron_AoBoC(eye, eye, gl, &ref_xyz_nd[2 * vol_nodes], nodes, 1,
                      nodes, 1, nodes, 1);
    /* mpi_abort("Ordering == 1"); */
#endif
#if ORDERING == 2
    /* linalg_kron_AoBoC(eye, gl, eye, &ref_xyz_nd[0], nodes, 1, nodes, 1, nodes, */
    /*                   1); */
    /* linalg_kron_AoBoC(eye, eye, gl, &ref_xyz_nd[vol_nodes], nodes, 1, nodes, 1, */
    /*                   nodes, 1); */
    /* linalg_kron_AoBoC(gl, eye, eye, &ref_xyz_nd[2 * vol_nodes], nodes, 1, */
    /*                   nodes, 1, nodes, 1); */
    mpi_abort("Ordering == 2");
#endif
#if ORDERING == 3
    linalg_kron_AoBoC(eye, eye, gl, &ref_xyz_nd[0], nodes, 1, nodes, 1, nodes,
                      1);
    linalg_kron_AoBoC(eye, gl, eye, &ref_xyz_nd[vol_nodes], nodes, 1, nodes, 1,
                      nodes, 1);
    linalg_kron_AoBoC(gl, eye, eye, &ref_xyz_nd[2 * vol_nodes], nodes, 1,
                      nodes, 1, nodes, 1);
#endif
  } else
    mpi_abort("ERROR 1: d4est_operators_build_xyz_nd");

  P4EST_FREE(eye);
}

/* double* d4est_operators_fetch_xyz_nd(d4est_operators_t* d4est_ops, int dim, int deg, */
/*                             int dir) { */
/*   /\* printf("dim,deg,dir=%d,%d,%d\n",dim,deg,dir); *\/ */
/*   mpi_assert(deg < d4est_ops->d4est_operators_max_storage && (dim == 3 || dim == 2) && */
/*              (dir >= 0 && dir < dim)); */
/*   int* last_stride = &(d4est_ops->d4est_operators_ref_xyz_nd_table_last_stride); */
/*   if (d4est_ops->d4est_operators_ref_xyz_nd_table[dim * d4est_ops->d4est_operators_max_storage + deg] != */
/*       -1) { */
/*     return &( */
/*         d4est_ops->d4est_operators_ref_xyz_nd[d4est_ops->d4est_operators_ref_xyz_nd_table */
/*                                       [dim * d4est_ops->d4est_operators_max_storage + deg] + */
/*                                   dir * d4est_operators_get_nodes(dim, deg)]); */
/*   } else { */
/*     int size = (P4EST_DIM)*d4est_operators_get_nodes(dim, deg); */
/*     d4est_ops->d4est_operators_ref_xyz_nd = */
/*         P4EST_REALLOC(d4est_ops->d4est_operators_ref_xyz_nd, double, *last_stride + size); */

/*     d4est_operators_build_xyz_nd(d4est_ops, &(d4est_ops->d4est_operators_ref_xyz_nd)[*last_stride], dim, */
/*                         deg); */
/*     d4est_ops->d4est_operators_ref_xyz_nd_table[dim * d4est_ops->d4est_operators_max_storage + deg] = */
/*         *last_stride; */
/*     *last_stride += size; */
/*     return &(d4est_ops->d4est_operators_ref_xyz_nd)[*last_stride - size + */
/*                                         dir * d4est_operators_get_nodes(dim, deg)]; */
/*   } */
/*   /\* return NULL; *\/ */
/* } */


double* d4est_operators_fetch_xyz_nd(d4est_operators_t* d4est_ops, int dim, int deg,
                            int dir) {

  mpi_assert(deg < d4est_ops->d4est_operators_max_storage && (dim == 3 || dim == 2 || dim == 1) &&
             (dir >= 0 && dir < dim));

  
  if (dim == 1)
    return d4est_operators_fetch_GLL_nodes_1d(d4est_ops, deg);
  
  if (d4est_ops->d4est_operators_ref_xyz_nd_table[(dim-2) * d4est_ops->d4est_operators_max_storage + deg] != NULL){
    double* table = d4est_ops->d4est_operators_ref_xyz_nd_table[(dim-2) * d4est_ops->d4est_operators_max_storage + deg];
    return &table[dir * d4est_operators_get_nodes(dim, deg)];
  }
  
  else {
    int size = dim*d4est_operators_get_nodes(dim, deg);
    d4est_ops->d4est_operators_ref_xyz_nd_table[(dim-2)* d4est_ops->d4est_operators_max_storage + deg] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_ref_xyz_nd_table[(dim-2) * d4est_ops->d4est_operators_max_storage + deg];
    d4est_operators_build_xyz_nd(d4est_ops, op, dim, deg);
    return &op[dir * d4est_operators_get_nodes(dim, deg)];
  }
}



double* d4est_operators_fetch_Gauss_xyz_nd(d4est_operators_t* d4est_ops, int dim, int deg,
                            int dir) {

  mpi_assert(deg < d4est_ops->d4est_operators_max_storage && (dim == 3 || dim == 2 || dim == 1) &&
             (dir >= 0 && dir < dim));

  
  if (dim == 1)
    return d4est_operators_fetch_GL_nodes_1d(d4est_ops, deg);
  
  if (d4est_ops->d4est_operators_ref_Gauss_xyz_nd_table[(dim-2) * d4est_ops->d4est_operators_max_storage + deg] != NULL){
    double* table = d4est_ops->d4est_operators_ref_Gauss_xyz_nd_table[(dim-2) * d4est_ops->d4est_operators_max_storage + deg];
    return &table[dir * d4est_operators_get_nodes(dim, deg)];
  }
  
  else {
    int size = dim*d4est_operators_get_nodes(dim, deg);
    d4est_ops->d4est_operators_ref_Gauss_xyz_nd_table[(dim-2)* d4est_ops->d4est_operators_max_storage + deg] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_ref_Gauss_xyz_nd_table[(dim-2) * d4est_ops->d4est_operators_max_storage + deg];
    d4est_operators_build_Gauss_xyz_nd(d4est_ops, op, dim, deg);
    return &op[dir * d4est_operators_get_nodes(dim, deg)];
  }
}


static void d4est_operators_build_LIFT_1d(double* LIFT_1d, int deg) {
  memset(LIFT_1d, 0., sizeof(double) * 2 * (deg + 1));
  LIFT_1d[0] = 1.;
  LIFT_1d[2 * deg + 1] = 1.;
}

/* static double* d4est_operators_fetch_LIFT_1d(d4est_operators_t* d4est_ops, int deg) { */
  
/*   mpi_assert(deg < d4est_ops->d4est_operators_max_storage); */
/*   int* last_stride = &(d4est_ops->d4est_operators_LIFT_1d_table_last_stride); */

/*   if (d4est_ops->d4est_operators_LIFT_1d_table[deg] != -1) { */
/*     return &(d4est_ops->d4est_operators_LIFT_1d[d4est_ops->d4est_operators_LIFT_1d_table[deg]]); */
/*   } else { */
/*     int size = 2 * (deg + 1); */
/*     d4est_ops->d4est_operators_LIFT_1d = */
/*         P4EST_REALLOC(d4est_ops->d4est_operators_LIFT_1d, double, *last_stride + size); */
/*     d4est_operators_build_LIFT_1d(&(d4est_ops->d4est_operators_LIFT_1d)[*last_stride], deg); */
/*     d4est_ops->d4est_operators_LIFT_1d_table[deg] = *last_stride; */
/*     *last_stride += size; */
/*     return &(d4est_ops->d4est_operators_LIFT_1d)[*last_stride - size]; */
/*   } */
/* } */


static double* d4est_operators_fetch_LIFT_1d(d4est_operators_t* d4est_ops, int deg) {

  mpi_assert(deg < d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_LIFT_1d_table[deg] != NULL) {
    return d4est_ops->d4est_operators_LIFT_1d_table[deg];
  }
  else {
    int size = 2 * (deg + 1);
    d4est_ops->d4est_operators_LIFT_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_LIFT_1d_table[deg];
    d4est_operators_build_LIFT_1d(op, deg);
    return op;
  }
}

void d4est_operators_apply_LIFT(d4est_operators_t* d4est_ops, double* in, int dim, int deg,
                       int face, double* out) {
  mpi_assert(dim == 2 || dim == 3);
  mpi_assert(face < 2 * (dim));

  double* lift_1d = d4est_operators_fetch_LIFT_1d(d4est_ops, deg);

  /* int nodes = util_pow_int(deg+1, DIM); */

  int dir;  /* x <-> 0; y <-> 1; z <-> 2  */
  int side; /* left <-> 0; right <-> 1 */
  int nodes = deg + 1;

  if (face == 0) {
    dir = 0;
    side = 0;
  } else if (face == 1) {
    dir = 0;
    side = 1;
  } else if (face == 2) {
    dir = 1;
    side = 0;
  } else if (face == 3) {
    dir = 1;
    side = 1;
  } else if (face == 4) {
    dir = 2;
    side = 0;
  } else if (face == 5) {
    dir = 2;
    side = 1;
  } else {
    dir = -1;
    side = -1;
    mpi_abort("ERROR 0: d4est_operators_lift_boundary_vec");
  }

  if (dim == 2){
#if (ORDERING) == 1
    if (dir == 0)
      linalg_kron_VECoIx_SQR(out, &lift_1d[side * nodes], in, nodes);
    else if (dir == 1)
      linalg_kron_IoVECx_SQR(out, &lift_1d[side * nodes], in, nodes);
    else
      mpi_abort("dim == 2 so dir == 0 OR 1");
    /* mpi_abort("ORDERING == 1"); */
#endif
    
#if (ORDERING) == 3
    if (dir == 0)
      linalg_kron_IoVECx_SQR(out, &lift_1d[side * nodes], in, nodes);
    else if (dir == 1)
      linalg_kron_VECoIx_SQR(out, &lift_1d[side * nodes], in, nodes);
    else 
      mpi_abort("dim == 2 so dir == 0 OR 1");
#endif
    
  }


  else if (dim == 3){

#if ORDERING == 1
  /* mpi_abort("ORDERING == 1"); */
  if (dir == 0)
    linalg_kron_VECoIoIx_SQR(out, &lift_1d[side * nodes], in, nodes);
  else if (dir == 1)
    linalg_kron_IoVECoIx_SQR(out, &lift_1d[side * nodes], in, nodes);
  else if (dir == 2)
    linalg_kron_IoIoVECx_SQR(out, &lift_1d[side * nodes], in, nodes);
  else
    mpi_abort("DIM = 3 so DIR = 0,1,2");
#endif

#if ORDERING == 2
  mpi_abort("ORDERING == 2");
  /* else if (dim == 3 && dir == 0) */
  /*   linalg_kron_IoVECoIx_SQR(out, &lift_1d[side * nodes], in, nodes); */
  /* else if (dim == 3 && dir == 1) */
  /*   linalg_kron_IoIoVECx_SQR(out, &lift_1d[side * nodes], in, nodes); */
  /* else if (dim == 3 && dir == 2) */
  /*   linalg_kron_VECoIoIx_SQR(out, &lift_1d[side * nodes], in, nodes); */
#endif

  
#if ORDERING == 3
  if (dir == 0)
    linalg_kron_IoIoVECx_SQR(out, &lift_1d[side * nodes], in, nodes);
  else if (dir == 1)
    linalg_kron_IoVECoIx_SQR(out, &lift_1d[side * nodes], in, nodes);
  else if (dir == 2)
    linalg_kron_VECoIoIx_SQR(out, &lift_1d[side * nodes], in, nodes);
  else
    mpi_abort("DIM = 3 so DIR = 0,1,2");
#endif
  
  }

  else {
    mpi_abort("DIM not supported");
  }
}

void d4est_operators_apply_slicer(d4est_operators_t* d4est_ops, double* in, int dim,
                         int face, int deg, double* out) {
  mpi_assert(face < 2 * (dim));
  /* int nodes = util_pow_int(deg+1, DIM); */

  double* slicer_1d = d4est_operators_fetch_LIFT_1d(d4est_ops, deg);

  int dir;  /* x <-> 0; y <-> 1; z <-> 2  */
  int side; /* left <-> 0; right <-> 1 */

  if (face == 0) {
    dir = 0;
    side = 0;
  } else if (face == 1) {
    dir = 0;
    side = 1;
  } else if (face == 2) {
    dir = 1;
    side = 0;
  } else if (face == 3) {
    dir = 1;
    side = 1;
  } else if (face == 4) {
    dir = 2;
    side = 0;
  } else if (face == 5) {
    dir = 2;
    side = 1;
  } else {
    dir = -1;
    side = -1;
    mpi_abort("ERROR 0: d4est_operators_lift_boundary_vec");
  }

  /* util_print_matrix(slicer_1d, (deg+1), 1, "slicer_1d = ", 0); */

  if (dim == 2){
#if (ORDERING) == 1
    /* mpi_abort("ORDERING == 1"); */
    if (dir == 0)
      linalg_kron_VEC_TRANSoIx_SQR(out, &slicer_1d[side * (deg + 1)], in,
                                 deg + 1);
    else if (dir == 1)
      linalg_kron_IoVEC_TRANSx_SQR(out, &slicer_1d[side * (deg + 1)], in,
                                 deg + 1);
    else
      mpi_abort("DIM == 2, so dir == 0,1");
#endif
#if (ORDERING) == 3
    if (dir == 0)
      linalg_kron_IoVEC_TRANSx_SQR(out, &slicer_1d[side * (deg + 1)], in,
                                 deg + 1);
    else if (dir == 1)
      linalg_kron_VEC_TRANSoIx_SQR(out, &slicer_1d[side * (deg + 1)], in,
                                 deg + 1);
    else
      mpi_abort("DIM == 2, so dir == 0,1");
#endif
  }

  else if (dim == 3){
#if ORDERING == 1
  /* mpi_abort("ORDERING == 1"); */
  if (dir == 0)
    linalg_kron_VEC_TRANSoIoIx_SQR(out, &slicer_1d[side * (deg + 1)], in,
                                   deg + 1);
  else if (dir == 1)
    linalg_kron_IoVEC_TRANSoIx_SQR(out, &slicer_1d[side * (deg + 1)], in,
                                   deg + 1);
  else if (dir == 2)
    linalg_kron_IoIoVEC_TRANSx_SQR(out, &slicer_1d[side * (deg + 1)], in,
                                   deg + 1);
  else
    mpi_abort("DIM == 3 so DIR=0,1 or 2");
#endif
  
#if ORDERING == 2
  mpi_abort("ORDERING == 2");
  /* else if (dim == 3 && dir == 0) */
  /*   linalg_kron_IoVEC_TRANSoIx_SQR(out, &slicer_1d[side * (deg + 1)], in, */
  /*                                  deg + 1); */
  /* else if (dim == 3 && dir == 1) */
  /*   linalg_kron_IoIoVEC_TRANSx_SQR(out, &slicer_1d[side * (deg + 1)], in, */
  /*                                  deg + 1); */
  /* else if (dim == 3 && dir == 2) */
  /*   linalg_kron_VEC_TRANSoIoIx_SQR(out, &slicer_1d[side * (deg + 1)], in, */
  /*                                  deg + 1); */
#endif
  
#if ORDERING == 3
  if (dir == 0)
    linalg_kron_IoIoVEC_TRANSx_SQR(out, &slicer_1d[side * (deg + 1)], in, deg + 1);
  else if (dir == 1)
    linalg_kron_IoVEC_TRANSoIx_SQR(out, &slicer_1d[side * (deg + 1)], in, deg + 1);
  else if (dir == 2)
    linalg_kron_VEC_TRANSoIoIx_SQR(out, &slicer_1d[side * (deg + 1)], in,
                                   deg + 1);
  else
    mpi_abort("DIM == 3 so DIR=0,1 or 2");
#endif
  
  }
  else {
    mpi_abort("DIM must be 2 or 3");
  }
}

void d4est_operators_get_normal(int face, int dim, double* n) {
  mpi_assert((face < 2 * dim) && (dim == 2 || dim == 3));

  double* nx = &n[0];
  double* ny = &n[1];
  double* nz = NULL;

  *nx = 0.;
  *ny = 0.;

  if (dim == 3) {
    nz = &n[2];
    *nz = 0.;
  }
  if (face == 0) {
    *nx = -1.0;
  } else if (face == 1) {
    *nx = 1.0;
  } else if (face == 2) {
    *ny = -1.;
  } else if (face == 3) {
    *ny = 1.;
  } else if (face == 4) {
    *nz = -1.;
  } else if (face == 5) {
    *nz = 1.;
  } else {
    mpi_abort("ERROR: d4est_operators_get_normal");
  }
}

double d4est_operators_jacobi(double r, double alpha, double beta, int N) {
  int i;
  double J0, J1, J2 = 0;
  double h1, anew, bnew;
  double gamma0 = (pow(2.0, alpha + beta + 1) / (alpha + beta + 1.0)) *
                  tgamma(alpha + 1.0) * tgamma(beta + 1.0) /
                  tgamma(alpha + beta + 1.0);
  double gamma1 =
      (alpha + 1.0) * ((beta + 1.0) / (alpha + beta + 3.0)) * gamma0;
  double aold = (2.0 / (2.0 + alpha + beta)) *
                sqrt((alpha + 1.0) * (beta + 1.0) / (alpha + beta + 3.0));

  J0 = 1.0 / sqrt(gamma0);
  if (N == 0) return J0;

  J1 = ((alpha + beta + 2.0) * r / 2.0 + (alpha - beta) / 2) / sqrt(gamma1);
  if (N == 1) return J1;

  // Calculate new polynomials from Jacobi polynomial recurrence relation
  for (i = 0; i < N - 1; i++) {
    h1 = 2.0 * (i + 1) + alpha + beta;
    anew = (2.0 / (h1 + 2.0)) *
           sqrt(((i + 1) + 1.0) * ((i + 1) + 1.0 + alpha + beta) *
                ((i + 1) + 1.0 + alpha) *
                (((i + 1) + 1.0 + beta) / (h1 + 1.0)) / (h1 + 3.0));
    bnew = -((alpha * alpha - beta * beta) / h1) / (h1 + 2.0);
    J2 = (1.0 / anew) * (-aold * J0 + (r - bnew) * J1);

    /* cycle */
    aold = anew;
    J0 = J1;
    J1 = J2;
  }
  return J2;
}

double d4est_operators_gradjacobi(double r, double alpha, double beta, int N) {
  if (N == 0)
    return 0;
  else
    return sqrt(N * (N + alpha + beta + 1.0)) *
           d4est_operators_jacobi(r, alpha + 1.0, beta + 1.0, N - 1);
}

void d4est_operators_vandermonde_1d(double* v1d, double* lobatto_nodes, int degree) {
  int i, j, rows, cols;
  rows = cols = degree + 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      v1d[i * cols + j] = d4est_operators_jacobi(lobatto_nodes[i], 0., 0., j);
}

void d4est_operators_grad_vandermonde_1d(double* dr_v1d, double* lobatto_nodes,
                                int degree) {
  int i, j, rows, cols;
  rows = cols = degree + 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      dr_v1d[i * cols + j] = d4est_operators_gradjacobi(lobatto_nodes[i], 0., 0., j);
}

void d4est_operators_convert_nodal_to_modal(d4est_operators_t* d4est_ops, double* in,
                                   int dim, int deg, double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int nodes = deg + 1;
  double* invVij_1d = P4EST_ALLOC(double, nodes* nodes);

  /* TODO: probably could use build function for inv_Vij */
  d4est_operators_build_Vij_1d(d4est_ops, invVij_1d, deg);
  linalg_invert(invVij_1d, nodes);

  if (dim == 1)
    linalg_matvec_plus_vec(1.0, invVij_1d, in, 0., out, nodes, nodes);
  else if (dim == 2) {
    /* linalg_kron_A1A2x_SQR(out, invVij_1d, invVij_1d, in, nodes, nodes); */
    linalg_kron_A1A2x_NONSQR(out, invVij_1d, invVij_1d, in, nodes, nodes, nodes,
                             nodes);
  } else if (dim == 3) {
    /* linalg_kron_A1A2A3x_SQR(out, invVij_1d, invVij_1d, invVij_1d, in, nodes,
     * nodes, nodes); */
    linalg_kron_A1A2A3x_NONSQR(out, invVij_1d, invVij_1d, invVij_1d, in, nodes,
                               nodes, nodes, nodes, nodes, nodes);
  } else {
    mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension.");
  }

  P4EST_FREE(invVij_1d);
}

double d4est_operators_rtox(double r, double xl, double h) {
  return h * (r + 1.) / 2. + xl;
}

void d4est_operators_rtox_array(double* r, double xl, double h, double* x, int nodes) {
  int i;
  for (i = 0; i < nodes; i++) {
    x[i] = d4est_operators_rtox(r[i], xl, h);
  }
}

static void d4est_operators_build_p_prolong_transpose_1d(d4est_operators_t* d4est_ops,
                                                double* p_prolong_transpose_1d,
                                                int degH, int degh) {
  double* p_prolong_1d = d4est_operators_fetch_p_prolong_1d(d4est_ops, degH, degh);
  linalg_mat_transpose_nonsqr(p_prolong_1d, p_prolong_transpose_1d, (degh + 1),
                              (degH + 1));
}





static void d4est_operators_build_hp_prolong_transpose_1d(
    d4est_operators_t* d4est_ops, double* hp_prolong_transpose_1d, int degH,
    int degh) {
  double* hp_prolong_1d = d4est_operators_fetch_hp_prolong_1d(d4est_ops, degH, degh);
  /* linalg_mat_transpose_nonsqr(hp_prolong_1d, hp_prolong_transpose_1d,
   * 2*(degh+1), (degH+1)); */
  linalg_mat_transpose_nonsqr(hp_prolong_1d, hp_prolong_transpose_1d,
                              (degh + 1), (degH + 1));
  linalg_mat_transpose_nonsqr(&hp_prolong_1d[(degh + 1) * (degH + 1)],
                              &hp_prolong_transpose_1d[(degh + 1) * (degH + 1)],
                              (degh + 1), (degH + 1));
}

static double* d4est_operators_fetch_hp_prolong_transpose_1d(d4est_operators_t* d4est_ops,
                                                    int degH, int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < d4est_ops->d4est_operators_max_storage);
  if (d4est_ops->d4est_operators_hp_prolong_transpose_1d_table[d4est_ops->d4est_operators_max_storage *degH + degh] != NULL) {
    return d4est_ops->d4est_operators_hp_prolong_transpose_1d_table[d4est_ops->d4est_operators_max_storage *degH + degh];
  }
  else {
    int size = 2 * (degh + 1) * (degH + 1);
    d4est_ops->d4est_operators_hp_prolong_transpose_1d_table[d4est_ops->d4est_operators_max_storage *degH + degh] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_hp_prolong_transpose_1d_table[d4est_ops->d4est_operators_max_storage *degH + degh];
    d4est_operators_build_hp_prolong_transpose_1d(d4est_ops, op, degH, degh);
    return op;
  }
}





void d4est_operators_build_p_prolong_1d_inverse(d4est_operators_t* d4est_ops,
                                      double* p_prolong_1d_inverse, int degH,
                                      int degh)
{
  double* p_prolong_1d = d4est_operators_fetch_p_prolong_1d(d4est_ops, degH, degh);
  linalg_leftinverse(p_prolong_1d, p_prolong_1d_inverse, degh + 1, degH + 1);
}





double* d4est_operators_fetch_p_prolong_transpose_1d(d4est_operators_t* d4est_ops,
                                                   int degH, int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_p_prolong_transpose_1d_table[d4est_ops->d4est_operators_max_storage *degH + degh] != NULL) {
    return d4est_ops->d4est_operators_p_prolong_transpose_1d_table[d4est_ops->d4est_operators_max_storage *degH + degh];
  }
  else {
    int size = (degh + 1) * (degH + 1);
    d4est_ops->d4est_operators_p_prolong_transpose_1d_table[d4est_ops->d4est_operators_max_storage *degH + degh] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_p_prolong_transpose_1d_table[d4est_ops->d4est_operators_max_storage *degH + degh];
    d4est_operators_build_p_prolong_transpose_1d(d4est_ops, op, degH, degh);
    return op;
  }
}


void d4est_operators_build_p_prolong_transpose_1d_inverse(d4est_operators_t* d4est_ops,
                                                double* p_prolong_transpose_1d_inverse,
                                                int degH, int degh) {
  double* p_prolong_transpose_1d = d4est_operators_fetch_p_prolong_transpose_1d(d4est_ops, degH, degh);
  linalg_leftinverse(p_prolong_transpose_1d, p_prolong_transpose_1d_inverse, degH + 1, degh + 1);
}

 double* d4est_operators_fetch_p_prolong_transpose_1d_inverse(d4est_operators_t* d4est_ops,
                                                   int degH, int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_p_prolong_transpose_1d_inverse_table[d4est_ops->d4est_operators_max_storage *degH + degh] != NULL) {
    return d4est_ops->d4est_operators_p_prolong_transpose_1d_inverse_table[d4est_ops->d4est_operators_max_storage *degH + degh];
  }
  else {
    int size = (degh + 1) * (degH + 1);
    d4est_ops->d4est_operators_p_prolong_transpose_1d_inverse_table[d4est_ops->d4est_operators_max_storage *degH + degh] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_p_prolong_transpose_1d_inverse_table[d4est_ops->d4est_operators_max_storage *degH + degh];
    d4est_operators_build_p_prolong_transpose_1d_inverse(d4est_ops, op, degH, degh);
    return op;
  }
}


void d4est_operators_apply_hp_prolong_transpose(d4est_operators_t* d4est_ops, double* in,
                                       int* degh, int dim, int degH,
                                       double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int children = (1 << dim);
  int c;
  int nodesH = d4est_operators_get_nodes(dim, degH);
  double* hp_prolong_transpose_1d;

  int stride = 0;
  double* tmp = P4EST_ALLOC(double, nodesH);
  linalg_fill_vec(out, 0., nodesH);

  for (c = 0; c < children; c++) {
    hp_prolong_transpose_1d =
        d4est_operators_fetch_hp_prolong_transpose_1d(d4est_ops, degH, degh[c]);
    /* int nodesh = d4est_operators_get_nodes(dim,degh[c]); */
    d4est_operators_hp_apply_nd_prolong_transpose_with_ptr(
        tmp, degH, &in[stride], degh[c], dim, c, hp_prolong_transpose_1d);
    /* printf("child = %d\n",c); */
    /* util_print_matrix(tmp, nodesH, 1, "tmp = ", 0); */
    /* util_print_matrix(&uh[stride], nodesh, 1, "uh[stride] = ", 0); */
    linalg_vec_axpy(1.0, tmp, out, nodesH);
    /* util_print_matrix(uH, nodesH, 1, "uH = ", 0); */
    stride += d4est_operators_get_nodes(dim, degh[c]);
  }

  P4EST_FREE(tmp);
}

void d4est_operators_apply_p_prolong_transpose(d4est_operators_t* d4est_ops, double* in,
                                      int degh, int dim, int degH,
                                      double* out) {
  double* p_prolong_transpose_1d =
      d4est_operators_fetch_p_prolong_transpose_1d(d4est_ops, degH, degh);

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  /* prolong matrix is identity if degrees are the same */
  if (degh == degH) {
    linalg_copy_1st_to_2nd(in, out, d4est_operators_get_nodes(dim, degh));
    return;
  }

  if (dim == 1)
    linalg_matvec_plus_vec(1.0, p_prolong_transpose_1d, in, 0., out, nodesH,
                           nodesh);
  else if (dim == 2) {
    linalg_kron_A1A2x_NONSQR(out, p_prolong_transpose_1d,
                             p_prolong_transpose_1d, in, nodesH, nodesh, nodesH,
                             nodesh);
  } else if (dim == 3) {
    linalg_kron_A1A2A3x_NONSQR(out, p_prolong_transpose_1d,
                               p_prolong_transpose_1d, p_prolong_transpose_1d,
                               in, nodesH, nodesh, nodesH, nodesh, nodesH,
                               nodesh);
  } else {
    mpi_abort("ERROR: d4est_operators_apply_p_prolong_transpose\n");
  }
}

void d4est_operators_project_side_onto_mortar_space(d4est_operators_t* d4est_ops,
                                           double* in_side, int faces_side,
                                           int* deg_side, double* out_mortar,
                                           int faces_mortar, int* deg_mortar) {
  /* most-common so first */
  /* p-prolong to mortar space */
  if (faces_side == 1 && faces_mortar == 1) {
    d4est_operators_apply_p_prolong(d4est_ops, in_side, deg_side[0], (P4EST_DIM)-1,
                           deg_mortar[0], out_mortar);
  }

  /* hp-prolong to mortar space. */
  else if (faces_side == 1 && faces_mortar == (P4EST_HALF)) {
    d4est_operators_apply_hp_prolong(d4est_ops, in_side, deg_side[0], (P4EST_DIM)-1,
                            &deg_mortar[0], out_mortar);
  }

  /* if faces_side = (P4EST_HALF), then faces_mortar == (P4EST_HALF) */
  /* here we loop through faces and p-prolong */
  else if (faces_side == (P4EST_HALF) && faces_mortar == (P4EST_HALF)) {
    int i;
    int stride_side = 0;
    int stride_mortar = 0;
    for (i = 0; i < faces_side; i++) {
      d4est_operators_apply_p_prolong(d4est_ops, &in_side[stride_side], deg_side[i],
                             (P4EST_DIM)-1, deg_mortar[i],
                             &out_mortar[stride_mortar]);

      stride_side += d4est_operators_get_nodes((P4EST_DIM)-1, deg_side[i]);
      stride_mortar += d4est_operators_get_nodes((P4EST_DIM)-1, deg_mortar[i]);
    }
  } else {
    mpi_abort("ERROR: d4est_operators_project_side_onto_mortar_space");
  }
}

void d4est_operators_project_mass_mortar_onto_side(d4est_operators_t* dgmath,
                                          double* in_mortar, int faces_mortar,
                                          int* deg_mortar, double* out_side,
                                          int faces_side, int* deg_side) {

  /* most-common so first */
  /* p-prolong to mortar space */
  if (faces_side == 1 && faces_mortar == 1) {
    d4est_operators_apply_p_prolong_transpose(dgmath, in_mortar, deg_mortar[0], (P4EST_DIM)-1,
                            deg_side[0], out_side);
  }

  /* hp-prolong to mortar space. */
  else if (faces_side == 1 && faces_mortar == (P4EST_HALF)) {
    d4est_operators_apply_hp_prolong_transpose(dgmath, in_mortar, &deg_mortar[0], (P4EST_DIM)-1,
                             deg_side[0], out_side);
  }

  /* if faces_side = (P4EST_HALF), then faces_mortar == (P4EST_HALF) */
  /* here we loop through faces and p-prolong */
  else if (faces_side == (P4EST_HALF) && faces_mortar == (P4EST_HALF)) {
    int i;
    int stride_side = 0;
    int stride_mortar = 0;
    for (i = 0; i < faces_side; i++) {
      d4est_operators_apply_p_prolong_transpose(dgmath, &in_mortar[stride_mortar], deg_mortar[i],
                              (P4EST_DIM)-1, deg_side[i],
                              &out_side[stride_side]);

      stride_side += d4est_operators_get_nodes((P4EST_DIM)-1, deg_side[i]);
      stride_mortar += d4est_operators_get_nodes((P4EST_DIM)-1, deg_mortar[i]);
    }
  }

  else {
    mpi_abort("ERROR: d4est_operators_project_mass_mortar_onto_side_space");
  }
}


void d4est_operators_project_mortar_onto_side(d4est_operators_t* d4est_ops,
                                     double* in_mortar, int faces_mortar,
                                     int* deg_mortar, double* out_side,
                                     int faces_side, int* deg_side) {
  /* most-common so first */
  /* p-prolong to mortar space */
  if (faces_side == 1 && faces_mortar == 1) {
    d4est_operators_apply_p_restrict(d4est_ops, in_mortar, deg_mortar[0], (P4EST_DIM)-1,
                            deg_side[0], out_side);
  }

  /* hp-prolong to mortar space. */
  else if (faces_side == 1 && faces_mortar == (P4EST_HALF)) {
    d4est_operators_apply_hp_restrict(d4est_ops, in_mortar, &deg_mortar[0], (P4EST_DIM)-1,
                             deg_side[0], out_side);
  }

  /* if faces_side = (P4EST_HALF), then faces_mortar == (P4EST_HALF) */
  /* here we loop through faces and p-prolong */
  else if (faces_side == (P4EST_HALF) && faces_mortar == (P4EST_HALF)) {
    int i;
    int stride_side = 0;
    int stride_mortar = 0;
    for (i = 0; i < faces_side; i++) {
      d4est_operators_apply_p_restrict(d4est_ops, &in_mortar[stride_mortar], deg_mortar[i],
                              (P4EST_DIM)-1, deg_side[i],
                              &out_side[stride_side]);

      stride_side += d4est_operators_get_nodes((P4EST_DIM)-1, deg_side[i]);
      stride_mortar += d4est_operators_get_nodes((P4EST_DIM)-1, deg_mortar[i]);
    }
  }

  else {
    mpi_abort("ERROR: d4est_operators_project_side_onto_mortar_space");
  }
}

static
void d4est_operators_build_hp_restrict_interp_1d_aux(
    int degh, int degH, int c, double* inv_v1d_trans_degh,
    double* lobatto_nodes_degH, double* hp_restrict_interp_matrix_1d) {
  mpi_assert(degH <= degh);
  mpi_assert(c == 0 || c == 1);

  int n, i;
  int nodesH = degH + 1;
  int nodesh = degh + 1;
  double* phi_modal = (double*)P4EST_ALLOC(double, nodesh);
  double r;

  for (n = 0; n < nodesH; n++) {
    r = lobatto_nodes_degH[n];

    if (c == 0)
      r = 2 * r + 1;
    else
      r = 2 * r - 1;

    for (i = 0; i < nodesh; i++) {
      phi_modal[i] = d4est_operators_jacobi(r, 0, 0, i);
    }

    if (n == degH / 2 && degH % 2 == 0)
      linalg_matvec_plus_vec(.5, inv_v1d_trans_degh, phi_modal, 0.0,
                             &hp_restrict_interp_matrix_1d[n * nodesh], nodesh,
                             nodesh);
    else if (r < -1. || r > 1.)
      linalg_fill_vec(&hp_restrict_interp_matrix_1d[n * nodesh], 0., nodesh);
    else
      linalg_matvec_plus_vec(1., inv_v1d_trans_degh, phi_modal, 0.0,
                             &hp_restrict_interp_matrix_1d[n * nodesh], nodesh,
                             nodesh);
  }

  P4EST_FREE(phi_modal);
}

static void d4est_operators_build_hp_restrict_interp_1d(d4est_operators_t* d4est_ops,
                                               double* hp_restrict_interp_1d,
                                               int degH, int degh) {
  int edge_nodes_degH = degH + 1;
  int edge_nodes_degh = degh + 1;

  double* v1d = (double*)P4EST_ALLOC(double, edge_nodes_degh* edge_nodes_degh);
  double* inv_v1d =
      (double*)P4EST_ALLOC(double, edge_nodes_degh* edge_nodes_degh);
  double* inv_v1d_trans_degh =
      (double*)P4EST_ALLOC(double, edge_nodes_degh* edge_nodes_degh);

  d4est_operators_build_Vij_1d(d4est_ops, v1d, degh);
  linalg_invert_and_copy(v1d, inv_v1d, edge_nodes_degh);
  linalg_mat_transpose(inv_v1d, inv_v1d_trans_degh, edge_nodes_degh);

  double* lobatto_nodes_degH = d4est_operators_fetch_GLL_nodes_1d(d4est_ops, degH);
  d4est_operators_build_hp_restrict_interp_1d_aux(degh, degH, 0, inv_v1d_trans_degh,
                                         lobatto_nodes_degH,
                                         &hp_restrict_interp_1d[0]);
  d4est_operators_build_hp_restrict_interp_1d_aux(
      degh, degH, 1, inv_v1d_trans_degh, lobatto_nodes_degH,
      &hp_restrict_interp_1d[edge_nodes_degH * edge_nodes_degh]);

  /* util_print_matrix(hp_restrict_interp_1d, edge_nodes_degH, edge_nodes_degh,
   * "hp_restrict_interp_1d_LOWER = ", 0); */

  /* util_print_matrix(&hp_restrict_interp_1d[edge_nodes_degH*edge_nodes_degh],
   * edge_nodes_degH, edge_nodes_degh, "hp_restrict_interp_1d_UPPER = ", 0); */

  P4EST_FREE(inv_v1d);
  P4EST_FREE(inv_v1d_trans_degh);
  P4EST_FREE(v1d);
}

static double* d4est_operators_fetch_hp_restrict_interp_1d(d4est_operators_t* d4est_ops,
                                                  int degH, int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_hp_restrict_interp_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh] != NULL) {
    return d4est_ops->d4est_operators_hp_restrict_interp_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh];
  }
  else {
    int size = 2 * (degh + 1) * (degH + 1);
    d4est_ops->d4est_operators_hp_restrict_interp_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_hp_restrict_interp_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh];
    d4est_operators_build_hp_restrict_interp_1d(d4est_ops, op, degH, degh);
    return op;
  }
}

void d4est_operators_apply_hp_restrict_interp(d4est_operators_t* d4est_ops, double* in,
                                     int* degh, int dim, int degH,
                                     double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int children = (1 << dim);
  int c;
  int nodesH = d4est_operators_get_nodes(dim, degH);
  double* hp_restrict_1d;

  int stride = 0;
  double* tmp = P4EST_ALLOC(double, nodesH);
  linalg_fill_vec(out, 0., nodesH);

  for (c = 0; c < children; c++) {
    hp_restrict_1d = d4est_operators_fetch_hp_restrict_interp_1d(d4est_ops, degH, degh[c]);
    /* int nodesh = d4est_operators_get_nodes(dim,degh[c]); */

    d4est_operators_hp_apply_nd_restrict_with_ptr(tmp, degH, &in[stride], degh[c], dim,
                                         c, hp_restrict_1d);
    /* printf("child = %d\n",c); */
    /* util_print_matrix(tmp, nodesH, 1, "tmp = ", 0); */
    /* util_print_matrix(&uh[stride], nodesh, 1, "uh[stride] = ", 0); */
    linalg_vec_axpy(1.0, tmp, out, nodesH);
    /* util_print_matrix(uH, nodesH, 1, "uH = ", 0); */
    stride += d4est_operators_get_nodes(dim, degh[c]);
  }

  P4EST_FREE(tmp);
}

static void d4est_operators_build_p_restrict_interp_1d_aux(
    int degh, int degH, double* inv_v1d_trans_degh, double* lobatto_nodes_degH,
    double* hp_restrict_interp_matrix_1d) {
  int n, i;
  int nodesH = degH + 1;
  int nodesh = degh + 1;
  double* phi_modal = (double*)P4EST_ALLOC(double, nodesh);
  double r;

  for (n = 0; n < nodesH; n++) {
    r = lobatto_nodes_degH[n];
    for (i = 0; i < nodesh; i++) {
      phi_modal[i] = d4est_operators_jacobi(r, 0, 0, i);
    }
    linalg_matvec_plus_vec(1., inv_v1d_trans_degh, phi_modal, 0.0,
                           &hp_restrict_interp_matrix_1d[n * nodesh], nodesh,
                           nodesh);
  }

  P4EST_FREE(phi_modal);
}

static void d4est_operators_build_p_restrict_interp_1d(d4est_operators_t* d4est_ops,
                                              double* hp_restrict_interp_1d,
                                              int degH, int degh) {
  mpi_assert(degH <= degh);
  /* int edge_nodes_degH = degH + 1; */
  int edge_nodes_degh = degh + 1;

  double* v1d = (double*)P4EST_ALLOC(double, edge_nodes_degh* edge_nodes_degh);
  double* inv_v1d =
      (double*)P4EST_ALLOC(double, edge_nodes_degh* edge_nodes_degh);
  double* inv_v1d_trans_degh =
      (double*)P4EST_ALLOC(double, edge_nodes_degh* edge_nodes_degh);

  d4est_operators_build_Vij_1d(d4est_ops, v1d, degh);
  linalg_invert_and_copy(v1d, inv_v1d, edge_nodes_degh);
  linalg_mat_transpose(inv_v1d, inv_v1d_trans_degh, edge_nodes_degh);

  double* lobatto_nodes_degH = d4est_operators_fetch_GLL_nodes_1d(d4est_ops, degH);
  d4est_operators_build_p_restrict_interp_1d_aux(degh, degH, inv_v1d_trans_degh,
                                        lobatto_nodes_degH,
                                        &hp_restrict_interp_1d[0]);

  P4EST_FREE(inv_v1d);
  P4EST_FREE(inv_v1d_trans_degh);
  P4EST_FREE(v1d);
}

static double* d4est_operators_fetch_p_restrict_interp_1d(d4est_operators_t* d4est_ops,
                                                 int degH, int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < d4est_ops->d4est_operators_max_storage);

  if (d4est_ops->d4est_operators_p_restrict_interp_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh] != NULL) {
    return d4est_ops->d4est_operators_p_restrict_interp_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh];
  }
  else {
    int size = (degh + 1) * (degH + 1);
    d4est_ops->d4est_operators_p_restrict_interp_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_p_restrict_interp_1d_table[d4est_ops->d4est_operators_max_storage * degH + degh];
    d4est_operators_build_p_restrict_interp_1d(d4est_ops, op, degH, degh);
    return op;
  }
}

void d4est_operators_apply_p_restrict_interp(d4est_operators_t* d4est_ops, double* in,
                                    int degh, int dim, int degH, double* out) {
  double* p_restrict_1d = d4est_operators_fetch_p_restrict_interp_1d(d4est_ops, degH, degh);

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  /* restrict matrix should be the identity if degrees are the same */
  if (degh == degH) {
    linalg_copy_1st_to_2nd(in, out, d4est_operators_get_nodes(dim, degh));
    return;
  }

  if (dim == 1)
    linalg_matvec_plus_vec(1.0, p_restrict_1d, in, 0., out, nodesH, nodesh);
  else if (dim == 2) {
    linalg_kron_A1A2x_NONSQR(out, p_restrict_1d, p_restrict_1d, in, nodesH,
                             nodesh, nodesH, nodesh);
  } else if (dim == 3) {
    linalg_kron_A1A2A3x_NONSQR(out, p_restrict_1d, p_restrict_1d, p_restrict_1d,
                               in, nodesH, nodesh, nodesH, nodesh, nodesH,
                               nodesh);
  } else {
    mpi_abort("ERROR: d4est_operators_restrict_nd_U\n");
  }
}

//         +----------------------+
//         |			  |
//         |			  |	     b
//         |			  |	     |
//         |			  |	     |
//         |	   face f	  |	     |
//         |	       	       	  |  	     |
//         |			  |	     /---------	a
//         |			  |	    /
//         |			  |	   /
//         |			  |	  c
//         |----------------------+



void d4est_operators_get_face_info(int f, d4est_operators_face_info_t* face_info) {

  if (f == 0) {
    face_info->a = 1;
    face_info->b = 2;
    face_info->sgn = -1.;
    face_info->c = 0;
  } else if (f == 1) {
    face_info->a = 1;
    face_info->b = 2;
    face_info->sgn = 1.;
    face_info->c = 0;
  } else if (f == 2) {
    face_info->a = 0;
    face_info->b = 2;
    face_info->sgn = -1.;
    face_info->c = 1;
  } else if (f == 3) {
    face_info->a = 0;
    face_info->b = 2;
    face_info->sgn = 1.;
    face_info->c = 1;
  } else if (f == 4) {
    face_info->a = 0;
    face_info->b = 1;
    face_info->sgn = -1.;
    face_info->c = 2;
  } else if (f==5) {
    face_info->a = 0;
    face_info->b = 1;
    face_info->sgn = 1.;
    face_info->c = 2;
  }
  else {
    mpi_abort("face info error: only f > 5 or f < 0 in 3d");
  }

}

/* TODO: OPTIMIZE BY PUTTING IN A degh = degH check */
void d4est_operators_compute_M_vec_wo_aliasing(d4est_operators_t* d4est_ops,
                                      double* xyz_rst[(P4EST_DIM)][P4EST_DIM],
                                      double* vec, int degH, double* MJvec,
                                     int degh) {
#if (P4EST_DIM)==3
  int volume_nodes_degH = d4est_operators_get_nodes((P4EST_DIM), degH);
  int volume_nodes_degh = d4est_operators_get_nodes((P4EST_DIM), degh);

  double* xyz_rst_prolonged[(P4EST_DIM)][P4EST_DIM];
  double* vec_prolonged = P4EST_ALLOC(double, volume_nodes_degh);
  double* Jvec_prolonged = P4EST_ALLOC(double, volume_nodes_degh);
  double* Jvec = P4EST_ALLOC(double, volume_nodes_degH);

  int i, j;
  for (i = 0; i < (P4EST_DIM); i++) {
    for (j = 0; j < (P4EST_DIM); j++) {
      xyz_rst_prolonged[i][j] = P4EST_ALLOC(double, volume_nodes_degh);
      d4est_operators_apply_p_prolong(d4est_ops, &(xyz_rst[i][j][0]), degH,
                             (P4EST_DIM), degh, &(xyz_rst_prolonged[i][j][0]));
    }
  }

  d4est_operators_apply_p_prolong(d4est_ops, vec, degH, (P4EST_DIM), degh,
                         vec_prolonged);

  for (i = 0; i < volume_nodes_degh; i++) {
    double xr = xyz_rst_prolonged[0][0][i];
    double xs = xyz_rst_prolonged[0][1][i];
    double xt = xyz_rst_prolonged[0][2][i];

    double yr = xyz_rst_prolonged[1][0][i];
    double ys = xyz_rst_prolonged[1][1][i];
    double yt = xyz_rst_prolonged[1][2][i];

    double zr = xyz_rst_prolonged[2][0][i];
    double zs = xyz_rst_prolonged[2][1][i];
    double zt = xyz_rst_prolonged[2][2][i];

    double J = xr * (ys * zt - zs * yt) - yr * (xs * zt - zs * xt) +
               zr * (xs * yt - ys * xt);

    Jvec_prolonged[i] = J * vec_prolonged[i];
  }

  d4est_operators_apply_p_restrict(d4est_ops, Jvec_prolonged, degh, (P4EST_DIM),
                          degH, Jvec);

  d4est_operators_apply_Mij(d4est_ops, Jvec, (P4EST_DIM), degH, MJvec);

  for (i = 0; i < (P4EST_DIM); i++) {
    for (j = 0; j < (P4EST_DIM); j++) {
      P4EST_FREE(xyz_rst_prolonged[i][j]);
    }
  }
  
  P4EST_FREE(vec_prolonged);
  P4EST_FREE(Jvec_prolonged);
  P4EST_FREE(Jvec);
#else
  mpi_abort("d4est_operators_compute_M_vec_wo_aliasing currently only supports DIM=3");
#endif
}


void d4est_operators_compute_invM_vec_wo_aliasing(d4est_operators_t* d4est_ops,
                                         double* xyz_rst[(P4EST_DIM)][P4EST_DIM],
                                         double* vec, int degH, double* invMvec,
                                         int degh) {
#if (P4EST_DIM)==3 
  /* int volume_nodes_degH = d4est_operators_get_nodes((P4EST_DIM), degH); */
  int volume_nodes_degh = d4est_operators_get_nodes((P4EST_DIM), degh);

  double* xyz_rst_prolonged[(P4EST_DIM)][P4EST_DIM];
  double* invMvec_prolonged = P4EST_ALLOC(double, volume_nodes_degh);
  double* invMvec_divJ_prolonged = P4EST_ALLOC(double, volume_nodes_degh);

  int i, j;
  for (i = 0; i < (P4EST_DIM); i++) {
    for (j = 0; j < (P4EST_DIM); j++) {
      xyz_rst_prolonged[i][j] = P4EST_ALLOC(double, volume_nodes_degh);
      d4est_operators_apply_p_prolong(d4est_ops, &(xyz_rst[i][j][0]), degH,
                             (P4EST_DIM), degh, &(xyz_rst_prolonged[i][j][0]));
    }
  }

  d4est_operators_apply_invMij(d4est_ops, vec, (P4EST_DIM), degH, invMvec);
  d4est_operators_apply_p_prolong(d4est_ops, invMvec, degH, (P4EST_DIM), degh,
                         invMvec_prolonged);

 
  
  for (i = 0; i < volume_nodes_degh; i++) {
    double xr = xyz_rst_prolonged[0][0][i];
    double xs = xyz_rst_prolonged[0][1][i];
    double xt = xyz_rst_prolonged[0][2][i];

    double yr = xyz_rst_prolonged[1][0][i];
    double ys = xyz_rst_prolonged[1][1][i];
    double yt = xyz_rst_prolonged[1][2][i];

    double zr = xyz_rst_prolonged[2][0][i];
    double zs = xyz_rst_prolonged[2][1][i];
    double zt = xyz_rst_prolonged[2][2][i];

    double J = xr * (ys * zt - zs * yt) - yr * (xs * zt - zs * xt) +
               zr * (xs * yt - ys * xt);

    invMvec_divJ_prolonged[i] = invMvec_prolonged[i]/J;
  }

  d4est_operators_apply_p_restrict(d4est_ops, invMvec_divJ_prolonged, degh, (P4EST_DIM),
                          degH, invMvec);

  for (i = 0; i < (P4EST_DIM); i++) {
    for (j = 0; j < (P4EST_DIM); j++) {
      P4EST_FREE(xyz_rst_prolonged[i][j]);
    }
  }
   
  P4EST_FREE(invMvec_divJ_prolonged);
  P4EST_FREE(invMvec_prolonged);
#else
  mpi_abort(" d4est_operators_compute_invM_vec_wo_aliasing currently only supports DIM=3 ");
#endif
}

static void d4est_operators_build_FLIP_1d(double* FLIP_1d, int deg) {
  memset(FLIP_1d, 0., sizeof(double) * (deg + 1) * (deg + 1));
  int offset;
  int i;
  for (i = 0; i < deg + 1; i++) {
    offset = (i+1)*deg;
    FLIP_1d[offset] = 1.;
  }
}

static double* d4est_operators_fetch_FLIP_1d(d4est_operators_t* d4est_ops, int deg) {
  mpi_assert(deg < d4est_ops->d4est_operators_max_storage);
  if (d4est_ops->d4est_operators_FLIP_1d_table[deg] != NULL) {
    return d4est_ops->d4est_operators_FLIP_1d_table[deg];
  }
  else {
    int size = (deg + 1) * (deg + 1);
    d4est_ops->d4est_operators_FLIP_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_FLIP_1d_table[deg];
    d4est_operators_build_FLIP_1d(op, deg);
    return op;
  }
}

void d4est_operators_apply_FLIP(d4est_operators_t* d4est_ops, double* in, int dim, int deg,
                       int dir, double* out) {
  mpi_assert(dim == 1 || dim == 2);
  double* flip_1d = d4est_operators_fetch_FLIP_1d(d4est_ops, deg);
  int nodes = deg + 1;
  if (dim == 1)
    linalg_matvec_plus_vec(1.0, flip_1d, in, 0., out, nodes, nodes);
  else if (dim == 2) {
#if (ORDERING)==1
    if (dir == 0) linalg_kron_MAToIx_SQR(out, flip_1d, in, nodes);
    if (dir == 1) linalg_kron_IoMATx_SQR(out, flip_1d, in, nodes);
    /* mpi_abort("ORDERING==1"); */
#endif
#if (ORDERING)==3
    if (dir == 0) linalg_kron_IoMATx_SQR(out, flip_1d, in, nodes);
    if (dir == 1) linalg_kron_MAToIx_SQR(out, flip_1d, in, nodes);
#endif
    if (dir == 2)
      linalg_kron_A1A2x_NONSQR(out, flip_1d, flip_1d, in, nodes, nodes, nodes,
                               nodes);
  } else {
    mpi_abort("ERROR: FLIP not supported in this dimension atm.");
  }
}

int d4est_operators_reorient_face_order
(
 int face_dim,
 int f_m,
 int f_p,
 int o,
 int i
)
{
  mpi_assert((face_dim == 1 && o < 2 && i < 2) || (face_dim == 2 && o < 4 && i < 4));

  if (face_dim == 1){
    if (o == 1){
      return (i == 0) ? 1 : 0;
    }
    else {
      return i;
    }
  }
  else if (face_dim == 2){
    int perm = d4est_operators_p8est_code_to_perm[d4est_operators_p8est_FToF_code[f_m][f_p]][o];
    return d4est_operators_p8est_perm_to_order[perm][i];
  }
  else {
    mpi_abort("FACE DIM == 1, 2\n");
    return -1;
  }
}



void d4est_operators_reorient_face_data
(
 d4est_operators_t* d4est_ops,
 double* in,
 int face_dim,
 int deg,
 int o, 
 int f_m,
 int f_p,
 double* out
)
{
  /* dim == 1 and o == 4 hasn't been unit tested therefore check for them */
  mpi_assert((face_dim == 2 && o < 4) || (face_dim == 1 && o < 2));

  if (o == 4)
    printf("[WARNING]: o == 4 and this orientation has not been tested");
  
  int nodes = d4est_operators_get_nodes(face_dim, deg);
  
  if (face_dim == 1){
    if (o == 1)
      d4est_operators_apply_FLIP
        (
         d4est_ops,
         in,
         face_dim,
         deg,
         0,
         out
        );
    else
      linalg_copy_1st_to_2nd(in, out, nodes);
    return;
  }
  double* tmp0 = P4EST_ALLOC(double, nodes);
  double* tmp1 = P4EST_ALLOC(double, nodes);
  
  int ftransform [9];
  p4est_expand_face_transform
    (
     ((f_m <= f_p) ? f_m : f_p), /* dominant face */
     (P4EST_FACES)*o + ((f_m <= f_p) ? f_p : f_m), /* slave face */
     ftransform
    );
  
  int flip0 = ftransform[6];
  int flip1 = ftransform[7];
  int aligned = ((ftransform[1] - ftransform[0])
                *(ftransform[4] - ftransform[3]) > 0);

  if (flip0 == 1){
    d4est_operators_apply_FLIP
      (
       d4est_ops,
       in,
       face_dim,
       deg,
       0,
       tmp0
      );
  }
  else {
    P4EST_FREE(tmp0);
    tmp0 = in;
  }

  if (flip1 == 1){
    d4est_operators_apply_FLIP
      (
       d4est_ops,
       tmp0,
       face_dim,
       deg,
       1,
       tmp1
      );
  }
  else {
    P4EST_FREE(tmp1);
    tmp1 = tmp0;
  }
  
  if (aligned == 0){
    linalg_mat_transpose(tmp1, out, deg+1);
  }
  else {
    linalg_copy_1st_to_2nd(tmp1, out, nodes);
  }

  if (tmp0 != in)
    P4EST_FREE(tmp0);
  if(tmp1 != tmp0)
    P4EST_FREE(tmp1);
}


/* d4est_rst_t */
/* d4est_operators_get_rst_points */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  int deg, */
/*  int dim, */
/*  quadrature_type_t type */
/* ) */
/* { */
/*   d4est_rst_t rst; */
/*   double* rsttmp [3] = {NULL}; */
  
/*   if (type == QUAD_GAUSS){ */
/*     for (int d = 0; d < dim; d++){ */
/*       rsttmp[d] = d4est_operators_fetch_Gauss_xyz_nd */
/*                   ( */
/*                    d4est_ops, */
/*                    dim, */
/*                    deg, */
/*                    d */
/*                   ); */
/*     } */
/*   } */
/*   else if (type == QUAD_LOBATTO) { */
/*     for (int d = 0; d < dim; d++){ */
/*       rsttmp[d] = d4est_operators_fetch_xyz_nd */
/*                   ( */
/*                    d4est_ops, */
/*                    dim, */
/*                    deg, */
/*                    d */
/*                   ); */
/*     } */
/*   } */
/*   else { */
/*     mpi_abort("type == Lobatto or Gauss"); */
/*   } */

/*   rst.r = rsttmp[0]; */
/*   rst.s = rsttmp[1]; */
/*   rst.t = rsttmp[2]; */
  
/*   return rst; */
/* } */


int d4est_operators_corner_to_node
(
 int dim,
 int deg,
 int corner
)
{
  mpi_assert(
             (
              (dim == 3 && corner < 8)
              || (dim == 2 && dim < 4)
             )
            );
  
  int nodes_1d = deg+1;
  int nodes_2d = nodes_1d*nodes_1d;
  int nodes_3d = nodes_1d*nodes_1d*nodes_1d; 

#if (ORDERING == 1)
  if (dim == 2)  {
    if (corner == 0)
      return 0;
    else if (corner == 1)
      return nodes_2d - nodes_1d;
    else if (corner == 2)
      return nodes_1d - 1;
    else if (corner == 3)
      return nodes_2d - 1;
    else{
      mpi_abort("[ERROR]: corner < 4\n");
      return -1;
    }
  }

  else if (dim == 3){
    if (corner == 0)
      return 0;
    else if (corner == 1)
      return nodes_3d - nodes_2d;
    else if (corner == 2)
      return nodes_2d - nodes_1d;
    else if (corner == 3)
      return nodes_3d - nodes_1d;
    else if (corner == 4)
        return nodes_1d - 1;
    else if (corner == 5)
      return nodes_3d - nodes_2d + nodes_1d - 1;
    else if (corner == 6)
      return nodes_2d - 1;
    else if (corner == 7)
      return nodes_3d - 1;
    else{
      mpi_abort("[ERROR]: corner < 8\n");
      return -1;
    }
  }

  else {
    mpi_abort("[ERROR]: DIM == 3 or DIM == 2\n");
    return -1;
  }
#endif
#if (ORDERING == 3)
    if (dim == 2)  {
    if (corner == 0)
      return 0;
    else if (corner == 1) 
      return nodes_1d - 1;
    else if (corner == 2)
      return nodes_2d - nodes_1d;
    else if (corner == 3)
      return nodes_2d - 1;
    else{
      mpi_abort("[ERROR]: corner < 4\n");
      return -1;
    }
  }

  else if (dim == 3){
    if (corner == 0)
      return 0;
    else if (corner == 1)
        return nodes_1d - 1;
    else if (corner == 2) 
      return nodes_2d - nodes_1d;
    else if (corner == 3)
      return nodes_2d - 1;
    else if (corner == 4)
      return nodes_3d - nodes_2d;
    else if (corner == 5)
      return nodes_3d - nodes_2d + nodes_1d - 1;
    else if (corner == 6)
      return nodes_3d - nodes_1d;
    else if (corner == 7)
      return nodes_3d - 1;
    else{
      mpi_abort("[ERROR]: corner < 8\n");
      return -1;
    }
  }

  else {
    mpi_abort("[ERROR]: DIM == 3 or DIM == 2\n");
    return -1;
  }
#endif
}


static void d4est_operators_build_vtk_interp_1d(double* vtk_interp_1d, int deg) {
  memset(vtk_interp_1d, 0., sizeof(double) * (deg + 1) * (deg) * 2);
  /* left corners */
  for (int i = 0; i < deg; i++) {
    for (int j = 0; j < deg + 1; j++) {
      if (i == j)
        vtk_interp_1d[i*(deg+1) + j] = 1.;
    }
  }

  /* right corners */
  int stride = (deg + 1) * deg;
  for (int i = 0; i < deg; i++) {
    for (int j = 0; j < deg + 1; j++) {
      if (i + 1 == j)
        vtk_interp_1d[stride + i*(deg+1) + j] = 1.;
    }
  }
}

static double* d4est_operators_fetch_vtk_interp_1d(d4est_operators_t* d4est_ops, int deg) {
  if (d4est_ops->d4est_operators_vtk_interp_1d_table[deg-1] != NULL) {
    return d4est_ops->d4est_operators_vtk_interp_1d_table[deg-1];
  }
  else {
    int size = (deg) * (deg + 1) * 2;
    d4est_ops->d4est_operators_vtk_interp_1d_table[deg-1] = P4EST_ALLOC(double, size);
    double* op = d4est_ops->d4est_operators_vtk_interp_1d_table[deg-1];
    d4est_operators_build_vtk_interp_1d(op, deg);
    return op;
  }
}

void d4est_operators_apply_vtk_interp
(
 d4est_operators_t* d4est_ops,
 double* in,
 int dim,
 int deg,
 int c,
 double* out
){
  mpi_assert((c < (1 << dim)));

  double* vtk_interp_1d = d4est_operators_fetch_vtk_interp_1d(d4est_ops, deg);
  
  int nodes_in = deg + 1;
  int nodes_out = deg;

  if (dim == 1)
    linalg_matvec_plus_vec(1.0, &vtk_interp_1d[c * nodes_in * nodes_out], in,
                           0., out, nodes_out, nodes_in);
  else if (dim == 2) {
    int cx = d4est_operators_is_child_left_or_right(c, 0);
    int cy = d4est_operators_is_child_left_or_right(c, 1);
    linalg_kron_A1A2x_NONSQR(out, &vtk_interp_1d[cx * nodes_in * nodes_out],
                             &vtk_interp_1d[cy * nodes_in * nodes_out], in,
                             nodes_out, nodes_in, nodes_out, nodes_in);
  } else if (dim == 3) {
    int cx = d4est_operators_is_child_left_or_right(c, 0);
    int cy = d4est_operators_is_child_left_or_right(c, 1);
    int cz = d4est_operators_is_child_left_or_right(c, 2);
    linalg_kron_A1A2A3x_NONSQR(out, &vtk_interp_1d[cx * nodes_in * nodes_out],
                               &vtk_interp_1d[cy * nodes_in * nodes_out],
                               &vtk_interp_1d[cz * nodes_in * nodes_out], in,
                               nodes_out, nodes_in, nodes_out, nodes_in, nodes_out, nodes_in);
  } else {
    mpi_abort("Dim = 1 or 2 or 3 in vtk_interp");
  }
}

static void d4est_operators_build_vtk_rst
(
 d4est_operators_t* d4est_ops,
 double* vtk_rst,
 int dim, int deg
)
{
  int deg_dim = (dim == 2) ? deg*deg : deg*deg*deg;
  double* temp = P4EST_ALLOC(double, deg_dim);
  int children = (dim == 2) ? 4 : 8;
  
  for (int c = 0; c < children; c++){
    for (int dir = 0; dir < dim; dir++){
      double* rst = d4est_operators_fetch_xyz_nd(d4est_ops, dim, deg, dir);
      d4est_operators_apply_vtk_interp(d4est_ops, rst, dim, deg, c, temp);
      for (int n = 0; n < deg_dim; n++){
        vtk_rst[dir + c*3 + n*children*3] = temp[n];
      }
    }
  }  

  /* if DIM is 2, zero out z coordinate */
  if (dim == 2){
    for (int c = 0; c < children; c++){
      for (int n = 0; n < deg_dim; n++){
        vtk_rst[2 + c*3 + n*children*3] = 0.;
      }
    }
  }
  P4EST_FREE(temp);
}

void d4est_operators_convert_nodal_to_vtk
(
 d4est_operators_t* d4est_ops,
 double* vec,
 int dim,
 int deg,
 double* vtk_vec /* should be (chiildren)*(deg)^dim */
)
{
  int deg_dim = (dim == 2) ? deg*deg : deg*deg*deg;
  double* temp = P4EST_ALLOC(double, deg_dim);
  int children = (dim == 2) ? 4 : 8;
  
  for (int c = 0; c < children; c++){
    d4est_operators_apply_vtk_interp(d4est_ops, vec, dim, deg, c, temp);
    for (int n = 0; n < deg_dim; n++){
      vtk_vec[c + n*children] = temp[n];
    }
  }

  P4EST_FREE(temp);
}


double* d4est_operators_fetch_vtk_rst
(
 d4est_operators_t* d4est_ops,
 int deg,
 int dim
)
{
  mpi_assert(deg < d4est_ops->d4est_operators_max_storage);

  double* table_ptr = (dim==3) ? d4est_ops->d4est_operators_vtk_rst_3d_table[deg] : d4est_ops->d4est_operators_vtk_rst_2d_table[deg];
  
  if (table_ptr != NULL) {
    return table_ptr;
  }
  else {
    int children = (dim == 2) ? 4 : 8;
    int deg_dim = (dim == 2) ? deg*deg : deg*deg*deg;
    int size = 3*deg_dim*(children);

    if (dim == 3){
      d4est_ops->d4est_operators_vtk_rst_3d_table[deg] = P4EST_ALLOC(double, size);
    }
    else if (dim == 2){
      d4est_ops->d4est_operators_vtk_rst_2d_table[deg] = P4EST_ALLOC(double, size);
    }
    else{
      mpi_abort("[D4EST_ERROR]: Not a support dimension");
    }
    
    double* op = (dim==3) ? d4est_ops->d4est_operators_vtk_rst_3d_table[deg] : d4est_ops->d4est_operators_vtk_rst_2d_table[deg];
    d4est_operators_build_vtk_rst(d4est_ops, op, dim, deg);
    return op;
  }
}


/* void */
/* d4est_operators_compute_derivative_on_quadrature_points */
/* ( */
/*  double* vec, */
/*  double* rst_xyz_quad [(P4EST_DIM)][(P4EST_DIM)], */
/*  double* dvec [(P4EST_DIM)], */
/*  int deg_lobatto, */
/*  int deg_quad, */
/*  quadrature_type_t quad_type, */
/*  d4est_operators_t* d4est_ops */
/* ) */
/* { */

/*   int volume_nodes_quad = d4est_operators_get_nodes((P4EST_DIM), deg_quad); */
/*   int volume_nodes_lobatto = d4est_operators_get_nodes((P4EST_DIM), deg_lobatto); */

/*   double* dvec_di_prolonged = P4EST_ALLOC(double, volume_nodes_quad); */
/*   double* dvec_di_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*   double* dvec_di_lobatto = P4EST_ALLOC(double, volume_nodes_lobatto); */

/*   for (int j = 0; j < (P4EST_DIM); j++){ */
/*     for (int k = 0; k < volume_nodes_quad; k++){ */
/*       dvec[j][k] = 0.; */
/*     } */
/*   } */
  
/*   for (int i = 0; i < (P4EST_DIM); i++){ */

/*     d4est_operators_apply_Dij(d4est_ops, vec, (P4EST_DIM), deg_lobatto, i, dvec_di_lobatto); */

/*     d4est_operators_interp(d4est_ops, dvec_di_lobatto, QUAD_LOBATTO, */
/*                   deg_lobatto, dvec_di_quad, quad_type, deg_quad, (P4EST_DIM)); */

/*     for (int j = 0; j < (P4EST_DIM); j++){ */
/*       for (int k = 0; k < volume_nodes_quad; k++){ */
/*        dvec[j][k] += rst_xyz_quad[i][j][k]*dvec_di_quad[k]; */
/*       } */
/*     } */
/*   } */

/*   P4EST_FREE(dvec_di_prolonged); */
/*   P4EST_FREE(dvec_di_quad); */
/*   P4EST_FREE(dvec_di_lobatto); */
/* } */
