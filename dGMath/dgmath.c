#include "dgmath.h"
#include "dgmath_nodes_and_weights.h"
#include <assert.h>
#include <math.h>
#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../pXest/pXest.h"

static const int dgmath_p8est_FToF_code[6][6] = {
  {0, 1, 1, 0, 0, 1}, {2, 0, 0, 1, 1, 0}, {2, 0, 0, 1, 1, 0},
  {0, 2, 2, 0, 0, 1}, {0, 2, 2, 0, 0, 1}, {2, 0, 0, 2, 2, 0}};

static const int dgmath_p8est_code_to_perm[3][4] = {
  {1, 2, 5, 6}, {0, 3, 4, 7}, {0, 4, 3, 7}};

static const int dgmath_p8est_perm_to_order[8][4] = {
  {0, 1, 2, 3}, {0, 2, 1, 3}, {1, 0, 3, 2}, {1, 3, 0, 2},
  {2, 0, 3, 1}, {2, 3, 0, 1}, {3, 1, 2, 0}, {3, 2, 1, 0}};
  

#define ORDERING 3

dgmath_jit_dbase_t* dgmath_jit_dbase_init() {
  dgmath_jit_dbase_t* dgbase = P4EST_ALLOC(dgmath_jit_dbase_t, 1);
  int deg_max = dgbase->dgmath_max_storage = (MAX_DEGREE);
  
  dgbase->dgmath_ref_Mij_1d_table = P4EST_ALLOC(double*, deg_max);
  dgbase->dgmath_ref_invMij_1d_table = P4EST_ALLOC(double*, deg_max);
  dgbase->dgmath_ref_invVij_1d_table = P4EST_ALLOC(double*, deg_max);
  dgbase->dgmath_ref_invGaussVij_1d_table = P4EST_ALLOC(double*, deg_max);

  dgbase->dgmath_GLL_nodes_1d_table = P4EST_ALLOC(double*, deg_max);
  dgbase->dgmath_GLL_weights_1d_table = P4EST_ALLOC(double*, deg_max);
  dgbase->dgmath_GL_nodes_1d_table = P4EST_ALLOC(double*, deg_max);
  dgbase->dgmath_GL_weights_1d_table = P4EST_ALLOC(double*, deg_max);
  
  dgbase->dgmath_ref_Dij_1d_table = P4EST_ALLOC(double*, deg_max);
  dgbase->dgmath_LIFT_1d_table = P4EST_ALLOC(double*, deg_max);
  dgbase->dgmath_ref_xyz_nd_table = P4EST_ALLOC(double*, deg_max*(2)); /* 2-> DIM = 2 or 3  */
  dgbase->dgmath_ref_Gauss_xyz_nd_table = P4EST_ALLOC(double*, deg_max*(2)); /* 2-> DIM = 2 or 3  */
  dgbase->dgmath_FLIP_1d_table = P4EST_ALLOC(double*, deg_max);
  dgbase->dgmath_hp_prolong_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  dgbase->dgmath_p_prolong_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  dgbase->dgmath_p_prolong_1d_inverse_table = P4EST_ALLOC(double*, deg_max * deg_max);
  dgbase->dgmath_hp_restrict_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  dgbase->dgmath_p_restrict_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  dgbase->dgmath_hp_prolong_transpose_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  dgbase->dgmath_p_prolong_transpose_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);  
  dgbase->dgmath_p_prolong_transpose_1d_inverse_table = P4EST_ALLOC(double*, deg_max * deg_max);  
  dgbase->dgmath_ref_GaussVij_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  dgbase->dgmath_ref_GaussDij_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);

  dgbase->dgmath_ref_GLL_to_GL_interp_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  dgbase->dgmath_ref_GLL_to_GL_interp_1d_inverse_table = P4EST_ALLOC(double*, deg_max * deg_max);
  dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_inverse_table = P4EST_ALLOC(double*, deg_max * deg_max);
  dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  /* dgbase->dgmath_ref_GL_to_GLL_interp_trans_1d_table = P4EST_ALLOC(double*, deg_max * deg_max); */
  
  dgbase->dgmath_hp_restrict_interp_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);
  dgbase->dgmath_p_restrict_interp_1d_table = P4EST_ALLOC(double*, deg_max * deg_max);

  int i;
  for (i = 0; i < dgbase->dgmath_max_storage * dgbase->dgmath_max_storage; i++){
    dgbase->dgmath_hp_prolong_1d_table[i] = NULL;
    dgbase->dgmath_p_prolong_1d_table[i] = NULL;
    dgbase->dgmath_p_prolong_1d_inverse_table[i] = NULL;
    dgbase->dgmath_hp_restrict_1d_table[i] = NULL;
    dgbase->dgmath_hp_prolong_transpose_1d_table[i] = NULL;
    dgbase->dgmath_p_restrict_1d_table[i] = NULL;
    dgbase->dgmath_p_prolong_transpose_1d_table[i] = NULL;
    dgbase->dgmath_p_prolong_transpose_1d_inverse_table[i] = NULL;
    dgbase->dgmath_hp_restrict_interp_1d_table[i] = NULL;
    dgbase->dgmath_p_restrict_interp_1d_table[i] = NULL;
    dgbase->dgmath_ref_GaussVij_1d_table[i] = NULL;
    dgbase->dgmath_ref_GaussDij_1d_table[i] = NULL;
    dgbase->dgmath_ref_GLL_to_GL_interp_1d_table[i] = NULL;
    dgbase->dgmath_ref_GLL_to_GL_interp_1d_inverse_table[i] = NULL;
    dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_table[i] = NULL;
    dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_inverse_table[i] = NULL;
    if (i < dgbase->dgmath_max_storage){
      dgbase->dgmath_ref_Mij_1d_table[i] = NULL;
      dgbase->dgmath_GLL_nodes_1d_table[i] = NULL;
      dgbase->dgmath_GLL_weights_1d_table[i] = NULL;
      dgbase->dgmath_GL_nodes_1d_table[i] = NULL;
      dgbase->dgmath_GL_weights_1d_table[i] = NULL;
      dgbase->dgmath_ref_invMij_1d_table[i] = NULL;
      dgbase->dgmath_ref_invVij_1d_table[i] = NULL;
      dgbase->dgmath_ref_invGaussVij_1d_table[i] = NULL;
      dgbase->dgmath_ref_Dij_1d_table[i] = NULL;
      dgbase->dgmath_LIFT_1d_table[i] = NULL;
      dgbase->dgmath_FLIP_1d_table[i] = NULL;
    }
    if (i < dgbase->dgmath_max_storage * (2)){
      dgbase->dgmath_ref_xyz_nd_table[i] = NULL;    
      dgbase->dgmath_ref_Gauss_xyz_nd_table[i] = NULL;
    }
  }

  return dgbase;
}

void dgmath_jit_dbase_destroy(dgmath_jit_dbase_t* dgbase) {

  int i;
  for (i = 0; i < dgbase->dgmath_max_storage; i++) {
    P4EST_FREE(dgbase->dgmath_ref_Mij_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_GLL_nodes_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_GLL_weights_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_GL_nodes_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_GL_weights_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_ref_invMij_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_ref_invVij_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_ref_invGaussVij_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_ref_Dij_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_LIFT_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_FLIP_1d_table[i]);
  }
  P4EST_FREE(dgbase->dgmath_ref_Mij_1d_table);
  P4EST_FREE(dgbase->dgmath_GLL_nodes_1d_table);
  P4EST_FREE(dgbase->dgmath_GLL_weights_1d_table);
  P4EST_FREE(dgbase->dgmath_GL_nodes_1d_table);
  P4EST_FREE(dgbase->dgmath_GL_weights_1d_table);
  P4EST_FREE(dgbase->dgmath_ref_invMij_1d_table);
  P4EST_FREE(dgbase->dgmath_ref_invVij_1d_table);
  P4EST_FREE(dgbase->dgmath_ref_invGaussVij_1d_table);
  P4EST_FREE(dgbase->dgmath_ref_Dij_1d_table);
  P4EST_FREE(dgbase->dgmath_LIFT_1d_table);
  P4EST_FREE(dgbase->dgmath_FLIP_1d_table);

  
  for (i = 0; i < dgbase->dgmath_max_storage * dgbase->dgmath_max_storage; i++) {
    P4EST_FREE(dgbase->dgmath_hp_prolong_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_hp_prolong_transpose_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_p_prolong_transpose_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_p_prolong_transpose_1d_inverse_table[i]);
    P4EST_FREE(dgbase->dgmath_p_prolong_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_p_prolong_1d_inverse_table[i]);
    P4EST_FREE(dgbase->dgmath_hp_restrict_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_hp_restrict_interp_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_p_restrict_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_p_restrict_interp_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_ref_GaussVij_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_ref_GaussDij_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_ref_GLL_to_GL_interp_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_ref_GLL_to_GL_interp_1d_inverse_table[i]);
    P4EST_FREE(dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_table[i]);
    P4EST_FREE(dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_inverse_table[i]);
  }
  P4EST_FREE(dgbase->dgmath_hp_prolong_1d_table);
  P4EST_FREE(dgbase->dgmath_hp_prolong_transpose_1d_table);
  P4EST_FREE(dgbase->dgmath_p_prolong_transpose_1d_table);
  P4EST_FREE(dgbase->dgmath_p_prolong_transpose_1d_inverse_table);
  P4EST_FREE(dgbase->dgmath_hp_restrict_1d_table);
  P4EST_FREE(dgbase->dgmath_hp_restrict_interp_1d_table);
  P4EST_FREE(dgbase->dgmath_p_restrict_1d_table);
  P4EST_FREE(dgbase->dgmath_p_restrict_interp_1d_table);
  P4EST_FREE(dgbase->dgmath_ref_GaussVij_1d_table);
  P4EST_FREE(dgbase->dgmath_ref_GaussDij_1d_table);
  P4EST_FREE(dgbase->dgmath_ref_GLL_to_GL_interp_1d_table);
  P4EST_FREE(dgbase->dgmath_ref_GLL_to_GL_interp_1d_inverse_table);
  P4EST_FREE(dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_table);
  P4EST_FREE(dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_inverse_table);
  P4EST_FREE(dgbase->dgmath_p_prolong_1d_table);
  P4EST_FREE(dgbase->dgmath_p_prolong_1d_inverse_table);
  
  for (i = 0; i < dgbase->dgmath_max_storage * (2); i++) {
    P4EST_FREE(dgbase->dgmath_ref_xyz_nd_table[i]);
    P4EST_FREE(dgbase->dgmath_ref_Gauss_xyz_nd_table[i]);
  }
  P4EST_FREE(dgbase->dgmath_ref_Gauss_xyz_nd_table);
  P4EST_FREE(dgbase->dgmath_ref_xyz_nd_table);
  P4EST_FREE(dgbase);
}

int dgmath_get_nodes(int dim, int deg) {
  mpi_assert(dim > 0 && dim < 4 && deg > 0);
  int i;
  int volume_nodes = 1;
  for (i = 0; i < dim; i++) volume_nodes *= deg + 1;
  return volume_nodes;
}

int dgmath_is_child_left_or_right(int c, int dir) {
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
    mpi_abort("Error: dgmath_is_child_left_or_right\n");
    return -1;
  }
  /* printf("Error: dgmath_is_child_left_or_right\n"); */
}

void dgmath_child_to_parent_ref_coords(int c, double* r) {
  *r *= .5;
  if (c == 0) {
    *r -= .5;
  } else if (c == 1) {
    *r += .5;
  } else {
    mpi_abort("Error: dgmath_child_to_parent_ref_coords\n");
    /* printf("ERROR: Boundary translation\n"); */
  }
}


void dgmath_build_Vij_1d(dgmath_jit_dbase_t* dgbase, double* Vij_1d,
                                int deg) {
  int i, j, rows, cols;
  const double* lobatto_nodes = dgmath_fetch_GLL_nodes_1d(dgbase, deg);
  rows = cols = deg + 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      Vij_1d[i * cols + j] = dgmath_jacobi(lobatto_nodes[i], 0., 0., j);
}



static void dgmath_build_invVij_1d(dgmath_jit_dbase_t* dgbase,
                                   double* invVij_1d, int deg)
{
  dgmath_build_Vij_1d(dgbase, invVij_1d, deg);
  linalg_invert(invVij_1d, deg + 1);
}


double* dgmath_fetch_invVij_1d(dgmath_jit_dbase_t* dgbase, int deg) {
  mpi_assert(deg < dgbase->dgmath_max_storage);
  if (dgbase->dgmath_ref_invVij_1d_table[deg] != NULL) {
    return dgbase->dgmath_ref_invVij_1d_table[deg];
  }
  else {
    int size = (deg + 1) * (deg + 1);
    dgbase->dgmath_ref_invVij_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_ref_invVij_1d_table[deg];
    dgmath_build_invVij_1d(dgbase, op, deg);
    return op;
  }
}

void dgmath_hp_apply_nd_prolong_with_ptr(double* Uh, int degh, double* UH,
                                         int degH, int dim, int c,
                                         double* hp_prolong_matrix_1d) {
  mpi_assert((degH <= degh) && (c < (1 << dim)));
  int nodesh = degh + 1;
  int nodesH = degH + 1;

  if (dim == 1) {
    linalg_matvec_plus_vec(1.0, &hp_prolong_matrix_1d[c * nodesh * nodesH], UH,
                           0., Uh, nodesh, nodesH);
  } else if (dim == 2) {
    int cx = dgmath_is_child_left_or_right(c, 0);
    int cy = dgmath_is_child_left_or_right(c, 1);
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
    int cx = dgmath_is_child_left_or_right(c, 0);
    int cy = dgmath_is_child_left_or_right(c, 1);
    int cz = dgmath_is_child_left_or_right(c, 2);
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
    mpi_abort("ERROR: dgmath_prolong_nd_U\n");
  }
}


void dgmath_build_ref_GaussVij_1d
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* ref_Vij_Gauss_1d,
 int Lobatto_degree,
 int Gauss_degree
)
{
  double* Gauss_nodes = dgmath_fetch_GL_nodes_1d(dgmath_jit_dbase, Gauss_degree);
  
  int rows = Gauss_degree + 1;
  int cols = Lobatto_degree + 1;
  
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      ref_Vij_Gauss_1d[i * cols + j] = dgmath_jacobi(Gauss_nodes[i], 0., 0., j);

  // util_print_matrix(Gauss_nodes, Gauss_degree + 1, 1, "Gauss_nodes = ", 0);
  
}

double* dgmath_fetch_ref_GaussVij_1d(dgmath_jit_dbase_t* dgbase, int deg_Lobatto,
                                         int deg_Gauss) {
  mpi_assert(deg_Gauss <= dgbase->dgmath_max_storage);
  mpi_assert(deg_Lobatto <= dgbase->dgmath_max_storage);

  if (dgbase->dgmath_ref_GaussVij_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] != NULL) {
    return dgbase->dgmath_ref_GaussVij_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss];
  }
  else {
    int size = (deg_Gauss + 1) * (deg_Lobatto + 1);
    dgbase->dgmath_ref_GaussVij_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_ref_GaussVij_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss];
    dgmath_build_ref_GaussVij_1d(dgbase, op, deg_Lobatto, deg_Gauss);
    return op;
  }
}

void dgmath_build_ref_invGaussVij_1d
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* ref_invGaussVij_1d,
 int Gauss_degree
)
{
  double* GaussVij = dgmath_fetch_ref_GaussVij_1d(dgmath_jit_dbase, Gauss_degree, Gauss_degree);
  linalg_invert_and_copy(GaussVij, ref_invGaussVij_1d, Gauss_degree + 1);
}


double* dgmath_fetch_invGaussVij_1d(dgmath_jit_dbase_t* dgbase, int deg) {
  mpi_assert(deg < dgbase->dgmath_max_storage);
  if (dgbase->dgmath_ref_invGaussVij_1d_table[deg] != NULL) {
    return dgbase->dgmath_ref_invGaussVij_1d_table[deg];
  }
  else {
    int size = (deg + 1) * (deg + 1);
    dgbase->dgmath_ref_invGaussVij_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_ref_invGaussVij_1d_table[deg];
    dgmath_build_ref_invGaussVij_1d(dgbase, op, deg);
    return op;
  }
}

void dgmath_build_ref_GLL_to_GL_interp_1d
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* ref_GLL_to_GL_interp_1d,
 int Lobatto_degree,
 int Gauss_degree
)
{
  double* Gauss_nodes = dgmath_fetch_GL_nodes_1d(dgmath_jit_dbase, Gauss_degree);
  double* ref_GaussVij = P4EST_ALLOC(double, (Gauss_degree + 1)*(Lobatto_degree+1));
  
  int rows = Gauss_degree + 1;
  int cols = Lobatto_degree + 1;
  
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      ref_GaussVij[i * cols + j] = dgmath_jacobi(Gauss_nodes[i], 0., 0., j);

  double* invVij = dgmath_fetch_invVij_1d(dgmath_jit_dbase, Lobatto_degree);
  linalg_mat_multiply(ref_GaussVij, invVij, ref_GLL_to_GL_interp_1d, Gauss_degree + 1, Lobatto_degree + 1, Lobatto_degree + 1);

  P4EST_FREE(ref_GaussVij);
}



/* void dgmath_build_ref_GL_to_GLL_interp_1d */
/* ( */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase, */
/*  double* ref_GL_to_GLL_interp_1d, */
/*  int Lobatto_degree, */
/*  int Gauss_degree */
/* ) */
/* { */
/*   double* Lobatto_nodes = dgmath_fetch_GLL_nodes_1d(dgmath_jit_dbase, Gauss_degree); */
/*   double* ref_LobattoVij = P4EST_ALLOC(double, (Gauss_degree + 1)*(Lobatto_degree+1)); */
  
/*   int cols = Gauss_degree + 1; */
/*   int rows = Lobatto_degree + 1; */
  
/*   for (int i = 0; i < rows; i++) */
/*     for (int j = 0; j < cols; j++) */
/*       ref_LobattoVij[i * cols + j] = dgmath_jacobi(Lobatto_nodes[i], 0., 0., j); */

/*   double* invGaussVij = dgmath_fetch_invGaussVij_1d(dgmath_jit_dbase, Gauss_degree); */
/*   linalg_mat_multiply(ref_LobattoVij, invGaussVij, ref_GL_to_GLL_interp_1d, Lobatto_degree + 1, Gauss_degree + 1, Gauss_degree + 1); */

/*   P4EST_FREE(ref_LobattoVij); */
/* } */



double*
dgmath_fetch_ref_GLL_to_GL_interp_1d(dgmath_jit_dbase_t* dgbase, int deg_Lobatto,
                                         int deg_Gauss) {
  mpi_assert(deg_Gauss <= dgbase->dgmath_max_storage);
  mpi_assert(deg_Lobatto <= dgbase->dgmath_max_storage);

  if (dgbase->dgmath_ref_GLL_to_GL_interp_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] != NULL) {
    return dgbase->dgmath_ref_GLL_to_GL_interp_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss];
  }
  else {
    int size = (deg_Gauss + 1) * (deg_Lobatto + 1);
    dgbase->dgmath_ref_GLL_to_GL_interp_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_ref_GLL_to_GL_interp_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss];
    dgmath_build_ref_GLL_to_GL_interp_1d(dgbase, op, deg_Lobatto, deg_Gauss);
    return op;
  }
}

double* dgmath_fetch_ref_GLL_to_GL_interp_1d_inverse(dgmath_jit_dbase_t* dgbase, int deg_Lobatto,
                                         int deg_Gauss) {
  mpi_assert(deg_Gauss <= dgbase->dgmath_max_storage);
  mpi_assert(deg_Lobatto <= dgbase->dgmath_max_storage);

  if (dgbase->dgmath_ref_GLL_to_GL_interp_1d_inverse_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] != NULL) {
    return dgbase->dgmath_ref_GLL_to_GL_interp_1d_inverse_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss];
  }
  else {
    int size = (deg_Gauss + 1) * (deg_Lobatto + 1);
    dgbase->dgmath_ref_GLL_to_GL_interp_1d_inverse_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_ref_GLL_to_GL_interp_1d_inverse_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss];
    dgmath_build_ref_GLL_to_GL_interp_1d_inverse(dgbase, op, deg_Lobatto, deg_Gauss);
    return op;
  }
}


void dgmath_build_ref_GLL_to_GL_interp_1d_inverse
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* ref_GLL_to_GL_interp_1d_inverse,
 int Lobatto_degree,
 int Gauss_degree
)
{
  double* ref_GLL_to_GL_interp_1d = dgmath_fetch_ref_GLL_to_GL_interp_1d(dgmath_jit_dbase, Lobatto_degree, Gauss_degree);
  linalg_leftinverse(ref_GLL_to_GL_interp_1d, ref_GLL_to_GL_interp_1d_inverse, Gauss_degree + 1, Lobatto_degree + 1);  
}

/* double* dgmath_fetch_ref_GL_to_GLL_interp_1d(dgmath_jit_dbase_t* dgbase, int deg_Lobatto, */
/*                                          int deg_Gauss) { */
/*   mpi_assert(deg_Gauss <= dgbase->dgmath_max_storage); */
/*   mpi_assert(deg_Lobatto <= dgbase->dgmath_max_storage); */

/*   if (dgbase->dgmath_ref_GL_to_GLL_interp_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] != NULL) { */
/*     return dgbase->dgmath_ref_GL_to_GLL_interp_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss]; */
/*   } */
/*   else { */
/*     int size = (deg_Gauss + 1) * (deg_Lobatto + 1); */
/*     dgbase->dgmath_ref_GL_to_GLL_interp_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size); */
/*     double* op = dgbase->dgmath_ref_GL_to_GLL_interp_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss]; */
/*     dgmath_build_ref_GL_to_GLL_interp_1d(dgbase, op, deg_Lobatto, deg_Gauss); */
/*     return op; */
/*   } */
/* } */




/* void dgmath_build_ref_GL_to_GLL_interp_trans_1d */
/* ( */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase, */
/*  double* ref_GL_to_GLL_interp_trans_1d, */
/*  int Lobatto_degree, */
/*  int Gauss_degree */
/* ) */
/* { */
/*   double* ref_GL_to_GLL_interp_1d = dgmath_fetch_ref_GL_to_GLL_interp_1d(dgmath_jit_dbase, Lobatto_degree, Gauss_degree); */
/*   linalg_mat_transpose_nonsqr(ref_GL_to_GLL_interp_1d, ref_GL_to_GLL_interp_trans_1d, Lobatto_degree + 1, Gauss_degree + 1); */
/* } */



double* dgmath_fetch_ref_GLL_to_GL_interp_trans_1d(dgmath_jit_dbase_t* dgbase, int deg_Lobatto,
                                         int deg_Gauss) {
  mpi_assert(deg_Gauss <= dgbase->dgmath_max_storage);
  mpi_assert(deg_Lobatto <= dgbase->dgmath_max_storage);

  if (dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] != NULL) {
    return dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss];
  }
  else {
    int size = (deg_Gauss + 1) * (deg_Lobatto + 1);
    dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss];
    dgmath_build_ref_GLL_to_GL_interp_trans_1d(dgbase, op, deg_Lobatto, deg_Gauss);
    return op;
  }
}


double* dgmath_fetch_ref_GLL_to_GL_interp_trans_1d_inverse(dgmath_jit_dbase_t* dgbase, int deg_Lobatto,
                                         int deg_Gauss) {
  mpi_assert(deg_Gauss <= dgbase->dgmath_max_storage);
  mpi_assert(deg_Lobatto <= dgbase->dgmath_max_storage);

  if (dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_inverse_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] != NULL) {
    return dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_inverse_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss];
  }
  else {
    int size = (deg_Gauss + 1) * (deg_Lobatto + 1);
    dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_inverse_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_ref_GLL_to_GL_interp_trans_1d_inverse_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss];
    dgmath_build_ref_GLL_to_GL_interp_trans_1d_inverse(dgbase, op, deg_Lobatto, deg_Gauss);
    return op;
  }
}


void dgmath_build_ref_GLL_to_GL_interp_trans_1d
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* ref_GLL_to_GL_interp_trans_1d,
 int Lobatto_degree,
 int Gauss_degree
)
{
  double* ref_GLL_to_GL_interp_1d = dgmath_fetch_ref_GLL_to_GL_interp_1d(dgmath_jit_dbase, Lobatto_degree, Gauss_degree);
  linalg_mat_transpose_nonsqr(ref_GLL_to_GL_interp_1d, ref_GLL_to_GL_interp_trans_1d, Gauss_degree + 1, Lobatto_degree + 1);
}


void dgmath_build_ref_GLL_to_GL_interp_trans_1d_inverse
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* ref_GLL_to_GL_interp_trans_1d_inverse,
 int Lobatto_degree,
 int Gauss_degree
)
{
  double* ref_GLL_to_GL_interp_trans_1d = dgmath_fetch_ref_GLL_to_GL_interp_trans_1d(dgmath_jit_dbase, Lobatto_degree, Gauss_degree);
  linalg_leftinverse(ref_GLL_to_GL_interp_trans_1d, ref_GLL_to_GL_interp_trans_1d_inverse, Lobatto_degree + 1, Gauss_degree + 1);
}


/* double* dgmath_fetch_ref_GL_to_GLL_interp_trans_1d(dgmath_jit_dbase_t* dgbase, int deg_Lobatto, */
/*                                          int deg_Gauss) { */
/*   mpi_assert(deg_Gauss <= dgbase->dgmath_max_storage); */
/*   mpi_assert(deg_Lobatto <= dgbase->dgmath_max_storage); */

/*   if (dgbase->dgmath_ref_GL_to_GLL_interp_trans_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] != NULL) { */
/*     return dgbase->dgmath_ref_GL_to_GLL_interp_trans_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss]; */
/*   } */
/*   else { */
/*     int size = (deg_Gauss + 1) * (deg_Lobatto + 1); */
/*     dgbase->dgmath_ref_GL_to_GLL_interp_trans_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size); */
/*     double* op = dgbase->dgmath_ref_GL_to_GLL_interp_trans_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss]; */
/*     dgmath_build_ref_GL_to_GLL_interp_trans_1d(dgbase, op, deg_Lobatto, deg_Gauss); */
/*     return op; */
/*   } */
/* } */




void dgmath_interp_GLL_to_GL
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* u_Lobatto_in,
 int deg_Lobatto,
 int deg_Gauss,
 double* u_Gauss_out,
 int dim
)
{
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  double* GLL_to_GL_interp = dgmath_fetch_ref_GLL_to_GL_interp_1d(dgmath_jit_dbase, deg_Lobatto, deg_Gauss);

  int nodes_Lobatto = deg_Lobatto + 1;
  int nodes_Gauss = deg_Gauss + 1;
  
  if (dim == 1){
    linalg_matvec_plus_vec(1.0, GLL_to_GL_interp, u_Lobatto_in, 0., u_Gauss_out, nodes_Gauss, nodes_Lobatto);
  }
  else if (dim == 2) {
    linalg_kron_A1A2x_NONSQR(u_Gauss_out, GLL_to_GL_interp, GLL_to_GL_interp, u_Lobatto_in, nodes_Gauss, nodes_Lobatto,
                             nodes_Gauss, nodes_Lobatto);
  }
  else if (dim == 3){
    linalg_kron_A1A2A3x_NONSQR(u_Gauss_out, GLL_to_GL_interp, GLL_to_GL_interp, GLL_to_GL_interp, u_Lobatto_in,
                               nodes_Gauss, nodes_Lobatto, nodes_Gauss, nodes_Lobatto, nodes_Gauss, nodes_Lobatto);
  }
  else {
    mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension.");
  }
}



/* void dgmath_interp_GL_to_GLL */
/* ( */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase, */
/*  double* u_Gauss_in, */
/*  int deg_Gauss, */
/*  int deg_Lobatto, */
/*  double* u_Lobatto_out, */
/*  int dim */
/* ) */
/* { */
/*   mpi_assert(dim == 1 || dim == 2 || dim == 3); */
/*   double* GL_to_GLL_interp = dgmath_fetch_ref_GL_to_GLL_interp_1d(dgmath_jit_dbase, deg_Lobatto, deg_Gauss); */

/*   int nodes_Gauss = deg_Gauss + 1; */
/*   int nodes_Lobatto = deg_Lobatto + 1; */
  
/*   if (dim == 1){ */
/*     linalg_matvec_plus_vec(1.0, GL_to_GLL_interp, u_Gauss_in, 0., u_Lobatto_out, nodes_Lobatto, nodes_Gauss); */
/*   } */
/*   else if (dim == 2) { */
/*     linalg_kron_A1A2x_NONSQR(u_Lobatto_out, GL_to_GLL_interp, GL_to_GLL_interp, u_Gauss_in, nodes_Lobatto, nodes_Gauss, */
/*                              nodes_Lobatto, nodes_Gauss); */
/*   } */
/*   else if (dim == 3){ */
/*     linalg_kron_A1A2A3x_NONSQR(u_Lobatto_out, GL_to_GLL_interp, GL_to_GLL_interp, GL_to_GLL_interp, u_Gauss_in, */
/*                                nodes_Lobatto, nodes_Gauss, nodes_Lobatto, nodes_Gauss, nodes_Lobatto, nodes_Gauss); */
/*   } */
/*   else { */
/*     mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension."); */
/*   } */
/* } */


void dgmath_apply_curvedGaussMass
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* in,
 int deg_Lobatto,
 double* jac_Gauss,
 int deg_Gauss,
 int dim,
 double* out
){
  mpi_assert(dim == 1 || dim == 2 || dim == 3);

  /* for now assert this, can get rid of by p-prolonging then p-restricting */
  /* mpi_assert(deg_Lobatto == deg_Gauss); */
  int volume_nodes_Gauss = dgmath_get_nodes(dim, deg_Gauss);

  double* in_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
  double* w_j_in_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
  
  double* Gauss_weights = dgmath_fetch_GL_weights_1d(dgmath_jit_dbase, deg_Gauss);
  double* GLL_to_GL_interp_trans = dgmath_fetch_ref_GLL_to_GL_interp_trans_1d
                                   (dgmath_jit_dbase, deg_Lobatto, deg_Gauss);
  
  dgmath_interp_GLL_to_GL(dgmath_jit_dbase, in, deg_Lobatto, deg_Gauss, in_Gauss, dim);

  int nodes_Lobatto = deg_Lobatto + 1;
  int nodes_Gauss = deg_Gauss + 1;
  
  if (dim == 1){
    linalg_kron_vec_dot_xy(Gauss_weights, jac_Gauss, in_Gauss, nodes_Gauss, w_j_in_Gauss);
    linalg_matvec_plus_vec(1.0, GLL_to_GL_interp_trans, w_j_in_Gauss, 0., out, nodes_Lobatto, nodes_Gauss);
  }
  else if (dim == 2) {
    linalg_kron_vec_o_vec_dot_xy(Gauss_weights, jac_Gauss, in_Gauss, nodes_Gauss, w_j_in_Gauss);
    linalg_kron_A1A2x_NONSQR(out, GLL_to_GL_interp_trans, GLL_to_GL_interp_trans, w_j_in_Gauss, nodes_Lobatto, nodes_Gauss,
                             nodes_Lobatto, nodes_Gauss);
  } else if (dim == 3) {
    linalg_kron_vec_o_vec_o_vec_dot_xy(Gauss_weights, jac_Gauss, in_Gauss, nodes_Gauss, w_j_in_Gauss);
    linalg_kron_A1A2A3x_NONSQR(out, GLL_to_GL_interp_trans, GLL_to_GL_interp_trans, GLL_to_GL_interp_trans, w_j_in_Gauss,
                               nodes_Lobatto, nodes_Gauss, nodes_Lobatto, nodes_Gauss, nodes_Lobatto, nodes_Gauss);
  } else {    
    P4EST_FREE(w_j_in_Gauss);
    P4EST_FREE(in_Gauss);
    mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension.");
  }  
  P4EST_FREE(w_j_in_Gauss);
  P4EST_FREE(in_Gauss);
}

void dgmath_apply_curvedLobattoMass
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* in,
 int deg_Lobatto,
 double* jac_Lobatto_integ,
 int deg_Lobatto_integ,
 int dim,
 double* out
){
  mpi_assert(dim == 1 || dim == 2 || dim == 3);

  /* for now assert this, can get rid of by p-prolonging then p-restricting */
  /* mpi_assert(deg_Lobatto == deg_Lobatto_integ); */
  int volume_nodes_Lobatto_integ = dgmath_get_nodes(dim, deg_Lobatto_integ);

  double* in_Lobatto_integ = P4EST_ALLOC(double, volume_nodes_Lobatto_integ);
  double* w_j_in_Lobatto_integ = P4EST_ALLOC(double, volume_nodes_Lobatto_integ);
  
  double* Lobatto_integ_weights = dgmath_fetch_GLL_weights_1d(dgmath_jit_dbase, deg_Lobatto_integ);
  double* p_prolong_trans = dgmath_fetch_p_prolong_transpose_1d
                                   (dgmath_jit_dbase, deg_Lobatto, deg_Lobatto_integ);
  
  dgmath_apply_p_prolong(dgmath_jit_dbase, in, deg_Lobatto, dim, deg_Lobatto_integ, in_Lobatto_integ);

  int nodes_Lobatto = deg_Lobatto + 1;
  int nodes_Lobatto_integ = deg_Lobatto_integ + 1;
  
  if (dim == 1){
    linalg_kron_vec_dot_xy(Lobatto_integ_weights, jac_Lobatto_integ, in_Lobatto_integ, nodes_Lobatto_integ, w_j_in_Lobatto_integ);
    linalg_matvec_plus_vec(1.0, p_prolong_trans, w_j_in_Lobatto_integ, 0., out, nodes_Lobatto, nodes_Lobatto_integ);
  }
  else if (dim == 2) {
    linalg_kron_vec_o_vec_dot_xy(Lobatto_integ_weights, jac_Lobatto_integ, in_Lobatto_integ, nodes_Lobatto_integ, w_j_in_Lobatto_integ);
    linalg_kron_A1A2x_NONSQR(out, p_prolong_trans, p_prolong_trans, w_j_in_Lobatto_integ, nodes_Lobatto, nodes_Lobatto_integ,
                             nodes_Lobatto, nodes_Lobatto_integ);
  } else if (dim == 3) {
    linalg_kron_vec_o_vec_o_vec_dot_xy(Lobatto_integ_weights, jac_Lobatto_integ, in_Lobatto_integ, nodes_Lobatto_integ, w_j_in_Lobatto_integ);
    linalg_kron_A1A2A3x_NONSQR(out, p_prolong_trans, p_prolong_trans, p_prolong_trans, w_j_in_Lobatto_integ,
                               nodes_Lobatto, nodes_Lobatto_integ, nodes_Lobatto, nodes_Lobatto_integ, nodes_Lobatto, nodes_Lobatto_integ);
  } else {    
    P4EST_FREE(w_j_in_Lobatto_integ);
    P4EST_FREE(in_Lobatto_integ);
    mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension.");
  }  
  P4EST_FREE(w_j_in_Lobatto_integ);
  P4EST_FREE(in_Lobatto_integ);
}


/* void dgmath_apply_curvedInverseLobattoMass */
/* ( */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase, */
/*  double* in, */
/*  int deg_Lobatto, */
/*  double* jac_Lobatto, */
/*  int dim, */
/*  double* out */
/* ) */
/* { */
/*   mpi_assert(dim == 1 || dim == 2 || dim == 3); */
/*   double* Lobatto_weights = dgmath_fetch_GLL_weights_1d(dgmath_jit_dbase, deg_Lobatto); */
/*   int nodes_Lobatto = deg_Lobatto + 1; */
  
/*   if (dim == 1){ */
/*     linalg_kron_oneover_vec_dot_oneover_x_dot_y(Lobatto_weights, jac_Lobatto, in, nodes_Lobatto, out); */
/*   } */
/*   else if (dim == 2) { */
/*     linalg_kron_oneover_vec_o_vec_dot_oneover_x_dot_y */
/*       ( */
/*        Lobatto_weights, */
/*        jac_Lobatto, */
/*        in, */
/*        nodes_Lobatto, */
/*        out */
/*       ); */
/*   } */
/*   else if (dim == 3) { */
/*     linalg_kron_oneover_vec_o_vec_o_vec_dot_oneover_x_dot_y */
/*       ( */
/*        Lobatto_weights, */
/*        jac_Lobatto, */
/*        in, */
/*        nodes_Lobatto, */
/*        out */
/*       ); */
/*   } */
/*   else { */
/*     mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension."); */
/*   } */
/* } */

/* void dgmath_apply_curvedLobattoMass */
/* ( */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase, */
/*  double* in, */
/*  int deg_Lobatto, */
/*  double* jac_Lobatto, */
/*  int dim, */
/*  double* out */
/* ) */
/* { */
/*   mpi_assert(dim == 1 || dim == 2 || dim == 3); */
/*   double* Lobatto_weights = dgmath_fetch_GLL_weights_1d(dgmath_jit_dbase, deg_Lobatto); */
/*   int nodes_Lobatto = deg_Lobatto + 1; */
  
/*   if (dim == 1){ */
/*     linalg_kron_vec_dot_xy(Lobatto_weights, jac_Lobatto, in, nodes_Lobatto, out); */
/*   } */
/*   else if (dim == 2) { */
/*     linalg_kron_vec_o_vec_dot_xy(Lobatto_weights, jac_Lobatto, in, nodes_Lobatto, out); */
/*   } */
/*   else if (dim == 3) { */
/*     linalg_kron_vec_o_vec_o_vec_dot_xy(Lobatto_weights, jac_Lobatto, in, nodes_Lobatto, out); */
/*   } */
/*   else {     */
/*     mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension."); */
/*   }   */
/* } */



void dgmath_apply_curvedGaussMass_onGaussNodeVec
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* in_Gauss,
 int deg_Lobatto,
 double* jac_Gauss,
 int deg_Gauss,
 int dim,
 double* out
){
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int volume_nodes_Gauss = dgmath_get_nodes(dim, deg_Gauss);

  /* double* in_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss); */
  double* w_j_in_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
  
  double* Gauss_weights = dgmath_fetch_GL_weights_1d(dgmath_jit_dbase, deg_Gauss);
  double* GLL_to_GL_interp_trans = dgmath_fetch_ref_GLL_to_GL_interp_trans_1d(dgmath_jit_dbase, deg_Lobatto, deg_Gauss);
  
  /* dgmath_interp_GLL_to_GL(dgmath_jit_dbase, in, deg_Lobatto, deg_Gauss, in_Gauss, dim); */

  int nodes_Lobatto = deg_Lobatto + 1;
  int nodes_Gauss = deg_Gauss + 1;
  
  if (dim == 1){
    linalg_kron_vec_dot_xy(Gauss_weights, jac_Gauss, in_Gauss, nodes_Gauss, w_j_in_Gauss);
    linalg_matvec_plus_vec(1.0, GLL_to_GL_interp_trans, w_j_in_Gauss, 0., out, nodes_Lobatto, nodes_Gauss);
  }
  else if (dim == 2) {
    linalg_kron_vec_o_vec_dot_xy(Gauss_weights, jac_Gauss, in_Gauss, nodes_Gauss, w_j_in_Gauss);
    linalg_kron_A1A2x_NONSQR(out, GLL_to_GL_interp_trans, GLL_to_GL_interp_trans, w_j_in_Gauss, nodes_Lobatto, nodes_Gauss,
                             nodes_Lobatto, nodes_Gauss);
  }
  else if (dim == 3) {
    linalg_kron_vec_o_vec_o_vec_dot_xy(Gauss_weights, jac_Gauss, in_Gauss, nodes_Gauss, w_j_in_Gauss);
    linalg_kron_A1A2A3x_NONSQR(out, GLL_to_GL_interp_trans, GLL_to_GL_interp_trans, GLL_to_GL_interp_trans, w_j_in_Gauss,
                               nodes_Lobatto, nodes_Gauss, nodes_Lobatto, nodes_Gauss, nodes_Lobatto, nodes_Gauss);
  }
  else {    
    P4EST_FREE(w_j_in_Gauss);
    P4EST_FREE(in_Gauss);
    mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension.");
  }  
  P4EST_FREE(w_j_in_Gauss);
}

void dgmath_apply_curvedLobattoMass_onLobattoIntegNodeVec
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* in_Lobatto_integ,
 int deg_Lobatto,
 double* jac_Lobatto_integ,
 int deg_Lobatto_integ,
 int dim,
 double* out
){
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int volume_nodes_Lobatto_integ = dgmath_get_nodes(dim, deg_Lobatto_integ);

  /* double* in_Lobatto_integ = P4EST_ALLOC(double, volume_nodes_Lobatto_integ); */
  double* w_j_in_Lobatto_integ = P4EST_ALLOC(double, volume_nodes_Lobatto_integ);
  
  double* Lobatto_integ_weights = dgmath_fetch_GLL_weights_1d(dgmath_jit_dbase, deg_Lobatto_integ);
  double* p_prolong_trans = dgmath_fetch_p_prolong_transpose_1d(dgmath_jit_dbase, deg_Lobatto, deg_Lobatto_integ);
  
  /* dgmath_interp_GLL_to_GL(dgmath_jit_dbase, in, deg_Lobatto, deg_Lobatto_integ, in_Lobatto_integ, dim); */

  int nodes_Lobatto = deg_Lobatto + 1;
  int nodes_Lobatto_integ = deg_Lobatto_integ + 1;
  
  if (dim == 1){
    linalg_kron_vec_dot_xy(Lobatto_integ_weights, jac_Lobatto_integ, in_Lobatto_integ, nodes_Lobatto_integ, w_j_in_Lobatto_integ);
    linalg_matvec_plus_vec(1.0, p_prolong_trans, w_j_in_Lobatto_integ, 0., out, nodes_Lobatto, nodes_Lobatto_integ);
  }
  else if (dim == 2) {
    linalg_kron_vec_o_vec_dot_xy(Lobatto_integ_weights, jac_Lobatto_integ, in_Lobatto_integ, nodes_Lobatto_integ, w_j_in_Lobatto_integ);
    linalg_kron_A1A2x_NONSQR(out, p_prolong_trans, p_prolong_trans, w_j_in_Lobatto_integ, nodes_Lobatto, nodes_Lobatto_integ,
                             nodes_Lobatto, nodes_Lobatto_integ);
  }
  else if (dim == 3) {
    linalg_kron_vec_o_vec_o_vec_dot_xy(Lobatto_integ_weights, jac_Lobatto_integ, in_Lobatto_integ, nodes_Lobatto_integ, w_j_in_Lobatto_integ);
    linalg_kron_A1A2A3x_NONSQR(out, p_prolong_trans, p_prolong_trans, p_prolong_trans, w_j_in_Lobatto_integ,
                               nodes_Lobatto, nodes_Lobatto_integ, nodes_Lobatto, nodes_Lobatto_integ, nodes_Lobatto, nodes_Lobatto_integ);
  }
  else {    
    P4EST_FREE(w_j_in_Lobatto_integ);
    P4EST_FREE(in_Lobatto_integ);
    mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension.");
  }  
  P4EST_FREE(w_j_in_Lobatto_integ);
}

void dgmath_compute_curvedInverseGaussMass
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 int deg_Lobatto,
 int deg_Gauss,
 int dim,
 double* jac_Gauss,
 double* invM
)
{
  int volume_nodes_Lobatto = dgmath_get_nodes(dim,deg_Lobatto);
  
  double* M = P4EST_ALLOC(double, volume_nodes_Lobatto*volume_nodes_Lobatto);
  double* u = P4EST_ALLOC_ZERO(double, volume_nodes_Lobatto);
  double* Mu = P4EST_ALLOC_ZERO(double, volume_nodes_Lobatto);

  for (int i = 0; i < volume_nodes_Lobatto; i++){
    u[i] = 1.;
    dgmath_apply_curvedGaussMass(dgmath_jit_dbase, u, deg_Lobatto, jac_Gauss, deg_Gauss, dim, Mu);
    linalg_set_column(M, Mu, i, volume_nodes_Lobatto, volume_nodes_Lobatto);
    u[i] = 0.;
  }
  linalg_invert_and_copy(M, invM, volume_nodes_Lobatto);

  P4EST_FREE(Mu);
  P4EST_FREE(u);
  P4EST_FREE(M);
}

void dgmath_compute_curvedGaussMass
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 int deg_Lobatto,
 int deg_Gauss,
 int dim,
 double* jac_Gauss,
 double* M
)
{
  int volume_nodes_Lobatto = dgmath_get_nodes(dim,deg_Lobatto);
  
  double* u = P4EST_ALLOC_ZERO(double, volume_nodes_Lobatto);
  double* Mu = P4EST_ALLOC_ZERO(double, volume_nodes_Lobatto);

  for (int i = 0; i < volume_nodes_Lobatto; i++){
    u[i] = 1.;
    dgmath_apply_curvedGaussMass(dgmath_jit_dbase, u, deg_Lobatto, jac_Gauss, deg_Gauss, dim, Mu);
    linalg_set_column(M, Mu, i, volume_nodes_Lobatto, volume_nodes_Lobatto);
    u[i] = 0.;
  }

  P4EST_FREE(Mu);
  P4EST_FREE(u);
}

/* void dgmath_apply_hp_prolong(dgmath_jit_dbase_t* dgbase, double* in, int degH, */
                             /* int dim, int* degh, double* out) { */

/* void dgmath_apply_p_prolong(dgmath_jit_dbase_t* dgbase, double* in, int degH, */
                            /* int dim, int degh, double* out) { */
/* static double* dgmath_fetch_hp_prolong_transpose_1d(dgmath_jit_dbase_t* dgbase, */
                                                    /* int degH, int degh) { */
/* P will be \sum_i((degh_i+1)*(degH+1))^3 size */



void dgmath_compute_prolong_matrix
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 int degH,
 int dim,
 int* degh,
 int children,
 double* prolong_mat
)
{
  int volume_nodes_h = 0;
  for (int i = 0; i < children; i++){
    volume_nodes_h += dgmath_get_nodes(dim, degh[i]);
  }
  int volume_nodes_H = dgmath_get_nodes(dim, degH);

  double* u = P4EST_ALLOC_ZERO(double, volume_nodes_H);
  double* Mu = P4EST_ALLOC_ZERO(double, volume_nodes_h);

  for (int i = 0; i < volume_nodes_H; i++){
    u[i] = 1.;
    if (children == (P4EST_CHILDREN)){
      dgmath_apply_hp_prolong(dgmath_jit_dbase, u, degH, dim, degh, Mu);
    }
    else {
      dgmath_apply_p_prolong(dgmath_jit_dbase, u, degH, dim, degh[0], Mu);
    }
    linalg_set_column(prolong_mat, Mu, i, volume_nodes_h, volume_nodes_H);
    u[i] = 0.;
  }

  P4EST_FREE(Mu);
  P4EST_FREE(u);
}


void dgmath_compute_PT_mat_P
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* in_mat, /* dimensions of the degH space (degH+1)^DIM x (degH+1)^DIM */
 int degH,
 int dim,
 int* degh,
 int children,
 double* out_mat /* dimensions of the degH space: (degH+1)^DIM x (degH+1)^DIM */
)
{
  int volume_nodes_h = 0;
  for (int i = 0; i < children; i++){
    volume_nodes_h += dgmath_get_nodes(dim, degh[i]);
  }
  int volume_nodes_H = dgmath_get_nodes(dim, degH);
  
  double* promat = P4EST_ALLOC(double, volume_nodes_h*volume_nodes_H);
  double* in_mat_x_promat = P4EST_ALLOC(double, volume_nodes_h*volume_nodes_H);
  double* proTmat = P4EST_ALLOC(double, volume_nodes_h*volume_nodes_H);

  dgmath_compute_prolong_matrix(dgmath_jit_dbase, degH, dim, degh, children, promat);
  linalg_mat_transpose_nonsqr(promat, proTmat, volume_nodes_h, volume_nodes_H);
  linalg_mat_multiply(in_mat, promat, in_mat_x_promat, volume_nodes_h, volume_nodes_h, volume_nodes_H);
  linalg_mat_multiply(proTmat, in_mat_x_promat, out_mat, volume_nodes_H, volume_nodes_h, volume_nodes_H);

  P4EST_FREE(in_mat_x_promat);
  P4EST_FREE(promat);
  P4EST_FREE(proTmat);
}

void dgmath_form_fofufofvlilj_matrix_Gaussnodes
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* u,
 double* v,
 int deg_Lobatto,
 double* xyz_Gauss [(P4EST_DIM)],
 double* jac_Gauss,
 int deg_Gauss,
 int dim,
 double* out,
 grid_fcn_ext_t fofu_fcn,
 void* fofu_ctx,
 grid_fcn_ext_t fofv_fcn,
 void* fofv_ctx
)
{
  if (fofu_fcn == NULL)
    fofu_fcn = identity_fcn;
  if (fofv_fcn == NULL)
    fofv_fcn = identity_fcn;
  
  double* u_Gauss = NULL;
  double* v_Gauss = NULL;

  int volume_nodes_Gauss = dgmath_get_nodes(dim, deg_Gauss);
  
  if (u != NULL){
    u_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase,u,deg_Lobatto,deg_Gauss,u_Gauss,dim);
  }
  if (v != NULL){
    v_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase,v,deg_Lobatto,deg_Gauss,v_Gauss,dim);
  }
  
  double* fofu_fofv_jac = P4EST_ALLOC(double, volume_nodes_Gauss);
  for (int i = 0; i < volume_nodes_Gauss; i++){
    fofu_fofv_jac[i] = jac_Gauss[i];
    if (u != NULL || fofu_fcn != identity_fcn){
      fofu_fofv_jac[i] *= fofu_fcn(xyz_Gauss[0][i],
                                   xyz_Gauss[1][i],
#if (P4EST_DIM)==3
                                   xyz_Gauss[2][i],
#endif
                                   u_Gauss[i],
                                   fofu_ctx);
    }
    if (v != NULL || fofv_fcn != identity_fcn){
      fofu_fofv_jac[i] *= fofv_fcn(xyz_Gauss[0][i],
                                   xyz_Gauss[1][i],
#if (P4EST_DIM)==3
                                   xyz_Gauss[2][i],
#endif
                                   v_Gauss[i],
                                   fofv_ctx);
    }
  }
  dgmath_compute_curvedGaussMass
    (
     dgmath_jit_dbase,
     deg_Lobatto,
     deg_Gauss,
     dim,
     fofu_fofv_jac,
     out
    );
  
  if (u != NULL){
    P4EST_FREE(u_Gauss);
  }
  if (v != NULL){
    P4EST_FREE(v_Gauss);
  }
  
}


void dgmath_apply_fofufofvlilj_Gaussnodes
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* vec,
 double* u,
 double* v,
 int deg_Lobatto,
 double* jac_Gauss,
 double* xyz_Gauss [(P4EST_DIM)],
 int deg_Gauss,
 int dim,
 double* Mvec,
 grid_fcn_ext_t fofu_fcn,
 void* fofu_ctx,
 grid_fcn_ext_t fofv_fcn,
 void* fofv_ctx
)
{
  if (fofu_fcn == NULL){
    fofu_fcn = identity_fcn;
  }
  if (fofv_fcn == NULL){
    fofv_fcn = identity_fcn;
  }
  
  double* u_Gauss = NULL;
  double* v_Gauss = NULL;

  int volume_nodes_Gauss = dgmath_get_nodes(dim, deg_Gauss);
  
  if (u != NULL){
    u_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase,u,deg_Lobatto,deg_Gauss,u_Gauss,dim);
  }
  if (v != NULL){
    v_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase,v,deg_Lobatto,deg_Gauss,v_Gauss,dim);
  }
  
  double* fofu_fofv_jac = P4EST_ALLOC(double, volume_nodes_Gauss);
  for (int i = 0; i < volume_nodes_Gauss; i++){
    fofu_fofv_jac[i] = jac_Gauss[i];
    if (u != NULL || fofu_fcn != identity_fcn){
      fofu_fofv_jac[i] *= fofu_fcn(xyz_Gauss[0][i],
                                   xyz_Gauss[1][i],
#if (P4EST_DIM)==3
                                   xyz_Gauss[2][i],
#endif
                                   u_Gauss[i],
                                   fofu_ctx);
    }
    if (v != NULL || fofv_fcn != identity_fcn){
      fofu_fofv_jac[i] *= fofv_fcn(xyz_Gauss[0][i],
                                   xyz_Gauss[1][i],
#if (P4EST_DIM)==3
                                   xyz_Gauss[2][i],
#endif
                                   v_Gauss[i],
                                   fofv_ctx);
    }
  }

  dgmath_apply_curvedGaussMass
    (
     dgmath_jit_dbase,
     vec,
     deg_Lobatto,
     fofu_fofv_jac,
     deg_Gauss,
     dim,
     Mvec
    );
  
  P4EST_FREE(u_Gauss);
  P4EST_FREE(v_Gauss);
  P4EST_FREE(fofu_fofv_jac);
  
}

void dgmath_apply_fofufofvlj_Gaussnodes
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* u,
 double* v,
 int deg_Lobatto,
 double* jac_Gauss,
 double* xyz_Gauss [(P4EST_DIM)],
 int deg_Gauss,
 int dim,
 double* out,
 grid_fcn_ext_t fofu_fcn,
 void* fofu_ctx,
 grid_fcn_ext_t fofv_fcn,
 void* fofv_ctx
)
{
  if (fofu_fcn == NULL){
    fofu_fcn = identity_fcn;
  }
  if (fofv_fcn == NULL){
    fofv_fcn = identity_fcn;
  }
  
  double* u_Gauss = NULL;
  double* v_Gauss = NULL;

  int volume_nodes_Gauss = dgmath_get_nodes(dim, deg_Gauss);
  
  if (u != NULL){
    u_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase,u,deg_Lobatto,deg_Gauss,u_Gauss,dim);
  }
  if (v != NULL){
    v_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase,v,deg_Lobatto,deg_Gauss,v_Gauss,dim);
  }
  
  double* fofu_fofv = P4EST_ALLOC(double, volume_nodes_Gauss);
  for (int i = 0; i < volume_nodes_Gauss; i++){
    fofu_fofv[i] = 1.0;
    if (u != NULL || fofu_fcn != identity_fcn){
      fofu_fofv[i] *= fofu_fcn(xyz_Gauss[0][i],
                                   xyz_Gauss[1][i],
#if (P4EST_DIM)==3
                                   xyz_Gauss[2][i],
#endif
                                   u_Gauss[i],
                                   fofu_ctx);
    }
    if (v != NULL || fofv_fcn != identity_fcn){
      fofu_fofv[i] *= fofv_fcn(xyz_Gauss[0][i],
                                   xyz_Gauss[1][i],
#if (P4EST_DIM)==3
                                   xyz_Gauss[2][i],
#endif
                                   v_Gauss[i],
                                   fofv_ctx);
    }
  }

  dgmath_apply_curvedGaussMass_onGaussNodeVec
    (
     dgmath_jit_dbase,
     fofu_fofv,
     deg_Lobatto,
     jac_Gauss,
     deg_Gauss,
     dim,
     out
    );
  
  P4EST_FREE(u_Gauss);
  P4EST_FREE(v_Gauss);
  P4EST_FREE(fofu_fofv);
  
}


void dgmath_apply_curvedInverseGaussMass
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* in,
 int deg_Lobatto,
 double* jac_Gauss,
 int deg_Gauss,
 int dim,
 double* out
){
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int volume_nodes_Gauss = dgmath_get_nodes(dim, deg_Gauss);

  double* in_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
  double* one_over_w_j_in_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
  
  double* Gauss_weights = dgmath_fetch_GL_weights_1d(dgmath_jit_dbase, deg_Gauss);
  double* GLL_to_GL_interp_trans_inverse = dgmath_fetch_ref_GLL_to_GL_interp_trans_1d_inverse(dgmath_jit_dbase, deg_Lobatto, deg_Gauss);
  double* GLL_to_GL_interp_inverse = dgmath_fetch_ref_GLL_to_GL_interp_1d_inverse(dgmath_jit_dbase, deg_Lobatto, deg_Gauss);
  
  int nodes_Lobatto = deg_Lobatto + 1;
  int nodes_Gauss = deg_Gauss + 1;
  
  if (dim == 1){
    linalg_matvec_plus_vec(1.0, GLL_to_GL_interp_trans_inverse, in, 0., in_Gauss, nodes_Gauss, nodes_Lobatto);
    linalg_kron_oneover_vec_dot_oneover_x_dot_y(Gauss_weights, jac_Gauss, in_Gauss, nodes_Gauss, one_over_w_j_in_Gauss);
    linalg_matvec_plus_vec(1.0, GLL_to_GL_interp_inverse, one_over_w_j_in_Gauss, 0., out, nodes_Lobatto, nodes_Gauss);
  }
  
  else if (dim == 2){
    linalg_kron_A1A2x_NONSQR
      (
       in_Gauss,
       GLL_to_GL_interp_trans_inverse,
       GLL_to_GL_interp_trans_inverse,
       in,
       nodes_Gauss,
       nodes_Lobatto,
       nodes_Gauss,
       nodes_Lobatto
      );
    
    linalg_kron_oneover_vec_o_vec_dot_oneover_x_dot_y
      (
       Gauss_weights,
       jac_Gauss,
       in_Gauss,
       nodes_Gauss,
       one_over_w_j_in_Gauss
      );
    
    linalg_kron_A1A2x_NONSQR
      (
       out,
       GLL_to_GL_interp_inverse,
       GLL_to_GL_interp_inverse,
       one_over_w_j_in_Gauss,
       nodes_Lobatto,
       nodes_Gauss,
       nodes_Lobatto,
       nodes_Gauss
      );
  }
  
  else if (dim == 3) {
    linalg_kron_A1A2A3x_NONSQR(
                               in_Gauss,
                               GLL_to_GL_interp_trans_inverse,
                               GLL_to_GL_interp_trans_inverse,
                               GLL_to_GL_interp_trans_inverse,
                               in,
                               nodes_Gauss,
                               nodes_Lobatto,
                               nodes_Gauss,
                               nodes_Lobatto,
                               nodes_Gauss,
                               nodes_Lobatto
                              );
    
    linalg_kron_oneover_vec_o_vec_o_vec_dot_oneover_x_dot_y
      (
       Gauss_weights,
       jac_Gauss,
       in_Gauss,
       nodes_Gauss,
       one_over_w_j_in_Gauss
      );
    
    linalg_kron_A1A2A3x_NONSQR(
                               out,
                               GLL_to_GL_interp_inverse,
                               GLL_to_GL_interp_inverse,
                               GLL_to_GL_interp_inverse,
                               one_over_w_j_in_Gauss,
                               nodes_Lobatto,
                               nodes_Gauss,
                               nodes_Lobatto,
                               nodes_Gauss,
                               nodes_Lobatto,
                               nodes_Gauss
                              );
  } else {
    P4EST_FREE(one_over_w_j_in_Gauss);
    P4EST_FREE(in_Gauss);
    mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension.");
  }  
  P4EST_FREE(one_over_w_j_in_Gauss);
  P4EST_FREE(in_Gauss);
}



void dgmath_apply_curvedInverseLobattoMass
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* in,
 int deg_Lobatto,
 double* jac_Lobatto_integ,
 int deg_Lobatto_integ,
 int dim,
 double* out
){
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int volume_nodes_Lobatto_integ = dgmath_get_nodes(dim, deg_Lobatto_integ);

  double* in_Lobatto_integ = P4EST_ALLOC(double, volume_nodes_Lobatto_integ);
  double* one_over_w_j_in_Lobatto_integ = P4EST_ALLOC(double, volume_nodes_Lobatto_integ);
  
  double* Lobatto_integ_weights = dgmath_fetch_GLL_weights_1d(dgmath_jit_dbase, deg_Lobatto_integ);
  double* p_prolong_trans_inverse = dgmath_fetch_p_prolong_transpose_1d_inverse(dgmath_jit_dbase, deg_Lobatto, deg_Lobatto_integ);
  double* p_prolong_inverse = dgmath_fetch_p_prolong_1d_inverse(dgmath_jit_dbase, deg_Lobatto, deg_Lobatto_integ);
  
  int nodes_Lobatto = deg_Lobatto + 1;
  int nodes_Lobatto_integ = deg_Lobatto_integ + 1;
  
  if (dim == 1){
    linalg_matvec_plus_vec(1.0, p_prolong_trans_inverse, in, 0., in_Lobatto_integ, nodes_Lobatto_integ, nodes_Lobatto);
    linalg_kron_oneover_vec_dot_oneover_x_dot_y(Lobatto_integ_weights, jac_Lobatto_integ, in_Lobatto_integ, nodes_Lobatto_integ, one_over_w_j_in_Lobatto_integ);
    linalg_matvec_plus_vec(1.0, p_prolong_inverse, one_over_w_j_in_Lobatto_integ, 0., out, nodes_Lobatto, nodes_Lobatto_integ);
  }
  
  else if (dim == 2){
    linalg_kron_A1A2x_NONSQR
      (
       in_Lobatto_integ,
       p_prolong_trans_inverse,
       p_prolong_trans_inverse,
       in,
       nodes_Lobatto_integ,
       nodes_Lobatto,
       nodes_Lobatto_integ,
       nodes_Lobatto
      );
    
    linalg_kron_oneover_vec_o_vec_dot_oneover_x_dot_y
      (
       Lobatto_integ_weights,
       jac_Lobatto_integ,
       in_Lobatto_integ,
       nodes_Lobatto_integ,
       one_over_w_j_in_Lobatto_integ
      );
    
    linalg_kron_A1A2x_NONSQR
      (
       out,
       p_prolong_inverse,
       p_prolong_inverse,
       one_over_w_j_in_Lobatto_integ,
       nodes_Lobatto,
       nodes_Lobatto_integ,
       nodes_Lobatto,
       nodes_Lobatto_integ
      );
  }
  
  else if (dim == 3) {
    linalg_kron_A1A2A3x_NONSQR(
                               in_Lobatto_integ,
                               p_prolong_trans_inverse,
                               p_prolong_trans_inverse,
                               p_prolong_trans_inverse,
                               in,
                               nodes_Lobatto_integ,
                               nodes_Lobatto,
                               nodes_Lobatto_integ,
                               nodes_Lobatto,
                               nodes_Lobatto_integ,
                               nodes_Lobatto
                              );
    
    linalg_kron_oneover_vec_o_vec_o_vec_dot_oneover_x_dot_y
      (
       Lobatto_integ_weights,
       jac_Lobatto_integ,
       in_Lobatto_integ,
       nodes_Lobatto_integ,
       one_over_w_j_in_Lobatto_integ
      );
    
    linalg_kron_A1A2A3x_NONSQR(
                               out,
                               p_prolong_inverse,
                               p_prolong_inverse,
                               p_prolong_inverse,
                               one_over_w_j_in_Lobatto_integ,
                               nodes_Lobatto,
                               nodes_Lobatto_integ,
                               nodes_Lobatto,
                               nodes_Lobatto_integ,
                               nodes_Lobatto,
                               nodes_Lobatto_integ
                              );
  } else {
    P4EST_FREE(one_over_w_j_in_Lobatto_integ);
    P4EST_FREE(in_Lobatto_integ);
    mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension.");
  }  
  P4EST_FREE(one_over_w_j_in_Lobatto_integ);
  P4EST_FREE(in_Lobatto_integ);
}


/* void dgmath_interp_to_GL_nodes */
/* ( */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase, */
/*  double* u_Lobatto_in, */
/*  int deg_Lobatto, */
/*  int deg_Gauss, */
/*  double* u_Gauss_out, */
/*  int dim */
/* ) */
/* { */
/*   mpi_assert(dim == 1 || dim == 2 || dim == 3); */
/*   double* invVij = dgmath_fetch_invVij_1d(dgmath_jit_dbase, deg_Lobatto); */
/*   double* GaussVij = dgmath_fetch_ref_GaussVij_1d(dgmath_jit_dbase, deg_Lobatto, deg_Gauss); */
/*   double* temp = P4EST_ALLOC(double, dgmath_get_nodes((dim), deg_Lobatto)); */

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




void dgmath_hp_apply_nd_restrict_with_ptr(double* uH, int degH, double* uh,
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
    int cx = dgmath_is_child_left_or_right(c, 0);
    int cy = dgmath_is_child_left_or_right(c, 1);

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
    int cx = dgmath_is_child_left_or_right(c, 0);
    int cy = dgmath_is_child_left_or_right(c, 1);
    int cz = dgmath_is_child_left_or_right(c, 2);
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
    mpi_abort("ERROR: dgmath_restrict_nd_U\n");
  }
}

/* double* dgmath_fetch_GLL_1d(dgmath_jit_dbase_t* dgbase, int deg) { */
/*   /\* printf("deg,dgbase->degmath_max_storage = %d,%d" ,deg,dgbase->dgmath_max_storage); *\/ */
/*   mpi_assert(deg < dgbase->dgmath_max_storage); */
/*   if (deg > dgbase->dgmath_max_degree_used) */
/*     dgbase->dgmath_max_degree_used = deg; */
/*   return &(dgmath_GLLnodes[deg - 1][0]); */
/* } */

/* TODO: put this into the database */


/* TODO: put this into the database */
static void dgmath_build_dVij_1d(dgmath_jit_dbase_t* dgbase, double* dVij_1d,
                                 int deg) {
  int i, j, rows, cols;
  double* lobatto_nodes = dgmath_fetch_GLL_nodes_1d(dgbase, deg);
  rows = cols = deg + 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      dVij_1d[i * cols + j] = dgmath_gradjacobi(lobatto_nodes[i], 0., 0., j);
}

/* TODO: put this into the database */
static void dgmath_build_dGaussVij_1d(dgmath_jit_dbase_t* dgbase, double* dVij_1d,
                                      int Lobatto_degree, int Gauss_degree) {

  double* Gauss_nodes = dgmath_fetch_GL_nodes_1d(dgbase, Gauss_degree);
  
  int rows = Gauss_degree + 1;
  int cols = Lobatto_degree + 1;
  
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      dVij_1d[i * cols + j] = dgmath_gradjacobi(Gauss_nodes[i], 0., 0., j);
}

void dgmath_build_Mij_1d(dgmath_jit_dbase_t* dgbase, double* Mij_1d,
                                int deg) {
  int edge_nodes = deg + 1;
  double* v1d = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);
  double* v1d_trans = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);
  dgmath_build_Vij_1d(dgbase, v1d, deg);
  linalg_mat_transpose(v1d, v1d_trans, edge_nodes);
  linalg_mat_multiply(v1d, v1d_trans, Mij_1d, edge_nodes, edge_nodes,
                      edge_nodes);
  linalg_invert(Mij_1d, edge_nodes);
  P4EST_FREE(v1d_trans);
  P4EST_FREE(v1d);
}

void
dgmath_build_GLL_nodes_and_weights_1d(dgmath_jit_dbase_t* dgmath_jit_dbase, double* GLL_nodes, double* GLL_weights, int deg)
{
  dgmath_GLL_nodes_and_weights(deg+1, GLL_nodes, GLL_weights);
}


double* dgmath_fetch_GLL_nodes_1d(dgmath_jit_dbase_t* dgbase, int deg) {
  mpi_assert(deg < dgbase->dgmath_max_storage);
  if (dgbase->dgmath_GLL_nodes_1d_table[deg] != NULL) {
    return dgbase->dgmath_GLL_nodes_1d_table[deg];
  }
  else {
    int size = (deg + 1);
    dgbase->dgmath_GLL_nodes_1d_table[deg] = P4EST_ALLOC(double, size);
    dgbase->dgmath_GLL_weights_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_GLL_nodes_1d_table[deg];
    double* op1 = dgbase->dgmath_GLL_weights_1d_table[deg];
    dgmath_build_GLL_nodes_and_weights_1d(dgbase, op, op1, deg);
    return op;
  }
}

double* dgmath_fetch_GLL_weights_1d(dgmath_jit_dbase_t* dgbase, int deg) {
  mpi_assert(deg < dgbase->dgmath_max_storage);
  if (dgbase->dgmath_GLL_weights_1d_table[deg] != NULL) {
    return dgbase->dgmath_GLL_weights_1d_table[deg];
  }
  else {
    int size = (deg + 1);
    dgbase->dgmath_GLL_nodes_1d_table[deg] = P4EST_ALLOC(double, size);
    dgbase->dgmath_GLL_weights_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_GLL_nodes_1d_table[deg];
    double* op1 = dgbase->dgmath_GLL_weights_1d_table[deg];
    dgmath_build_GLL_nodes_and_weights_1d(dgbase, op, op1, deg);
    return op1;
  }
}

void
dgmath_build_GL_nodes_and_weights_1d(dgmath_jit_dbase_t* dgmath_jit_dbase, double* GL_nodes, double* GL_weights, int deg)
{
  dgmath_GL_nodes_and_weights(deg+1, GL_nodes, GL_weights);
}

double* dgmath_fetch_GL_nodes_1d(dgmath_jit_dbase_t* dgbase, int deg) {
  mpi_assert(deg < dgbase->dgmath_max_storage);

  if (dgbase->dgmath_GL_nodes_1d_table[deg] != NULL) {
    return dgbase->dgmath_GL_nodes_1d_table[deg];
  }
  else {
    int size = (deg + 1);
    dgbase->dgmath_GL_nodes_1d_table[deg] = P4EST_ALLOC(double, size);
    dgbase->dgmath_GL_weights_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_GL_nodes_1d_table[deg];
    double* op1 = dgbase->dgmath_GL_weights_1d_table[deg];
    dgmath_build_GL_nodes_and_weights_1d(dgbase, op, op1, deg);
    return op;
  }
}

double* dgmath_fetch_GL_weights_1d(dgmath_jit_dbase_t* dgbase, int deg) {
  mpi_assert(deg < dgbase->dgmath_max_storage);
  if (dgbase->dgmath_GL_weights_1d_table[deg] != NULL) {
    return dgbase->dgmath_GL_weights_1d_table[deg];
  }
  else {
    int size = (deg + 1);
    dgbase->dgmath_GL_nodes_1d_table[deg] = P4EST_ALLOC(double, size);
    dgbase->dgmath_GL_weights_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_GL_nodes_1d_table[deg];
    double* op1 = dgbase->dgmath_GL_weights_1d_table[deg];
    dgmath_build_GL_nodes_and_weights_1d(dgbase, op, op1, deg);
    return op1;
  }
}



/* double* dgmath_fetch_Mij_1d(dgmath_jit_dbase_t* dgbase, int deg) { */
/*   mpi_assert(deg < dgbase->dgmath_max_storage); */
/*   int* last_stride = &(dgbase->dgmath_ref_Mij_1d_table_last_stride); */
/*   if (dgbase->dgmath_ref_Mij_1d_table[deg] != -1) { */
/*     return &(dgbase->dgmath_ref_Mij_1d[dgbase->dgmath_ref_Mij_1d_table[deg]]); */
/*   } */

/*   else { */
/*     int size = (deg + 1) * (deg + 1); */
/*     dgbase->dgmath_ref_Mij_1d = */
/*         P4EST_REALLOC(dgbase->dgmath_ref_Mij_1d, double, *last_stride + size); */
/*     dgmath_build_Mij_1d(dgbase, &(dgbase->dgmath_ref_Mij_1d)[*last_stride], */
/*                         deg); */
/*     dgbase->dgmath_ref_Mij_1d_table[deg] = *last_stride; */
/*     *last_stride += size; */
/*     return &(dgbase->dgmath_ref_Mij_1d)[*last_stride - size]; */
/*   } */
/* } */

double* dgmath_fetch_Mij_1d(dgmath_jit_dbase_t* dgbase, int deg) {
  mpi_assert(deg < dgbase->dgmath_max_storage);
  if (dgbase->dgmath_ref_Mij_1d_table[deg] != NULL) {
    return dgbase->dgmath_ref_Mij_1d_table[deg];
  }
  else {
    int size = (deg + 1) * (deg + 1);
    dgbase->dgmath_ref_Mij_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_ref_Mij_1d_table[deg];
    dgmath_build_Mij_1d(dgbase, op, deg);
    return op;
  }
}




static void dgmath_build_invMij_1d(dgmath_jit_dbase_t* dgbase,
                                   double* invMij_1d, int deg) {
  double* Mij_1d = dgmath_fetch_Mij_1d(dgbase, deg);
  linalg_invert_and_copy(Mij_1d, invMij_1d, deg + 1);
}

static void dgmath_build_Dij_1d(dgmath_jit_dbase_t* dgbase, double* Dij_1d,
                                int deg) {
  int edge_nodes = deg + 1;
  double* v1d = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);
  double* inv_v1d = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);
  double* v1dr = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);

  dgmath_build_Vij_1d(dgbase, v1d, deg);
  dgmath_build_dVij_1d(dgbase, v1dr, deg);
  linalg_invert_and_copy(v1d, inv_v1d, edge_nodes);

  linalg_mat_multiply(v1dr, inv_v1d, Dij_1d, edge_nodes, edge_nodes,
                      edge_nodes);

  P4EST_FREE(v1dr);
  P4EST_FREE(inv_v1d);
  P4EST_FREE(v1d);
}


static void dgmath_build_GaussDij_1d(dgmath_jit_dbase_t* dgbase, double* Dij_1d,
                                     int deg_Lobatto, int deg_Gauss) {
  int nodes_Lobatto = deg_Lobatto + 1;
  int nodes_Gauss = deg_Gauss + 1;
  
  double* inv_v1d = (double*)P4EST_ALLOC(double, nodes_Lobatto*nodes_Lobatto);
  double* v1dr = (double*)P4EST_ALLOC(double, nodes_Gauss*nodes_Lobatto);

  dgmath_build_Vij_1d(dgbase, inv_v1d, deg_Lobatto);
  linalg_invert(inv_v1d, nodes_Lobatto);
  
  dgmath_build_dGaussVij_1d(dgbase, v1dr, deg_Lobatto, deg_Gauss);

  linalg_mat_multiply(v1dr, inv_v1d, Dij_1d, nodes_Gauss, nodes_Lobatto,
                      nodes_Lobatto);

  P4EST_FREE(v1dr);
  P4EST_FREE(inv_v1d);
}

static double* dgmath_fetch_invMij_1d(dgmath_jit_dbase_t* dgbase, int deg){
  mpi_assert(deg < dgbase->dgmath_max_storage);
  if (dgbase->dgmath_ref_invMij_1d_table[deg] != NULL) {
    return dgbase->dgmath_ref_invMij_1d_table[deg];
  }
  else {
    int size = (deg + 1) * (deg + 1);
    dgbase->dgmath_ref_invMij_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_ref_invMij_1d_table[deg];
    dgmath_build_invMij_1d(dgbase, op, deg);
    return op;
  }
}

void dgmath_apply_Mij(dgmath_jit_dbase_t* dgbase, double* in, int dim, int deg,
                      double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  double* mass_ref_1d = dgmath_fetch_Mij_1d(dgbase, deg);

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

void dgmath_apply_invMij(dgmath_jit_dbase_t* dgbase, double* in, int dim,
                         int deg, double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  double* inv_mass_ref_1d = dgmath_fetch_invMij_1d(dgbase, deg);

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

/* static double* dgmath_fetch_Dij_1d(dgmath_jit_dbase_t* dgbase, int deg) { */
/*   mpi_assert(deg < dgbase->dgmath_max_storage); */
/*   int* last_stride = &(dgbase->dgmath_ref_Dij_1d_table_last_stride); */
/*   if (dgbase->dgmath_ref_Dij_1d_table[deg] != -1) { */
/*     return &(dgbase->dgmath_ref_Dij_1d[dgbase->dgmath_ref_Dij_1d_table[deg]]); */
/*   } else { */
/*     int size = (deg + 1) * (deg + 1); */
/*     dgbase->dgmath_ref_Dij_1d = */
/*         P4EST_REALLOC(dgbase->dgmath_ref_Dij_1d, double, *last_stride + size); */

/*     dgmath_build_Dij_1d(dgbase, &(dgbase->dgmath_ref_Dij_1d)[*last_stride], */
/*                         deg); */
/*     dgbase->dgmath_ref_Dij_1d_table[deg] = *last_stride; */
/*     *last_stride += size; */
/*     return &(dgbase->dgmath_ref_Dij_1d)[*last_stride - size]; */
/*   } */
/* } */


static double* dgmath_fetch_Dij_1d(dgmath_jit_dbase_t* dgbase, int deg) {
  mpi_assert(deg < dgbase->dgmath_max_storage);
  if (dgbase->dgmath_ref_Dij_1d_table[deg] != NULL) {
    return dgbase->dgmath_ref_Dij_1d_table[deg];
  }
  else {
    int size = (deg + 1) * (deg + 1);
    dgbase->dgmath_ref_Dij_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_ref_Dij_1d_table[deg];
    dgmath_build_Dij_1d(dgbase, op, deg);
    return op;
  }
}

double* dgmath_fetch_GaussDij_1d(dgmath_jit_dbase_t* dgbase, int deg_Lobatto,
                                         int deg_Gauss) {
  mpi_assert(deg_Gauss <= dgbase->dgmath_max_storage);
  mpi_assert(deg_Lobatto <= dgbase->dgmath_max_storage);

  if (dgbase->dgmath_ref_GaussDij_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] != NULL) {
    return dgbase->dgmath_ref_GaussDij_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss];
  }
  else {
    int size = (deg_Gauss + 1) * (deg_Lobatto + 1);
    dgbase->dgmath_ref_GaussDij_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_ref_GaussDij_1d_table[dgbase->dgmath_max_storage * deg_Lobatto + deg_Gauss];
    dgmath_build_GaussDij_1d(dgbase, op, deg_Lobatto, deg_Gauss);
    return op;
  }
}

  void
    dgmath_apply_GaussDij(dgmath_jit_dbase_t* dgbase, double* in, int dim, int deg_Lobatto, int deg_Gauss, int dir, double* out)
{
  /* printf("dir = %d, dim = %d\n", dir, dim); */

    mpi_assert(dir < dim && dir >= 0 && (dim == 1 || dim == 3 || dim == 2));

    /* TODO get rid of this */
  /* mpi_assert(deg_Lobatto == deg_Gauss); */
  
  double* Dr_1d = dgmath_fetch_GaussDij_1d(dgbase, deg_Lobatto, deg_Gauss);
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
    mpi_abort("ERROR: dgmath_apply_GaussDij");
  }
}

void dgmath_apply_Dij(dgmath_jit_dbase_t* dgbase, double* in, int dim, int deg,
                      int dir, double* out) {
  mpi_assert(dir < dim && dir >= 0 && (dim == 1 || dim == 3 || dim == 2));
  double* Dr_1d = dgmath_fetch_Dij_1d(dgbase, deg);
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
    mpi_abort("ERROR: dgmath_apply_Dij");
  }
}


void dgmath_apply_Dij_transpose(dgmath_jit_dbase_t* dgbase, double* in, int dim, int deg,
                      int dir, double* out) {
  mpi_assert(dir < dim && dir >= 0 && (dim == 1 || dim == 3 || dim == 2));
  double* Dr_1d = dgmath_fetch_Dij_1d(dgbase, deg);

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
    mpi_abort("ERROR: dgmath_apply_Dij");
  }

  P4EST_FREE(Dr_1d_transpose);
}

static void dgmath_build_hp_prolong_1d_aux(int degH, int degh, int c,
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
    if (c == 0 || c == 1) dgmath_child_to_parent_ref_coords(c, &r);
    for (i = 0; i < nodesH; i++) {
      phi_modal[i] = dgmath_jacobi(r, 0, 0, i);
    }

    linalg_matvec_plus_vec(1.0, inv_v1d_trans_degH, phi_modal, 0.0,
                           &hp_prolong_matrix_1d[n * nodesH], nodesH, nodesH);
  }
  P4EST_FREE(phi_modal);
}

static void dgmath_build_hp_prolong_1d(dgmath_jit_dbase_t* dgbase,
                                       double* hp_prolong_1d, int degH,
                                       int degh) {
  int edge_nodes_degH = degH + 1;
  int edge_nodes_degh = degh + 1;

  double* v1d = (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH);
  double* inv_v1d =
      (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH);
  double* inv_v1d_trans_degH =
      (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH);

  dgmath_build_Vij_1d(dgbase, v1d, degH);
  linalg_invert_and_copy(v1d, inv_v1d, edge_nodes_degH);
  linalg_mat_transpose(inv_v1d, inv_v1d_trans_degH, edge_nodes_degH);
  double* lobatto_nodes_degh = dgmath_fetch_GLL_nodes_1d(dgbase, degh);

  dgmath_build_hp_prolong_1d_aux(degH, degh, 0, inv_v1d_trans_degH,
                                 lobatto_nodes_degh, &hp_prolong_1d[0]);
  dgmath_build_hp_prolong_1d_aux(
      degH, degh, 1, inv_v1d_trans_degH, lobatto_nodes_degh,
      &hp_prolong_1d[edge_nodes_degH * edge_nodes_degh]);

  P4EST_FREE(inv_v1d);
  P4EST_FREE(inv_v1d_trans_degH);
  P4EST_FREE(v1d);
}

void dgmath_build_p_prolong_1d(dgmath_jit_dbase_t* dgbase,
                                      double* p_prolong_1d, int degH,
                                      int degh) {
  /* int edge_nodes_degH = degH + 1; */
  /* /\* int edge_nodes_degh = degh + 1; *\/ */

  /* double* v1d = (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH); */
  /* double* inv_v1d = */
  /*     (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH); */
  /* double* inv_v1d_trans_degH = */
  /*     (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH); */

  /* dgmath_build_Vij_1d(dgbase, v1d, degH); */
  /* linalg_invert_and_copy(v1d, inv_v1d, edge_nodes_degH); */
  /* linalg_mat_transpose(inv_v1d, inv_v1d_trans_degH, edge_nodes_degH); */
  /* double* lobatto_nodes_degh = dgmath_fetch_GLL_nodes_1d(dgbase, degh); */

  /* dgmath_build_hp_prolong_1d_aux(degH, degh, -1, inv_v1d_trans_degH, */
  /*                                lobatto_nodes_degh, &p_prolong_1d[0]); */

  /* P4EST_FREE(inv_v1d); */
  /* P4EST_FREE(inv_v1d_trans_degH); */
  /* P4EST_FREE(v1d); */
  double* degh_nodes = dgmath_fetch_GLL_nodes_1d(dgbase, degh);
  double* invVij = dgmath_fetch_invVij_1d(dgbase, degH);
  double* Vij_degh = P4EST_ALLOC(double, (degH+1)*(degh+1));
  
  int rows = degh + 1;
  int cols = degH + 1;
  
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      Vij_degh[i * cols + j] = dgmath_jacobi(degh_nodes[i], 0., 0., j);

  linalg_mat_multiply(Vij_degh, invVij, p_prolong_1d, degh+1, degH+1, degH+1);
  P4EST_FREE(Vij_degh);
}




void dgmath_hp_apply_nd_prolong_transpose_with_ptr(
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
    int cx = dgmath_is_child_left_or_right(c, 0);
    int cy = dgmath_is_child_left_or_right(c, 1);
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
    int cx = dgmath_is_child_left_or_right(c, 0);
    int cy = dgmath_is_child_left_or_right(c, 1);
    int cz = dgmath_is_child_left_or_right(c, 2);
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
    mpi_abort("ERROR: dgmath_prolong_transpose_nd_U\n");
  }
}

static double* dgmath_fetch_hp_prolong_1d(dgmath_jit_dbase_t* dgbase, int degH,
                                          int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < dgbase->dgmath_max_storage);
  
  if (dgbase->dgmath_hp_prolong_1d_table[dgbase->dgmath_max_storage * degH + degh] != NULL) {
    return dgbase->dgmath_hp_prolong_1d_table[dgbase->dgmath_max_storage * degH + degh];
  }
  else {
    int size = 2 * (degh + 1) * (degH + 1);
    dgbase->dgmath_hp_prolong_1d_table[dgbase->dgmath_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_hp_prolong_1d_table[dgbase->dgmath_max_storage * degH + degh];
    dgmath_build_hp_prolong_1d(dgbase, op, degH, degh);
    return op;
  }
}

static double* dgmath_fetch_p_prolong_1d(dgmath_jit_dbase_t* dgbase, int degH,
                                         int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh <= dgbase->dgmath_max_storage);

  if (dgbase->dgmath_p_prolong_1d_table[dgbase->dgmath_max_storage * degH + degh] != NULL) {
    return dgbase->dgmath_p_prolong_1d_table[dgbase->dgmath_max_storage * degH + degh];
  }
  else {
    int size = (degh + 1) * (degH + 1);
    dgbase->dgmath_p_prolong_1d_table[dgbase->dgmath_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_p_prolong_1d_table[dgbase->dgmath_max_storage * degH + degh];
    dgmath_build_p_prolong_1d(dgbase, op, degH, degh);
    return op;
  }
}

double* dgmath_fetch_p_prolong_1d_inverse(dgmath_jit_dbase_t* dgbase, int degH,
                                         int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh <= dgbase->dgmath_max_storage);

  if (dgbase->dgmath_p_prolong_1d_inverse_table[dgbase->dgmath_max_storage * degH + degh] != NULL) {
    return dgbase->dgmath_p_prolong_1d_inverse_table[dgbase->dgmath_max_storage * degH + degh];
  }
  else {
    int size = (degh + 1) * (degH + 1);
    dgbase->dgmath_p_prolong_1d_inverse_table[dgbase->dgmath_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_p_prolong_1d_inverse_table[dgbase->dgmath_max_storage * degH + degh];
    dgmath_build_p_prolong_1d_inverse(dgbase, op, degH, degh);
    return op;
  }
}


void dgmath_apply_hp_prolong(dgmath_jit_dbase_t* dgbase, double* in, int degH,
                             int dim, int* degh, double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int children = (1 << dim);
  int c;
  double* hp_prolong_1d;

  int stride = 0;
  for (c = 0; c < children; c++) {
    hp_prolong_1d = dgmath_fetch_hp_prolong_1d(dgbase, degH, degh[c]);
    dgmath_hp_apply_nd_prolong_with_ptr(&out[stride], degh[c], in, degH, dim, c,
                                        hp_prolong_1d);
    stride += dgmath_get_nodes(dim, degh[c]);
  }
}

void dgmath_apply_p_prolong(dgmath_jit_dbase_t* dgbase, double* in, int degH,
                            int dim, int degh, double* out) {
  double* p_prolong_1d = dgmath_fetch_p_prolong_1d(dgbase, degH, degh);

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  /* if degrees are the same, p_prolong matrix should be the identity */
  if (degh == degH) {
    linalg_copy_1st_to_2nd(in, out, dgmath_get_nodes(dim, degh));
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
    mpi_abort("ERROR: dgmath_prolong_nd_U\n");
  }
}

void dgmath_build_hp_restrict_1d_aux(int degh, int degH,
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

static void dgmath_build_p_restrict_1d(dgmath_jit_dbase_t* dgbase,
                                       double* p_restrict_1d, int degH,
                                       int degh) {
  int edge_nodes_degH = degH + 1;
  int edge_nodes_degh = degh + 1;

  /* printf("degH = %d\n",degH); */
  /* printf("degh = %d\n",degh); */

  double* p_prolong_1d = dgmath_fetch_p_prolong_1d(dgbase, degH, degh);
  double* inv_Mij_degH =
      (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH);

  /* util_print_matrix(Mij_degH, edge_nodes_degH, edge_nodes_degH, "Mij_degH =
   * ", 0); */
  /* util_print_matrix(inv_Mij_degH, edge_nodes_degH, edge_nodes_degH,
   * "Mij_degH_inv = ", 0); */
  double* Mij_degH = dgmath_fetch_Mij_1d(dgbase, degH);
  linalg_invert_and_copy(Mij_degH, inv_Mij_degH, edge_nodes_degH);

  double* Mij_degh = dgmath_fetch_Mij_1d(dgbase, degh);
  dgmath_build_hp_restrict_1d_aux(degh, degH, &p_prolong_1d[0], Mij_degh,
                                  inv_Mij_degH, &p_restrict_1d[0]);

  linalg_vec_scale(2., p_restrict_1d, edge_nodes_degH * edge_nodes_degh);

  /* util_print_matrix(p_restrict_1d, edge_nodes_degH, edge_nodes_degh,
   * "p_restrict_1d[0] = ", 0); */
  /* util_print_matrix(&p_restrict_1d[edge_nodes_degH*edge_nodes_degh],
   * edge_nodes_degH, edge_nodes_degh, "p_restrict_1d[1] = ", 0); */

  P4EST_FREE(inv_Mij_degH);
}

 double* dgmath_fetch_p_restrict_1d(dgmath_jit_dbase_t* dgbase, int degH,
                                          int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < dgbase->dgmath_max_storage);

  if (dgbase->dgmath_p_restrict_1d_table[dgbase->dgmath_max_storage * degH + degh] != NULL) {
    return dgbase->dgmath_p_restrict_1d_table[dgbase->dgmath_max_storage * degH + degh];
  }
  else {
    int size = (degh + 1) * (degH + 1);
    dgbase->dgmath_p_restrict_1d_table[dgbase->dgmath_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_p_restrict_1d_table[dgbase->dgmath_max_storage * degH + degh];
    dgmath_build_p_restrict_1d(dgbase, op, degH, degh);
    return op;
  }
}

void dgmath_apply_p_restrict(dgmath_jit_dbase_t* dgbase, double* in, int degh,
                             int dim, int degH, double* out) {
  double* p_restrict_1d = dgmath_fetch_p_restrict_1d(dgbase, degH, degh);

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  /* restrict matrix should be the identity if degrees are the same */
  if (degh == degH) {
    linalg_copy_1st_to_2nd(in, out, dgmath_get_nodes(dim, degh));
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
    mpi_abort("ERROR: dgmath_restrict_nd_U\n");
  }
}

static void dgmath_build_hp_restrict_1d(dgmath_jit_dbase_t* dgbase,
                                        double* hp_restrict_1d, int degH,
                                        int degh) {
  int edge_nodes_degH = degH + 1;
  int edge_nodes_degh = degh + 1;

  /* printf("degH = %d\n",degH); */
  /* printf("degh = %d\n",degh); */

  double* hp_prolong_1d = dgmath_fetch_hp_prolong_1d(dgbase, degH, degh);

  double* inv_Mij_degH =
      (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH);

  /* TODO: replace this with a call to the database for inverse Mij */
  double* Mij_degH = dgmath_fetch_Mij_1d(dgbase, degH);
  linalg_invert_and_copy(Mij_degH, inv_Mij_degH, edge_nodes_degH);

  double* Mij_degh = dgmath_fetch_Mij_1d(dgbase, degh);
  dgmath_build_hp_restrict_1d_aux(degh, degH, &hp_prolong_1d[0], Mij_degh,
                                  inv_Mij_degH, &hp_restrict_1d[0]);

  dgmath_build_hp_restrict_1d_aux(
      degh, degH, &hp_prolong_1d[edge_nodes_degH * edge_nodes_degh], Mij_degh,
      inv_Mij_degH, &hp_restrict_1d[edge_nodes_degH * edge_nodes_degh]);

  P4EST_FREE(inv_Mij_degH);
}

static double* dgmath_fetch_hp_restrict_1d(dgmath_jit_dbase_t* dgbase, int degH,
                                           int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < dgbase->dgmath_max_storage);

  if (dgbase->dgmath_hp_restrict_1d_table[dgbase->dgmath_max_storage * degH + degh] != NULL) {
    return dgbase->dgmath_hp_restrict_1d_table[dgbase->dgmath_max_storage * degH + degh];
  }
  else {
    int size = 2 * (degh + 1) * (degH + 1);
    dgbase->dgmath_hp_restrict_1d_table[dgbase->dgmath_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_hp_restrict_1d_table[dgbase->dgmath_max_storage * degH + degh];
    dgmath_build_hp_restrict_1d(dgbase, op, degH, degh);
    return op;
  }
}

void dgmath_apply_hp_restrict(dgmath_jit_dbase_t* dgbase, double* in, int* degh,
                              int dim, int degH, double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int children = (1 << dim);
  int c;
  int nodesH = dgmath_get_nodes(dim, degH);
  double* hp_restrict_1d;

  int stride = 0;
  double* tmp = P4EST_ALLOC(double, nodesH);
  linalg_fill_vec(out, 0., nodesH);

  for (c = 0; c < children; c++) {
    hp_restrict_1d = dgmath_fetch_hp_restrict_1d(dgbase, degH, degh[c]);
    /* int nodesh = dgmath_get_nodes(dim,degh[c]); */

    dgmath_hp_apply_nd_restrict_with_ptr(tmp, degH, &in[stride], degh[c], dim,
                                         c, hp_restrict_1d);
    /* printf("child = %d\n",c); */
    /* util_print_matrix(tmp, nodesH, 1, "tmp = ", 0); */
    /* util_print_matrix(&uh[stride], nodesh, 1, "uh[stride] = ", 0); */
    linalg_vec_axpy(1.0, tmp, out, nodesH);
    /* util_print_matrix(uH, nodesH, 1, "uH = ", 0); */
    stride += dgmath_get_nodes(dim, degh[c]);
  }

  P4EST_FREE(tmp);
}

int dgmath_fetch_max_degree_used(dgmath_jit_dbase_t* dgbase) {
  return dgbase->dgmath_max_degree_used;
}

static void dgmath_build_xyz_nd(dgmath_jit_dbase_t* dgbase, double* ref_xyz_nd,
                                int dim, int deg) {
  int nodes = deg + 1;
  int vol_nodes = dgmath_get_nodes(dim, deg);
  double* lgl = dgmath_fetch_GLL_nodes_1d(dgbase, deg);
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
    mpi_abort("ERROR 1: dgmath_build_xyz_nd");

  P4EST_FREE(eye);
}


static void dgmath_build_Gauss_xyz_nd(dgmath_jit_dbase_t* dgbase, double* ref_xyz_nd,
                                int dim, int deg) {
  int nodes = deg + 1;
  int vol_nodes = dgmath_get_nodes(dim, deg);
  double* gl = dgmath_fetch_GL_nodes_1d(dgbase, deg);
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
    mpi_abort("ERROR 1: dgmath_build_xyz_nd");

  P4EST_FREE(eye);
}

/* double* dgmath_fetch_xyz_nd(dgmath_jit_dbase_t* dgbase, int dim, int deg, */
/*                             int dir) { */
/*   /\* printf("dim,deg,dir=%d,%d,%d\n",dim,deg,dir); *\/ */
/*   mpi_assert(deg < dgbase->dgmath_max_storage && (dim == 3 || dim == 2) && */
/*              (dir >= 0 && dir < dim)); */
/*   int* last_stride = &(dgbase->dgmath_ref_xyz_nd_table_last_stride); */
/*   if (dgbase->dgmath_ref_xyz_nd_table[dim * dgbase->dgmath_max_storage + deg] != */
/*       -1) { */
/*     return &( */
/*         dgbase->dgmath_ref_xyz_nd[dgbase->dgmath_ref_xyz_nd_table */
/*                                       [dim * dgbase->dgmath_max_storage + deg] + */
/*                                   dir * dgmath_get_nodes(dim, deg)]); */
/*   } else { */
/*     int size = (P4EST_DIM)*dgmath_get_nodes(dim, deg); */
/*     dgbase->dgmath_ref_xyz_nd = */
/*         P4EST_REALLOC(dgbase->dgmath_ref_xyz_nd, double, *last_stride + size); */

/*     dgmath_build_xyz_nd(dgbase, &(dgbase->dgmath_ref_xyz_nd)[*last_stride], dim, */
/*                         deg); */
/*     dgbase->dgmath_ref_xyz_nd_table[dim * dgbase->dgmath_max_storage + deg] = */
/*         *last_stride; */
/*     *last_stride += size; */
/*     return &(dgbase->dgmath_ref_xyz_nd)[*last_stride - size + */
/*                                         dir * dgmath_get_nodes(dim, deg)]; */
/*   } */
/*   /\* return NULL; *\/ */
/* } */


double* dgmath_fetch_xyz_nd(dgmath_jit_dbase_t* dgbase, int dim, int deg,
                            int dir) {

  mpi_assert(deg < dgbase->dgmath_max_storage && (dim == 3 || dim == 2 || dim == 1) &&
             (dir >= 0 && dir < dim));

  
  if (dim == 1)
    return dgmath_fetch_GLL_nodes_1d(dgbase, deg);
  
  if (dgbase->dgmath_ref_xyz_nd_table[(dim-2) * dgbase->dgmath_max_storage + deg] != NULL){
    double* table = dgbase->dgmath_ref_xyz_nd_table[(dim-2) * dgbase->dgmath_max_storage + deg];
    return &table[dir * dgmath_get_nodes(dim, deg)];
  }
  
  else {
    int size = dim*dgmath_get_nodes(dim, deg);
    dgbase->dgmath_ref_xyz_nd_table[(dim-2)* dgbase->dgmath_max_storage + deg] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_ref_xyz_nd_table[(dim-2) * dgbase->dgmath_max_storage + deg];
    dgmath_build_xyz_nd(dgbase, op, dim, deg);
    return &op[dir * dgmath_get_nodes(dim, deg)];
  }
}



double* dgmath_fetch_Gauss_xyz_nd(dgmath_jit_dbase_t* dgbase, int dim, int deg,
                            int dir) {

  mpi_assert(deg < dgbase->dgmath_max_storage && (dim == 3 || dim == 2 || dim == 1) &&
             (dir >= 0 && dir < dim));

  
  if (dim == 1)
    return dgmath_fetch_GL_nodes_1d(dgbase, deg);
  
  if (dgbase->dgmath_ref_Gauss_xyz_nd_table[(dim-2) * dgbase->dgmath_max_storage + deg] != NULL){
    double* table = dgbase->dgmath_ref_Gauss_xyz_nd_table[(dim-2) * dgbase->dgmath_max_storage + deg];
    return &table[dir * dgmath_get_nodes(dim, deg)];
  }
  
  else {
    int size = dim*dgmath_get_nodes(dim, deg);
    dgbase->dgmath_ref_Gauss_xyz_nd_table[(dim-2)* dgbase->dgmath_max_storage + deg] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_ref_Gauss_xyz_nd_table[(dim-2) * dgbase->dgmath_max_storage + deg];
    dgmath_build_Gauss_xyz_nd(dgbase, op, dim, deg);
    return &op[dir * dgmath_get_nodes(dim, deg)];
  }
}


static void dgmath_build_LIFT_1d(double* LIFT_1d, int deg) {
  memset(LIFT_1d, 0., sizeof(double) * 2 * (deg + 1));
  LIFT_1d[0] = 1.;
  LIFT_1d[2 * deg + 1] = 1.;
}

/* static double* dgmath_fetch_LIFT_1d(dgmath_jit_dbase_t* dgbase, int deg) { */
  
/*   mpi_assert(deg < dgbase->dgmath_max_storage); */
/*   int* last_stride = &(dgbase->dgmath_LIFT_1d_table_last_stride); */

/*   if (dgbase->dgmath_LIFT_1d_table[deg] != -1) { */
/*     return &(dgbase->dgmath_LIFT_1d[dgbase->dgmath_LIFT_1d_table[deg]]); */
/*   } else { */
/*     int size = 2 * (deg + 1); */
/*     dgbase->dgmath_LIFT_1d = */
/*         P4EST_REALLOC(dgbase->dgmath_LIFT_1d, double, *last_stride + size); */
/*     dgmath_build_LIFT_1d(&(dgbase->dgmath_LIFT_1d)[*last_stride], deg); */
/*     dgbase->dgmath_LIFT_1d_table[deg] = *last_stride; */
/*     *last_stride += size; */
/*     return &(dgbase->dgmath_LIFT_1d)[*last_stride - size]; */
/*   } */
/* } */


static double* dgmath_fetch_LIFT_1d(dgmath_jit_dbase_t* dgbase, int deg) {

  mpi_assert(deg < dgbase->dgmath_max_storage);

  if (dgbase->dgmath_LIFT_1d_table[deg] != NULL) {
    return dgbase->dgmath_LIFT_1d_table[deg];
  }
  else {
    int size = 2 * (deg + 1);
    dgbase->dgmath_LIFT_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_LIFT_1d_table[deg];
    dgmath_build_LIFT_1d(op, deg);
    return op;
  }
}

void dgmath_apply_LIFT(dgmath_jit_dbase_t* dgbase, double* in, int dim, int deg,
                       int face, double* out) {
  mpi_assert(dim == 2 || dim == 3);
  mpi_assert(face < 2 * (dim));

  double* lift_1d = dgmath_fetch_LIFT_1d(dgbase, deg);

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
    mpi_abort("ERROR 0: dgmath_lift_boundary_vec");
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

void dgmath_apply_slicer(dgmath_jit_dbase_t* dgbase, double* in, int dim,
                         int face, int deg, double* out) {
  mpi_assert(face < 2 * (dim));
  /* int nodes = util_pow_int(deg+1, DIM); */

  double* slicer_1d = dgmath_fetch_LIFT_1d(dgbase, deg);

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
    mpi_abort("ERROR 0: dgmath_lift_boundary_vec");
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

void dgmath_get_normal(int face, int dim, double* n) {
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
    mpi_abort("ERROR: dgmath_get_normal");
  }
}

double dgmath_jacobi(double r, double alpha, double beta, int N) {
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

double dgmath_gradjacobi(double r, double alpha, double beta, int N) {
  if (N == 0)
    return 0;
  else
    return sqrt(N * (N + alpha + beta + 1.0)) *
           dgmath_jacobi(r, alpha + 1.0, beta + 1.0, N - 1);
}

void dgmath_vandermonde_1d(double* v1d, double* lobatto_nodes, int degree) {
  int i, j, rows, cols;
  rows = cols = degree + 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      v1d[i * cols + j] = dgmath_jacobi(lobatto_nodes[i], 0., 0., j);
}

void dgmath_grad_vandermonde_1d(double* dr_v1d, double* lobatto_nodes,
                                int degree) {
  int i, j, rows, cols;
  rows = cols = degree + 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      dr_v1d[i * cols + j] = dgmath_gradjacobi(lobatto_nodes[i], 0., 0., j);
}

void dgmath_convert_nodal_to_modal(dgmath_jit_dbase_t* dgbase, double* in,
                                   int dim, int deg, double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int nodes = deg + 1;
  double* invVij_1d = P4EST_ALLOC(double, nodes* nodes);

  /* TODO: probably could use build function for inv_Vij */
  dgmath_build_Vij_1d(dgbase, invVij_1d, deg);
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

double dgmath_rtox(double r, double xl, double h) {
  return h * (r + 1.) / 2. + xl;
}

void dgmath_rtox_array(double* r, double xl, double h, double* x, int nodes) {
  int i;
  for (i = 0; i < nodes; i++) {
    x[i] = dgmath_rtox(r[i], xl, h);
  }
}

static void dgmath_build_p_prolong_transpose_1d(dgmath_jit_dbase_t* dgbase,
                                                double* p_prolong_transpose_1d,
                                                int degH, int degh) {
  double* p_prolong_1d = dgmath_fetch_p_prolong_1d(dgbase, degH, degh);
  linalg_mat_transpose_nonsqr(p_prolong_1d, p_prolong_transpose_1d, (degh + 1),
                              (degH + 1));
}





static void dgmath_build_hp_prolong_transpose_1d(
    dgmath_jit_dbase_t* dgbase, double* hp_prolong_transpose_1d, int degH,
    int degh) {
  double* hp_prolong_1d = dgmath_fetch_hp_prolong_1d(dgbase, degH, degh);
  /* linalg_mat_transpose_nonsqr(hp_prolong_1d, hp_prolong_transpose_1d,
   * 2*(degh+1), (degH+1)); */
  linalg_mat_transpose_nonsqr(hp_prolong_1d, hp_prolong_transpose_1d,
                              (degh + 1), (degH + 1));
  linalg_mat_transpose_nonsqr(&hp_prolong_1d[(degh + 1) * (degH + 1)],
                              &hp_prolong_transpose_1d[(degh + 1) * (degH + 1)],
                              (degh + 1), (degH + 1));
}

static double* dgmath_fetch_hp_prolong_transpose_1d(dgmath_jit_dbase_t* dgbase,
                                                    int degH, int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < dgbase->dgmath_max_storage);
  if (dgbase->dgmath_hp_prolong_transpose_1d_table[dgbase->dgmath_max_storage *degH + degh] != NULL) {
    return dgbase->dgmath_hp_prolong_transpose_1d_table[dgbase->dgmath_max_storage *degH + degh];
  }
  else {
    int size = 2 * (degh + 1) * (degH + 1);
    dgbase->dgmath_hp_prolong_transpose_1d_table[dgbase->dgmath_max_storage *degH + degh] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_hp_prolong_transpose_1d_table[dgbase->dgmath_max_storage *degH + degh];
    dgmath_build_hp_prolong_transpose_1d(dgbase, op, degH, degh);
    return op;
  }
}





void dgmath_build_p_prolong_1d_inverse(dgmath_jit_dbase_t* dgbase,
                                      double* p_prolong_1d_inverse, int degH,
                                      int degh)
{
  double* p_prolong_1d = dgmath_fetch_p_prolong_1d(dgbase, degH, degh);
  linalg_leftinverse(p_prolong_1d, p_prolong_1d_inverse, degh + 1, degH + 1);
}





double* dgmath_fetch_p_prolong_transpose_1d(dgmath_jit_dbase_t* dgbase,
                                                   int degH, int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < dgbase->dgmath_max_storage);

  if (dgbase->dgmath_p_prolong_transpose_1d_table[dgbase->dgmath_max_storage *degH + degh] != NULL) {
    return dgbase->dgmath_p_prolong_transpose_1d_table[dgbase->dgmath_max_storage *degH + degh];
  }
  else {
    int size = (degh + 1) * (degH + 1);
    dgbase->dgmath_p_prolong_transpose_1d_table[dgbase->dgmath_max_storage *degH + degh] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_p_prolong_transpose_1d_table[dgbase->dgmath_max_storage *degH + degh];
    dgmath_build_p_prolong_transpose_1d(dgbase, op, degH, degh);
    return op;
  }
}


void dgmath_build_p_prolong_transpose_1d_inverse(dgmath_jit_dbase_t* dgbase,
                                                double* p_prolong_transpose_1d_inverse,
                                                int degH, int degh) {
  double* p_prolong_transpose_1d = dgmath_fetch_p_prolong_transpose_1d(dgbase, degH, degh);
  linalg_leftinverse(p_prolong_transpose_1d, p_prolong_transpose_1d_inverse, degH + 1, degh + 1);
}

 double* dgmath_fetch_p_prolong_transpose_1d_inverse(dgmath_jit_dbase_t* dgbase,
                                                   int degH, int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < dgbase->dgmath_max_storage);

  if (dgbase->dgmath_p_prolong_transpose_1d_inverse_table[dgbase->dgmath_max_storage *degH + degh] != NULL) {
    return dgbase->dgmath_p_prolong_transpose_1d_inverse_table[dgbase->dgmath_max_storage *degH + degh];
  }
  else {
    int size = (degh + 1) * (degH + 1);
    dgbase->dgmath_p_prolong_transpose_1d_inverse_table[dgbase->dgmath_max_storage *degH + degh] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_p_prolong_transpose_1d_inverse_table[dgbase->dgmath_max_storage *degH + degh];
    dgmath_build_p_prolong_transpose_1d_inverse(dgbase, op, degH, degh);
    return op;
  }
}


void dgmath_apply_hp_prolong_transpose(dgmath_jit_dbase_t* dgbase, double* in,
                                       int* degh, int dim, int degH,
                                       double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int children = (1 << dim);
  int c;
  int nodesH = dgmath_get_nodes(dim, degH);
  double* hp_prolong_transpose_1d;

  int stride = 0;
  double* tmp = P4EST_ALLOC(double, nodesH);
  linalg_fill_vec(out, 0., nodesH);

  for (c = 0; c < children; c++) {
    hp_prolong_transpose_1d =
        dgmath_fetch_hp_prolong_transpose_1d(dgbase, degH, degh[c]);
    /* int nodesh = dgmath_get_nodes(dim,degh[c]); */
    dgmath_hp_apply_nd_prolong_transpose_with_ptr(
        tmp, degH, &in[stride], degh[c], dim, c, hp_prolong_transpose_1d);
    /* printf("child = %d\n",c); */
    /* util_print_matrix(tmp, nodesH, 1, "tmp = ", 0); */
    /* util_print_matrix(&uh[stride], nodesh, 1, "uh[stride] = ", 0); */
    linalg_vec_axpy(1.0, tmp, out, nodesH);
    /* util_print_matrix(uH, nodesH, 1, "uH = ", 0); */
    stride += dgmath_get_nodes(dim, degh[c]);
  }

  P4EST_FREE(tmp);
}

void dgmath_apply_p_prolong_transpose(dgmath_jit_dbase_t* dgbase, double* in,
                                      int degh, int dim, int degH,
                                      double* out) {
  double* p_prolong_transpose_1d =
      dgmath_fetch_p_prolong_transpose_1d(dgbase, degH, degh);

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  /* prolong matrix is identity if degrees are the same */
  if (degh == degH) {
    linalg_copy_1st_to_2nd(in, out, dgmath_get_nodes(dim, degh));
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
    mpi_abort("ERROR: dgmath_apply_p_prolong_transpose\n");
  }
}

void dgmath_project_side_onto_mortar_space(dgmath_jit_dbase_t* dgbase,
                                           double* in_side, int faces_side,
                                           int* deg_side, double* out_mortar,
                                           int faces_mortar, int* deg_mortar) {
  /* most-common so first */
  /* p-prolong to mortar space */
  if (faces_side == 1 && faces_mortar == 1) {
    dgmath_apply_p_prolong(dgbase, in_side, deg_side[0], (P4EST_DIM)-1,
                           deg_mortar[0], out_mortar);
  }

  /* hp-prolong to mortar space. */
  else if (faces_side == 1 && faces_mortar == (P4EST_HALF)) {
    dgmath_apply_hp_prolong(dgbase, in_side, deg_side[0], (P4EST_DIM)-1,
                            &deg_mortar[0], out_mortar);
  }

  /* if faces_side = (P4EST_HALF), then faces_mortar == (P4EST_HALF) */
  /* here we loop through faces and p-prolong */
  else if (faces_side == (P4EST_HALF) && faces_mortar == (P4EST_HALF)) {
    int i;
    int stride_side = 0;
    int stride_mortar = 0;
    for (i = 0; i < faces_side; i++) {
      dgmath_apply_p_prolong(dgbase, &in_side[stride_side], deg_side[i],
                             (P4EST_DIM)-1, deg_mortar[i],
                             &out_mortar[stride_mortar]);

      stride_side += dgmath_get_nodes((P4EST_DIM)-1, deg_side[i]);
      stride_mortar += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar[i]);
    }
  } else {
    mpi_abort("ERROR: dgmath_project_side_onto_mortar_space");
  }
}

void dgmath_project_mass_mortar_onto_side(dgmath_jit_dbase_t* dgmath,
                                          double* in_mortar, int faces_mortar,
                                          int* deg_mortar, double* out_side,
                                          int faces_side, int* deg_side) {

  /* most-common so first */
  /* p-prolong to mortar space */
  if (faces_side == 1 && faces_mortar == 1) {
    dgmath_apply_p_prolong_transpose(dgmath, in_mortar, deg_mortar[0], (P4EST_DIM)-1,
                            deg_side[0], out_side);
  }

  /* hp-prolong to mortar space. */
  else if (faces_side == 1 && faces_mortar == (P4EST_HALF)) {
    dgmath_apply_hp_prolong_transpose(dgmath, in_mortar, &deg_mortar[0], (P4EST_DIM)-1,
                             deg_side[0], out_side);
  }

  /* if faces_side = (P4EST_HALF), then faces_mortar == (P4EST_HALF) */
  /* here we loop through faces and p-prolong */
  else if (faces_side == (P4EST_HALF) && faces_mortar == (P4EST_HALF)) {
    int i;
    int stride_side = 0;
    int stride_mortar = 0;
    for (i = 0; i < faces_side; i++) {
      dgmath_apply_p_prolong_transpose(dgmath, &in_mortar[stride_mortar], deg_mortar[i],
                              (P4EST_DIM)-1, deg_side[i],
                              &out_side[stride_side]);

      stride_side += dgmath_get_nodes((P4EST_DIM)-1, deg_side[i]);
      stride_mortar += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar[i]);
    }
  }

  else {
    mpi_abort("ERROR: dgmath_project_mass_mortar_onto_side_space");
  }
}


void dgmath_project_mortar_onto_side(dgmath_jit_dbase_t* dgbase,
                                     double* in_mortar, int faces_mortar,
                                     int* deg_mortar, double* out_side,
                                     int faces_side, int* deg_side) {
  /* most-common so first */
  /* p-prolong to mortar space */
  if (faces_side == 1 && faces_mortar == 1) {
    dgmath_apply_p_restrict(dgbase, in_mortar, deg_mortar[0], (P4EST_DIM)-1,
                            deg_side[0], out_side);
  }

  /* hp-prolong to mortar space. */
  else if (faces_side == 1 && faces_mortar == (P4EST_HALF)) {
    dgmath_apply_hp_restrict(dgbase, in_mortar, &deg_mortar[0], (P4EST_DIM)-1,
                             deg_side[0], out_side);
  }

  /* if faces_side = (P4EST_HALF), then faces_mortar == (P4EST_HALF) */
  /* here we loop through faces and p-prolong */
  else if (faces_side == (P4EST_HALF) && faces_mortar == (P4EST_HALF)) {
    int i;
    int stride_side = 0;
    int stride_mortar = 0;
    for (i = 0; i < faces_side; i++) {
      dgmath_apply_p_restrict(dgbase, &in_mortar[stride_mortar], deg_mortar[i],
                              (P4EST_DIM)-1, deg_side[i],
                              &out_side[stride_side]);

      stride_side += dgmath_get_nodes((P4EST_DIM)-1, deg_side[i]);
      stride_mortar += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar[i]);
    }
  }

  else {
    mpi_abort("ERROR: dgmath_project_side_onto_mortar_space");
  }
}

static
void dgmath_build_hp_restrict_interp_1d_aux(
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
      phi_modal[i] = dgmath_jacobi(r, 0, 0, i);
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

static void dgmath_build_hp_restrict_interp_1d(dgmath_jit_dbase_t* dgbase,
                                               double* hp_restrict_interp_1d,
                                               int degH, int degh) {
  int edge_nodes_degH = degH + 1;
  int edge_nodes_degh = degh + 1;

  double* v1d = (double*)P4EST_ALLOC(double, edge_nodes_degh* edge_nodes_degh);
  double* inv_v1d =
      (double*)P4EST_ALLOC(double, edge_nodes_degh* edge_nodes_degh);
  double* inv_v1d_trans_degh =
      (double*)P4EST_ALLOC(double, edge_nodes_degh* edge_nodes_degh);

  dgmath_build_Vij_1d(dgbase, v1d, degh);
  linalg_invert_and_copy(v1d, inv_v1d, edge_nodes_degh);
  linalg_mat_transpose(inv_v1d, inv_v1d_trans_degh, edge_nodes_degh);

  double* lobatto_nodes_degH = dgmath_fetch_GLL_nodes_1d(dgbase, degH);
  dgmath_build_hp_restrict_interp_1d_aux(degh, degH, 0, inv_v1d_trans_degh,
                                         lobatto_nodes_degH,
                                         &hp_restrict_interp_1d[0]);
  dgmath_build_hp_restrict_interp_1d_aux(
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

static double* dgmath_fetch_hp_restrict_interp_1d(dgmath_jit_dbase_t* dgbase,
                                                  int degH, int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < dgbase->dgmath_max_storage);

  if (dgbase->dgmath_hp_restrict_interp_1d_table[dgbase->dgmath_max_storage * degH + degh] != NULL) {
    return dgbase->dgmath_hp_restrict_interp_1d_table[dgbase->dgmath_max_storage * degH + degh];
  }
  else {
    int size = 2 * (degh + 1) * (degH + 1);
    dgbase->dgmath_hp_restrict_interp_1d_table[dgbase->dgmath_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_hp_restrict_interp_1d_table[dgbase->dgmath_max_storage * degH + degh];
    dgmath_build_hp_restrict_interp_1d(dgbase, op, degH, degh);
    return op;
  }
}

void dgmath_apply_hp_restrict_interp(dgmath_jit_dbase_t* dgbase, double* in,
                                     int* degh, int dim, int degH,
                                     double* out) {
  mpi_assert(dim == 1 || dim == 2 || dim == 3);
  int children = (1 << dim);
  int c;
  int nodesH = dgmath_get_nodes(dim, degH);
  double* hp_restrict_1d;

  int stride = 0;
  double* tmp = P4EST_ALLOC(double, nodesH);
  linalg_fill_vec(out, 0., nodesH);

  for (c = 0; c < children; c++) {
    hp_restrict_1d = dgmath_fetch_hp_restrict_interp_1d(dgbase, degH, degh[c]);
    /* int nodesh = dgmath_get_nodes(dim,degh[c]); */

    dgmath_hp_apply_nd_restrict_with_ptr(tmp, degH, &in[stride], degh[c], dim,
                                         c, hp_restrict_1d);
    /* printf("child = %d\n",c); */
    /* util_print_matrix(tmp, nodesH, 1, "tmp = ", 0); */
    /* util_print_matrix(&uh[stride], nodesh, 1, "uh[stride] = ", 0); */
    linalg_vec_axpy(1.0, tmp, out, nodesH);
    /* util_print_matrix(uH, nodesH, 1, "uH = ", 0); */
    stride += dgmath_get_nodes(dim, degh[c]);
  }

  P4EST_FREE(tmp);
}

static void dgmath_build_p_restrict_interp_1d_aux(
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
      phi_modal[i] = dgmath_jacobi(r, 0, 0, i);
    }
    linalg_matvec_plus_vec(1., inv_v1d_trans_degh, phi_modal, 0.0,
                           &hp_restrict_interp_matrix_1d[n * nodesh], nodesh,
                           nodesh);
  }

  P4EST_FREE(phi_modal);
}

static void dgmath_build_p_restrict_interp_1d(dgmath_jit_dbase_t* dgbase,
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

  dgmath_build_Vij_1d(dgbase, v1d, degh);
  linalg_invert_and_copy(v1d, inv_v1d, edge_nodes_degh);
  linalg_mat_transpose(inv_v1d, inv_v1d_trans_degh, edge_nodes_degh);

  double* lobatto_nodes_degH = dgmath_fetch_GLL_nodes_1d(dgbase, degH);
  dgmath_build_p_restrict_interp_1d_aux(degh, degH, inv_v1d_trans_degh,
                                        lobatto_nodes_degH,
                                        &hp_restrict_interp_1d[0]);

  P4EST_FREE(inv_v1d);
  P4EST_FREE(inv_v1d_trans_degh);
  P4EST_FREE(v1d);
}

static double* dgmath_fetch_p_restrict_interp_1d(dgmath_jit_dbase_t* dgbase,
                                                 int degH, int degh) {
  mpi_assert(degH <= degh);
  mpi_assert(degh < dgbase->dgmath_max_storage);

  if (dgbase->dgmath_p_restrict_interp_1d_table[dgbase->dgmath_max_storage * degH + degh] != NULL) {
    return dgbase->dgmath_p_restrict_interp_1d_table[dgbase->dgmath_max_storage * degH + degh];
  }
  else {
    int size = (degh + 1) * (degH + 1);
    dgbase->dgmath_p_restrict_interp_1d_table[dgbase->dgmath_max_storage * degH + degh] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_p_restrict_interp_1d_table[dgbase->dgmath_max_storage * degH + degh];
    dgmath_build_p_restrict_interp_1d(dgbase, op, degH, degh);
    return op;
  }
}

void dgmath_apply_p_restrict_interp(dgmath_jit_dbase_t* dgbase, double* in,
                                    int degh, int dim, int degH, double* out) {
  double* p_restrict_1d = dgmath_fetch_p_restrict_interp_1d(dgbase, degH, degh);

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  /* restrict matrix should be the identity if degrees are the same */
  if (degh == degH) {
    linalg_copy_1st_to_2nd(in, out, dgmath_get_nodes(dim, degh));
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
    mpi_abort("ERROR: dgmath_restrict_nd_U\n");
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


/* only use for cross products or obtaining normal direction */
void dgmath_get_face_info(int f, dgmath_face_info_t* face_info) {

  if (f == 0) {
    face_info->a = 2;
    face_info->b = 1;
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
    face_info->a = 2;
    face_info->b = 0;
    face_info->sgn = 1.;
    face_info->c = 1;
  } else if (f == 4) {
    face_info->a = 1;
    face_info->b = 0;
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
void dgmath_compute_M_vec_wo_aliasing(dgmath_jit_dbase_t* dgmath_jit_dbase,
                                      double* xyz_rst[(P4EST_DIM)][P4EST_DIM],
                                      double* vec, int degH, double* MJvec,
                                     int degh) {
#if (P4EST_DIM)==3
  int volume_nodes_degH = dgmath_get_nodes((P4EST_DIM), degH);
  int volume_nodes_degh = dgmath_get_nodes((P4EST_DIM), degh);

  double* xyz_rst_prolonged[(P4EST_DIM)][P4EST_DIM];
  double* vec_prolonged = P4EST_ALLOC(double, volume_nodes_degh);
  double* Jvec_prolonged = P4EST_ALLOC(double, volume_nodes_degh);
  double* Jvec = P4EST_ALLOC(double, volume_nodes_degH);

  int i, j;
  for (i = 0; i < (P4EST_DIM); i++) {
    for (j = 0; j < (P4EST_DIM); j++) {
      xyz_rst_prolonged[i][j] = P4EST_ALLOC(double, volume_nodes_degh);
      dgmath_apply_p_prolong(dgmath_jit_dbase, &(xyz_rst[i][j][0]), degH,
                             (P4EST_DIM), degh, &(xyz_rst_prolonged[i][j][0]));
    }
  }

  dgmath_apply_p_prolong(dgmath_jit_dbase, vec, degH, (P4EST_DIM), degh,
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

  dgmath_apply_p_restrict(dgmath_jit_dbase, Jvec_prolonged, degh, (P4EST_DIM),
                          degH, Jvec);

  dgmath_apply_Mij(dgmath_jit_dbase, Jvec, (P4EST_DIM), degH, MJvec);

  for (i = 0; i < (P4EST_DIM); i++) {
    for (j = 0; j < (P4EST_DIM); j++) {
      P4EST_FREE(xyz_rst_prolonged[i][j]);
    }
  }
  
  P4EST_FREE(vec_prolonged);
  P4EST_FREE(Jvec_prolonged);
  P4EST_FREE(Jvec);
#else
  mpi_abort("dgmath_compute_M_vec_wo_aliasing currently only supports DIM=3");
#endif
}


void dgmath_compute_invM_vec_wo_aliasing(dgmath_jit_dbase_t* dgmath_jit_dbase,
                                         double* xyz_rst[(P4EST_DIM)][P4EST_DIM],
                                         double* vec, int degH, double* invMvec,
                                         int degh) {
#if (P4EST_DIM)==3 
  /* int volume_nodes_degH = dgmath_get_nodes((P4EST_DIM), degH); */
  int volume_nodes_degh = dgmath_get_nodes((P4EST_DIM), degh);

  double* xyz_rst_prolonged[(P4EST_DIM)][P4EST_DIM];
  double* invMvec_prolonged = P4EST_ALLOC(double, volume_nodes_degh);
  double* invMvec_divJ_prolonged = P4EST_ALLOC(double, volume_nodes_degh);

  int i, j;
  for (i = 0; i < (P4EST_DIM); i++) {
    for (j = 0; j < (P4EST_DIM); j++) {
      xyz_rst_prolonged[i][j] = P4EST_ALLOC(double, volume_nodes_degh);
      dgmath_apply_p_prolong(dgmath_jit_dbase, &(xyz_rst[i][j][0]), degH,
                             (P4EST_DIM), degh, &(xyz_rst_prolonged[i][j][0]));
    }
  }

  dgmath_apply_invMij(dgmath_jit_dbase, vec, (P4EST_DIM), degH, invMvec);
  dgmath_apply_p_prolong(dgmath_jit_dbase, invMvec, degH, (P4EST_DIM), degh,
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

  dgmath_apply_p_restrict(dgmath_jit_dbase, invMvec_divJ_prolonged, degh, (P4EST_DIM),
                          degH, invMvec);

  for (i = 0; i < (P4EST_DIM); i++) {
    for (j = 0; j < (P4EST_DIM); j++) {
      P4EST_FREE(xyz_rst_prolonged[i][j]);
    }
  }
   
  P4EST_FREE(invMvec_divJ_prolonged);
  P4EST_FREE(invMvec_prolonged);
#else
  mpi_abort(" dgmath_compute_invM_vec_wo_aliasing currently only supports DIM=3 ");
#endif
}

static void dgmath_build_FLIP_1d(double* FLIP_1d, int deg) {
  memset(FLIP_1d, 0., sizeof(double) * (deg + 1) * (deg + 1));
  int offset;
  int i;
  for (i = 0; i < deg + 1; i++) {
    offset = (i+1)*deg;
    FLIP_1d[offset] = 1.;
  }
}

static double* dgmath_fetch_FLIP_1d(dgmath_jit_dbase_t* dgbase, int deg) {
  mpi_assert(deg < dgbase->dgmath_max_storage);
  if (dgbase->dgmath_FLIP_1d_table[deg] != NULL) {
    return dgbase->dgmath_FLIP_1d_table[deg];
  }
  else {
    int size = (deg + 1) * (deg + 1);
    dgbase->dgmath_FLIP_1d_table[deg] = P4EST_ALLOC(double, size);
    double* op = dgbase->dgmath_FLIP_1d_table[deg];
    dgmath_build_FLIP_1d(op, deg);
    return op;
  }
}

void dgmath_apply_FLIP(dgmath_jit_dbase_t* dgbase, double* in, int dim, int deg,
                       int dir, double* out) {
  mpi_assert(dim == 1 || dim == 2);
  double* flip_1d = dgmath_fetch_FLIP_1d(dgbase, deg);
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

int dgmath_reorient_face_order
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
    int perm = dgmath_p8est_code_to_perm[dgmath_p8est_FToF_code[f_m][f_p]][o];
    return dgmath_p8est_perm_to_order[perm][i];
  }
  else {
    mpi_abort("FACE DIM == 1, 2\n");
    return -1;
  }
}



void dgmath_reorient_face_data
(
 dgmath_jit_dbase_t* dgbase,
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
  
  int nodes = dgmath_get_nodes(face_dim, deg);
  
  if (face_dim == 1){
    if (o == 1)
      dgmath_apply_FLIP
        (
         dgbase,
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
    dgmath_apply_FLIP
      (
       dgbase,
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
    dgmath_apply_FLIP
      (
       dgbase,
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


dgmath_rst_t
dgmath_get_rst_points
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 int deg,
 int dim,
 quadrature_type_t type
)
{
  dgmath_rst_t rst;
  double* rsttmp [3] = {NULL};
  
  if (type == GAUSS){
    for (int d = 0; d < dim; d++){
      rsttmp[d] = dgmath_fetch_Gauss_xyz_nd
                  (
                   dgmath_jit_dbase,
                   dim,
                   deg,
                   d
                  );
    }
  }
  else if (type == LOBATTO) {
    for (int d = 0; d < dim; d++){
      rsttmp[d] = dgmath_fetch_xyz_nd
                  (
                   dgmath_jit_dbase,
                   dim,
                   deg,
                   d
                  );
    }
  }
  else {
    mpi_abort("type == Lobatto or Gauss");
  }

  rst.r = rsttmp[0];
  rst.s = rsttmp[1];
  rst.t = rsttmp[2];
  
  return rst;
}


int dgmath_corner_to_node
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
