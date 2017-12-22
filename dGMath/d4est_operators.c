#include <pXest.h>
#include <d4est_operators.h>
#include <assert.h>
#include <d4est_linalg.h>
#include <d4est_kron.h>
#include <d4est_util.h>
#include <GL_and_GLL_nodes_and_weights.h>

d4est_operators_t* d4est_ops_init(int max_degree) {
  d4est_operators_t* d4est_ops = P4EST_ALLOC(d4est_operators_t, 1);
  d4est_ops->max_degree = max_degree;
  
  d4est_ops->mij_1d_table = P4EST_ALLOC(double*, max_degree);
  d4est_ops->invmij_1d_table = P4EST_ALLOC(double*, max_degree);
  d4est_ops->invvij_1d_table = P4EST_ALLOC(double*, max_degree);
  d4est_ops->lobatto_nodes_1d_table = P4EST_ALLOC(double*, max_degree);
  d4est_ops->lobatto_weights_1d_table = P4EST_ALLOC(double*, max_degree);
  d4est_ops->gauss_nodes_1d_table = P4EST_ALLOC(double*, max_degree);
  d4est_ops->gauss_weights_1d_table = P4EST_ALLOC(double*, max_degree);
  d4est_ops->dij_1d_table = P4EST_ALLOC(double*, max_degree);
  d4est_ops->dij_trans_1d_table = P4EST_ALLOC(double*, max_degree);
  d4est_ops->lift_1d_table = P4EST_ALLOC(double*, max_degree);
  d4est_ops->vtk_rst_2d_table = P4EST_ALLOC(double*, max_degree); 
  d4est_ops->vtk_rst_3d_table = P4EST_ALLOC(double*, max_degree); 
  d4est_ops->lobatto_rst_2d_table = P4EST_ALLOC(double*, max_degree); 
  d4est_ops->lobatto_rst_3d_table = P4EST_ALLOC(double*, max_degree); 
  d4est_ops->gauss_rst_2d_table = P4EST_ALLOC(double*, max_degree); 
  d4est_ops->gauss_rst_3d_table = P4EST_ALLOC(double*, max_degree); 
  d4est_ops->flip_1d_table = P4EST_ALLOC(double*, max_degree);
  d4est_ops->vtk_interp_1d_table = P4EST_ALLOC(double*, max_degree);
  
  d4est_ops->hp_prolong_1d_table = P4EST_ALLOC(double**, max_degree);
  d4est_ops->hp_restrict_1d_table = P4EST_ALLOC(double**, max_degree);
  d4est_ops->hp_restrict_interp_1d_table = P4EST_ALLOC(double**, max_degree);
  d4est_ops->p_restrict_1d_table = P4EST_ALLOC(double**, max_degree);
  d4est_ops->p_restrict_interp_1d_table = P4EST_ALLOC(double**, max_degree);
  d4est_ops->p_prolong_1d_table = P4EST_ALLOC(double**, max_degree);
  d4est_ops->hp_prolong_transpose_1d_table = P4EST_ALLOC(double**, max_degree);
  d4est_ops->p_prolong_transpose_1d_table = P4EST_ALLOC(double**, max_degree);
  d4est_ops->lobatto_to_gauss_interp_1d_table = P4EST_ALLOC(double**, max_degree);
  d4est_ops->lobatto_to_gauss_interp_trans_1d_table = P4EST_ALLOC(double**, max_degree);

  for (int i = 0; i < d4est_ops->max_degree; i++){
    d4est_ops->hp_prolong_1d_table[i] = P4EST_ALLOC(double*, max_degree);
    d4est_ops->hp_prolong_transpose_1d_table[i] = P4EST_ALLOC(double*, max_degree);
    d4est_ops->hp_restrict_1d_table[i] = P4EST_ALLOC(double*, max_degree);
    d4est_ops->hp_restrict_interp_1d_table[i] = P4EST_ALLOC(double*, max_degree);
    d4est_ops->p_restrict_1d_table[i] = P4EST_ALLOC(double*, max_degree);
    d4est_ops->p_restrict_interp_1d_table[i] = P4EST_ALLOC(double*, max_degree);
    d4est_ops->p_prolong_1d_table[i] = P4EST_ALLOC(double*, max_degree);
    d4est_ops->p_prolong_transpose_1d_table[i] = P4EST_ALLOC(double*, max_degree);
    d4est_ops->lobatto_to_gauss_interp_trans_1d_table[i] = P4EST_ALLOC(double*, max_degree);
    d4est_ops->lobatto_to_gauss_interp_1d_table[i] = P4EST_ALLOC(double*, max_degree);
    d4est_ops->vtk_interp_1d_table[i] = NULL;
    d4est_ops->mij_1d_table[i] = NULL;
    d4est_ops->lobatto_nodes_1d_table[i] = NULL;
    d4est_ops->lobatto_weights_1d_table[i] = NULL;
    d4est_ops->gauss_nodes_1d_table[i] = NULL;
    d4est_ops->gauss_weights_1d_table[i] = NULL;
    d4est_ops->invmij_1d_table[i] = NULL;
    d4est_ops->vtk_rst_3d_table[i] = NULL;
    d4est_ops->lobatto_rst_3d_table[i] = NULL;    
    d4est_ops->lobatto_rst_2d_table[i] = NULL;    
    d4est_ops->gauss_rst_2d_table[i] = NULL;
    d4est_ops->gauss_rst_3d_table[i] = NULL;
    d4est_ops->vtk_rst_2d_table[i] = NULL;
    d4est_ops->invvij_1d_table[i] = NULL;
    d4est_ops->dij_1d_table[i] = NULL;
    d4est_ops->dij_trans_1d_table[i] = NULL;
    d4est_ops->lift_1d_table[i] = NULL;
    d4est_ops->flip_1d_table[i] = NULL;

    for (int j = 0; j < d4est_ops->max_degree; j++){
      d4est_ops->hp_prolong_1d_table[i][j] = NULL;
      d4est_ops->hp_prolong_transpose_1d_table[i][j] = NULL;
      d4est_ops->p_prolong_transpose_1d_table[i][j] = NULL;
      d4est_ops->hp_restrict_1d_table[i][j] = NULL;
      d4est_ops->hp_restrict_interp_1d_table[i][j] = NULL;
      d4est_ops->p_restrict_1d_table[i][j] = NULL;
      d4est_ops->p_restrict_interp_1d_table[i][j] = NULL;
      d4est_ops->p_prolong_1d_table[i][j] = NULL;
      d4est_ops->lobatto_to_gauss_interp_trans_1d_table[i][j] = NULL;
      d4est_ops->lobatto_to_gauss_interp_1d_table[i][j] = NULL;
    }
  }
  
  return d4est_ops;
}

void d4est_ops_destroy(d4est_operators_t* d4est_ops) {

  for (int i = 0; i < d4est_ops->max_degree; i++) {
    P4EST_FREE(d4est_ops->mij_1d_table[i]);
    P4EST_FREE(d4est_ops->vtk_rst_3d_table[i]);
    P4EST_FREE(d4est_ops->vtk_rst_2d_table[i]);
    P4EST_FREE(d4est_ops->lobatto_nodes_1d_table[i]);
    P4EST_FREE(d4est_ops->lobatto_weights_1d_table[i]);
    P4EST_FREE(d4est_ops->gauss_nodes_1d_table[i]);
    P4EST_FREE(d4est_ops->gauss_weights_1d_table[i]);
    P4EST_FREE(d4est_ops->invmij_1d_table[i]);
    P4EST_FREE(d4est_ops->invvij_1d_table[i]);
    P4EST_FREE(d4est_ops->dij_1d_table[i]);
    P4EST_FREE(d4est_ops->dij_trans_1d_table[i]);
    P4EST_FREE(d4est_ops->lift_1d_table[i]);
    P4EST_FREE(d4est_ops->flip_1d_table[i]);
    P4EST_FREE(d4est_ops->vtk_interp_1d_table[i]);
    P4EST_FREE(d4est_ops->lobatto_rst_3d_table[i]);
    P4EST_FREE(d4est_ops->lobatto_rst_2d_table[i]);
    P4EST_FREE(d4est_ops->gauss_rst_2d_table[i]);
    P4EST_FREE(d4est_ops->gauss_rst_3d_table[i]);
  }
  P4EST_FREE(d4est_ops->mij_1d_table);
  P4EST_FREE(d4est_ops->vtk_rst_3d_table);
  P4EST_FREE(d4est_ops->vtk_rst_2d_table);
  P4EST_FREE(d4est_ops->vtk_interp_1d_table);
  P4EST_FREE(d4est_ops->lobatto_nodes_1d_table);
  P4EST_FREE(d4est_ops->lobatto_weights_1d_table);
  P4EST_FREE(d4est_ops->gauss_nodes_1d_table);
  P4EST_FREE(d4est_ops->gauss_weights_1d_table);
  P4EST_FREE(d4est_ops->invmij_1d_table);
  P4EST_FREE(d4est_ops->invvij_1d_table);
  P4EST_FREE(d4est_ops->dij_1d_table);
  P4EST_FREE(d4est_ops->dij_trans_1d_table);
  P4EST_FREE(d4est_ops->lift_1d_table);
  P4EST_FREE(d4est_ops->flip_1d_table);
  P4EST_FREE(d4est_ops->lobatto_rst_3d_table);
  P4EST_FREE(d4est_ops->lobatto_rst_2d_table);
  P4EST_FREE(d4est_ops->gauss_rst_2d_table);
  P4EST_FREE(d4est_ops->gauss_rst_3d_table);

  for (int i = 0; i < d4est_ops->max_degree; i++){
    for(int j = 0; j < d4est_ops->max_degree; j++){
      P4EST_FREE(d4est_ops->hp_prolong_1d_table[i][j]);
      P4EST_FREE(d4est_ops->hp_prolong_transpose_1d_table[i][j]);
      P4EST_FREE(d4est_ops->p_prolong_transpose_1d_table[i][j]);
      P4EST_FREE(d4est_ops->hp_restrict_1d_table[i][j]);
      P4EST_FREE(d4est_ops->hp_restrict_interp_1d_table[i][j]);
      P4EST_FREE(d4est_ops->p_prolong_1d_table[i][j]);
      P4EST_FREE(d4est_ops->p_restrict_1d_table[i][j]);
      P4EST_FREE(d4est_ops->p_restrict_interp_1d_table[i][j]);
      P4EST_FREE(d4est_ops->lobatto_to_gauss_interp_trans_1d_table[i][j]);
      P4EST_FREE(d4est_ops->lobatto_to_gauss_interp_1d_table[i][j]);
    }
    P4EST_FREE(d4est_ops->hp_prolong_1d_table[i]);
    P4EST_FREE(d4est_ops->hp_prolong_transpose_1d_table[i]);
    P4EST_FREE(d4est_ops->p_prolong_transpose_1d_table[i]);
    P4EST_FREE(d4est_ops->hp_restrict_1d_table[i]);
    P4EST_FREE(d4est_ops->hp_restrict_interp_1d_table[i]);
    P4EST_FREE(d4est_ops->p_prolong_1d_table[i]);
    P4EST_FREE(d4est_ops->p_restrict_1d_table[i]);
    P4EST_FREE(d4est_ops->p_restrict_interp_1d_table[i]);
    P4EST_FREE(d4est_ops->lobatto_to_gauss_interp_trans_1d_table[i]);
    P4EST_FREE(d4est_ops->lobatto_to_gauss_interp_1d_table[i]);
  }
  P4EST_FREE(d4est_ops->hp_prolong_1d_table);
  P4EST_FREE(d4est_ops->hp_prolong_transpose_1d_table);
  P4EST_FREE(d4est_ops->p_prolong_transpose_1d_table);
  P4EST_FREE(d4est_ops->hp_restrict_1d_table);
  P4EST_FREE(d4est_ops->hp_restrict_interp_1d_table);
  P4EST_FREE(d4est_ops->p_restrict_1d_table);
  P4EST_FREE(d4est_ops->p_restrict_interp_1d_table);
  P4EST_FREE(d4est_ops->p_prolong_1d_table);
  P4EST_FREE(d4est_ops->lobatto_to_gauss_interp_trans_1d_table);
  P4EST_FREE(d4est_ops->lobatto_to_gauss_interp_1d_table);
     
  P4EST_FREE(d4est_ops);
}

double*
d4est_operators_1index_fetch
(
 d4est_operators_t* d4est_ops,
 double** table,
 int deg,
 int size,
 void (*build_fcn)(d4est_operators_t*, double*, int)
)
{
  D4EST_ASSERT(deg < d4est_ops->max_degree && deg > 0);
  if (table[deg-1] != NULL) {
    return table[deg-1];
  }
  else {
    table[deg-1] = P4EST_ALLOC(double, size);
    double* op = table[deg-1];
    build_fcn(d4est_ops, op, deg);
    return op;
  }
}

double* d4est_operators_2index_fetch
(
 d4est_operators_t* d4est_ops,
 double*** table,
 int deg1,
 int deg2,
 int size,
 void (*build_fcn)(d4est_operators_t*, double*, int, int)
)
{
  D4EST_ASSERT(deg1 < d4est_ops->max_degree &&
             deg2 < d4est_ops->max_degree &&
             deg1 > 0 &&
             deg2 > 0);
  
  if (table[deg1-1][deg2-1] != NULL) {
    return table[deg1-1][deg2-1];
  }
  else {
    table[deg1-1][deg2-1] = P4EST_ALLOC(double, size);
    double* op = table[deg1-1][deg2-1];
    build_fcn(d4est_ops, op, deg1, deg2);
    return op;
  }
}


double* d4est_operators_1index_2d_3d_fetch
(
 d4est_operators_t* d4est_ops,
 int deg,
 int dim,
 int size,
 double** table_2d,
 double** table_3d,
 void (*build_fcn)(d4est_operators_t*, double*, int, int)
)
{
  D4EST_ASSERT(deg < d4est_ops->max_degree &&
             dim <= 3 &&
             deg > 0 &&
             dim > 0
            );
  double* table_ptr = (dim==3) ? table_3d[deg] : table_2d[deg];
  
  if (table_ptr != NULL) {
    return table_ptr;
  }
  else {
    if (dim == 3){
      table_3d[deg] = P4EST_ALLOC(double, size);
    }
    else if (dim == 2){
      table_2d[deg] = P4EST_ALLOC(double, size);
    }
    else{
      D4EST_ABORT("[D4EST_ERROR]: Not a support dimension");
    }
    
    double* op = (dim==3) ? table_3d[deg] : table_2d[deg];
    build_fcn(d4est_ops, op, dim, deg);
    return op;
  }
}


void d4est_operators_build_Vij_1d(d4est_operators_t* d4est_ops, double* Vij_1d,
                                int deg) {
  int i, j, rows, cols;
  const double* lobatto_nodes = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, deg);
  rows = cols = deg + 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      Vij_1d[i * cols + j] = d4est_lgl_jacobi(lobatto_nodes[i], 0., 0., j);
}

static void d4est_operators_build_invvij_1d(d4est_operators_t* d4est_ops,
                                   double* invvij_1d, int deg)
{
  d4est_operators_build_Vij_1d(d4est_ops, invvij_1d, deg);
  d4est_linalg_invert(invvij_1d, deg + 1);
}

double* d4est_operators_fetch_invvij_1d(d4est_operators_t* d4est_ops, int deg) {
  int size = (deg + 1) * (deg + 1);
  return d4est_operators_1index_fetch
    (
     d4est_ops,
     d4est_ops->invvij_1d_table,
     deg,
     size,
     d4est_operators_build_invvij_1d
    );     
}

void d4est_operators_hp_apply_nd_prolong_with_ptr(double* Uh, int degh, double* UH,
                                         int degH, int dim, int c,
                                         double* hp_prolong_matrix_1d) {
  D4EST_ASSERT((degH <= degh) && (c < (1 << dim)));
  int nodesh = degh + 1;
  int nodesH = degH + 1;

  if (dim == 1) {
    d4est_linalg_matvec_plus_vec(1.0, &hp_prolong_matrix_1d[c * nodesh * nodesH], UH,
                           0., Uh, nodesh, nodesH);
  } else if (dim == 2) {
    int cx = d4est_reference_is_child_left_or_right(c, 0);
    int cy = d4est_reference_is_child_left_or_right(c, 1);

       d4est_kron_A1A2x_nonsqr(Uh, &hp_prolong_matrix_1d[cy * nodesh * nodesH],
                           &hp_prolong_matrix_1d[cx * nodesh * nodesH], UH,
                           nodesh, nodesH, nodesh, nodesH);
    
  } else if (dim == 3) {
    int cx = d4est_reference_is_child_left_or_right(c, 0);
    int cy = d4est_reference_is_child_left_or_right(c, 1);
    int cz = d4est_reference_is_child_left_or_right(c, 2);


    d4est_kron_A1A2A3x_nonsqr(Uh,
                               &hp_prolong_matrix_1d[cz * nodesh * nodesH],
                               &hp_prolong_matrix_1d[cy * nodesh * nodesH],
                               &hp_prolong_matrix_1d[cx * nodesh * nodesH], UH,
                               nodesh, nodesH, nodesh, nodesH, nodesh, nodesH);

  } else {
    D4EST_ABORT("[D4EST_ERROR]: d4est_operators_prolong_nd_U\n");
  }
}

void d4est_operators_build_lobatto_to_gauss_interp_1d
(
 d4est_operators_t* d4est_ops,
 double* ref_lobatto_to_gauss_interp_1d,
 int lobatto_degree,
 int gauss_degree
)
{
  double* gauss_nodes = d4est_operators_fetch_gauss_nodes_1d(d4est_ops, gauss_degree);
  double* ref_gaussVij = P4EST_ALLOC(double, (gauss_degree + 1)*(lobatto_degree+1));
  
  int rows = gauss_degree + 1;
  int cols = lobatto_degree + 1;
  
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      ref_gaussVij[i * cols + j] = d4est_lgl_jacobi(gauss_nodes[i], 0., 0., j);

  double* invvij = d4est_operators_fetch_invvij_1d(d4est_ops, lobatto_degree);

  /* DEBUG_PRINT_ARR_DBL(invvij, (lobatto_degree+1)*(lobatto_degree+1)); */
  /* DEBUG_PRINT_ARR_DBL(ref_gaussVij, (lobatto_degree+1)*(gauss_degree+1)); */
  d4est_linalg_mat_multiply(ref_gaussVij, invvij, ref_lobatto_to_gauss_interp_1d, gauss_degree + 1, lobatto_degree + 1, lobatto_degree + 1);

  /* DEBUG_PRINT_ARR_DBL(ref_lobatto_to_gauss_interp_1d, (gauss_degree + 1)*(lobatto_degree + 1)); */
  
  P4EST_FREE(ref_gaussVij);
}

void d4est_operators_build_custom_lobatto_interp_1d
(
 d4est_operators_t* d4est_ops,
 double* custom_lobatto_interp_1d,
 int lobatto_degree,
 int custom_degree,
 double* custom_points /* size = custom_degree + 1 */
)
{
  double* customVij = P4EST_ALLOC(double, (custom_degree + 1)*(lobatto_degree+1));
  
  int rows = custom_degree + 1;
  int cols = lobatto_degree + 1;
  
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      customVij[i * cols + j] = d4est_lgl_jacobi(custom_points[i], 0., 0., j);

  double* invvij = d4est_operators_fetch_invvij_1d(d4est_ops, lobatto_degree);
  d4est_linalg_mat_multiply(customVij, invvij, custom_lobatto_interp_1d, custom_degree + 1, lobatto_degree + 1, lobatto_degree + 1);

  P4EST_FREE(customVij);
}



double*
d4est_operators_fetch_lobatto_to_gauss_interp_1d(d4est_operators_t* d4est_ops, int deg_lobatto,
                                         int deg_gauss) {

  int size = (deg_gauss + 1) * (deg_lobatto + 1);
  return d4est_operators_2index_fetch(d4est_ops,
                                      d4est_ops->lobatto_to_gauss_interp_1d_table,
                                      deg_lobatto,
                                      deg_gauss,
                                      size,
                                      d4est_operators_build_lobatto_to_gauss_interp_1d
                                     );
 
}

void d4est_operators_build_lobatto_to_gauss_interp_trans_1d
(
 d4est_operators_t* d4est_ops,
 double* ref_lobatto_to_gauss_interp_trans_1d,
 int lobatto_degree,
 int gauss_degree
)
{
  double* ref_lobatto_to_gauss_interp_1d = d4est_operators_fetch_lobatto_to_gauss_interp_1d(d4est_ops, lobatto_degree, gauss_degree);
  d4est_linalg_mat_transpose_nonsqr(ref_lobatto_to_gauss_interp_1d, ref_lobatto_to_gauss_interp_trans_1d, gauss_degree + 1, lobatto_degree + 1);
}

double* d4est_operators_fetch_lobatto_to_gauss_interp_trans_1d
(
 d4est_operators_t* d4est_ops,
 int deg_lobatto,
 int deg_gauss
)
{
  int size = (deg_gauss + 1) * (deg_lobatto + 1);
  return d4est_operators_2index_fetch(d4est_ops,
                                      d4est_ops->lobatto_to_gauss_interp_trans_1d_table,
                                      deg_lobatto,
                                      deg_gauss,
                                      size,
                                      d4est_operators_build_lobatto_to_gauss_interp_trans_1d
                                     );  
}

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
    volume_nodes_h += d4est_lgl_get_nodes(dim, degh[i]);
  }
  int volume_nodes_H = d4est_lgl_get_nodes(dim, degH);

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
    d4est_linalg_set_column(prolong_mat, Mu, i, volume_nodes_h, volume_nodes_H);
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
    volume_nodes_h[i] = d4est_lgl_get_nodes(dim, degh[i]);
    total_volume_nodes_h += volume_nodes_h[i];
    max_volume_nodes_h = (volume_nodes_h[i] > max_volume_nodes_h) ? volume_nodes_h[i] : max_volume_nodes_h;
  }

  int volume_nodes_H = d4est_lgl_get_nodes(dim, degH);

  d4est_linalg_fill_vec(PT_mat_P, 0., volume_nodes_H*volume_nodes_H);
  double* P = P4EST_ALLOC(double, total_volume_nodes_h*volume_nodes_H);
  double* PT = P4EST_ALLOC(double, total_volume_nodes_h*volume_nodes_H);
  double* mat_P_i = P4EST_ALLOC(double, max_volume_nodes_h*volume_nodes_H);
  double* PT_mat_P_i = P4EST_ALLOC(double, max_volume_nodes_h*volume_nodes_H);

  d4est_operators_compute_prolong_matrix(d4est_ops, degH, dim, degh, children, P);
  d4est_linalg_mat_transpose_nonsqr(P, PT, total_volume_nodes_h, volume_nodes_H);

  int stride_mat = 0;
  int stride_P = 0;
  for (int i = 0; i < children; i++){
    d4est_linalg_mat_multiply(&mat[stride_mat],
                        &P[stride_P],
                        mat_P_i,
                        volume_nodes_h[i],
                        volume_nodes_h[i],
                        volume_nodes_H);
    
    d4est_linalg_mat_multiply(&PT[stride_P], mat_P_i, PT_mat_P_i, volume_nodes_H, volume_nodes_h[i], volume_nodes_H);
    d4est_linalg_vec_axpy(1., PT_mat_P_i, PT_mat_P, volume_nodes_H*volume_nodes_H);
    
    stride_mat += volume_nodes_h[i]*volume_nodes_h[i];
    stride_P += volume_nodes_h[i]*volume_nodes_H;
  }

  P4EST_FREE(PT_mat_P_i);
  P4EST_FREE(mat_P_i);
  P4EST_FREE(PT);
  P4EST_FREE(P);
}

void d4est_operators_hp_apply_nd_restrict_with_ptr(double* uH, int degH, double* uh,
                                          int degh, int dim, int c,
                                          double* hp_restrict_matrix_1d) {

  D4EST_ASSERT((degH <= degh) && (c < (1 << dim)));

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  if (dim == 1)
    d4est_linalg_matvec_plus_vec(1.0, &hp_restrict_matrix_1d[c * nodesh * nodesH], uh,
                           0., uH, nodesH, nodesh);
  else if (dim == 2) {
    int cx = d4est_reference_is_child_left_or_right(c, 0);
    int cy = d4est_reference_is_child_left_or_right(c, 1);



    d4est_kron_A1A2x_nonsqr(uH, &hp_restrict_matrix_1d[cy * nodesh * nodesH],
                             &hp_restrict_matrix_1d[cx * nodesh * nodesH], uh,
                             nodesH, nodesh, nodesH, nodesh);

  } else if (dim == 3) {
    int cx = d4est_reference_is_child_left_or_right(c, 0);
    int cy = d4est_reference_is_child_left_or_right(c, 1);
    int cz = d4est_reference_is_child_left_or_right(c, 2);


    d4est_kron_A1A2A3x_nonsqr(uH, &hp_restrict_matrix_1d[cz * nodesh * nodesH],
                               &hp_restrict_matrix_1d[cy * nodesh * nodesH],
                               &hp_restrict_matrix_1d[cx * nodesh * nodesH], uh,
                               nodesH, nodesh, nodesH, nodesh, nodesH, nodesh);

    
  } else {
    D4EST_ABORT("[D4EST_ERROR]: d4est_operators_restrict_nd_U\n");
  }
}

static void d4est_operators_build_dVij_1d(d4est_operators_t* d4est_ops, double* dVij_1d,
                                 int deg) {
  int i, j, rows, cols;
  double* lobatto_nodes = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, deg);
  rows = cols = deg + 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      dVij_1d[i * cols + j] = d4est_lgl_gradjacobi(lobatto_nodes[i], 0., 0., j);
}

void d4est_operators_build_mij_1d(d4est_operators_t* d4est_ops, double* mij_1d,
                                int deg) {
  int edge_nodes = deg + 1;
  double* v1d = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);
  double* v1d_trans = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);
  d4est_operators_build_Vij_1d(d4est_ops, v1d, deg);
  d4est_linalg_mat_transpose(v1d, v1d_trans, edge_nodes);
  d4est_linalg_mat_multiply(v1d, v1d_trans, mij_1d, edge_nodes, edge_nodes,
                      edge_nodes);
  d4est_linalg_invert(mij_1d, edge_nodes);
  P4EST_FREE(v1d_trans);
  P4EST_FREE(v1d);
}


static void
d4est_operators_build_lobatto_nodes_1d(d4est_operators_t* d4est_ops, double* lobatto_nodes, int deg)
{
  double* lobatto_weights = P4EST_ALLOC(double, deg+1);
  d4est_operators_lobatto_nodes_and_weights(deg+1, lobatto_nodes, lobatto_weights);
  P4EST_FREE(lobatto_weights);
}

static void
d4est_operators_build_lobatto_weights_1d(d4est_operators_t* d4est_ops, double* lobatto_weights, int deg)
{
  double* lobatto_nodes = P4EST_ALLOC(double, deg+1);
  d4est_operators_lobatto_nodes_and_weights(deg+1, lobatto_nodes, lobatto_weights);
  P4EST_FREE(lobatto_nodes);
}


double* d4est_operators_fetch_lobatto_nodes_1d(d4est_operators_t* d4est_ops, int deg) {
  int size = (deg + 1);
  return d4est_operators_1index_fetch
    (
     d4est_ops,
     d4est_ops->lobatto_nodes_1d_table,
     deg,
     size,
     d4est_operators_build_lobatto_nodes_1d
    );      
}

double* d4est_operators_fetch_lobatto_weights_1d(d4est_operators_t* d4est_ops, int deg) {
  int size = (deg + 1);
  return d4est_operators_1index_fetch
    (
     d4est_ops,
     d4est_ops->lobatto_weights_1d_table,
     deg,
     size,
     d4est_operators_build_lobatto_weights_1d
    );      
}

static void
d4est_operators_build_gauss_nodes_1d(d4est_operators_t* d4est_ops, double* gauss_nodes, int deg)
{
  double* gauss_weights = P4EST_ALLOC(double, deg+1);
  d4est_operators_gauss_nodes_and_weights(deg+1, gauss_nodes, gauss_weights);
  P4EST_FREE(gauss_weights);
}

static void
d4est_operators_build_gauss_weights_1d(d4est_operators_t* d4est_ops, double* gauss_weights, int deg)
{
  double* gauss_nodes = P4EST_ALLOC(double, deg+1);
  d4est_operators_gauss_nodes_and_weights(deg+1, gauss_nodes, gauss_weights);
  P4EST_FREE(gauss_nodes);
}


double* d4est_operators_fetch_gauss_nodes_1d(d4est_operators_t* d4est_ops, int deg) {
   int size = (deg + 1);
  return d4est_operators_1index_fetch
    (
     d4est_ops,
     d4est_ops->gauss_nodes_1d_table,
     deg,
     size,
     d4est_operators_build_gauss_nodes_1d
    );      
}

double* d4est_operators_fetch_gauss_weights_1d(d4est_operators_t* d4est_ops, int deg) {
  int size = (deg + 1);
  return d4est_operators_1index_fetch
    (
     d4est_ops,
     d4est_ops->gauss_weights_1d_table,
     deg,
     size,
     d4est_operators_build_gauss_weights_1d
    );      
}

double*
d4est_operators_fetch_mij_1d
(
 d4est_operators_t* d4est_ops,
 int deg
)
{
  int size = (deg + 1) * (deg + 1);
  return d4est_operators_1index_fetch
    (
     d4est_ops,
     d4est_ops->mij_1d_table,
     deg,
     size,
     d4est_operators_build_mij_1d
    );    
}

static void d4est_operators_build_invmij_1d(d4est_operators_t* d4est_ops,
                                   double* invmij_1d, int deg) {
  double* mij_1d = d4est_operators_fetch_mij_1d(d4est_ops, deg);
  d4est_linalg_invert_and_copy(mij_1d, invmij_1d, deg + 1);
}

static void d4est_operators_build_dij_1d(d4est_operators_t* d4est_ops, double* dij_1d,
                                int deg) {
  int edge_nodes = deg + 1;
  double* v1d = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);
  double* inv_v1d = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);
  double* v1dr = (double*)P4EST_ALLOC(double, edge_nodes* edge_nodes);

  d4est_operators_build_Vij_1d(d4est_ops, v1d, deg);
  d4est_operators_build_dVij_1d(d4est_ops, v1dr, deg);
  d4est_linalg_invert_and_copy(v1d, inv_v1d, edge_nodes);

  d4est_linalg_mat_multiply(v1dr, inv_v1d, dij_1d, edge_nodes, edge_nodes,
                      edge_nodes);

  P4EST_FREE(v1dr);
  P4EST_FREE(inv_v1d);
  P4EST_FREE(v1d);
}




static double* d4est_operators_fetch_invmij_1d(d4est_operators_t* d4est_ops, int deg){

  int size = (deg + 1) * (deg + 1);
  return d4est_operators_1index_fetch
    (
     d4est_ops,
     d4est_ops->invmij_1d_table,
     deg,
     size,
     d4est_operators_build_invmij_1d
    );    
  
}

void d4est_operators_apply_mij(d4est_operators_t* d4est_ops, double* in, int dim, int deg,
                      double* out) {
  D4EST_ASSERT(dim == 1 || dim == 2 || dim == 3);
  double* mass_1d = d4est_operators_fetch_mij_1d(d4est_ops, deg);

  int nodes = deg + 1;
  if (dim == 1)
    d4est_linalg_matvec_plus_vec(1.0, mass_1d, in, 0., out, nodes, nodes);
  else if (dim == 2) {
    d4est_kron_A1A2x_nonsqr(out, mass_1d, mass_1d, in, nodes, nodes,
                             nodes, nodes);
  } else if (dim == 3) {
    d4est_kron_A1A2A3x_nonsqr(out, mass_1d, mass_1d, mass_1d, in,
                               nodes, nodes, nodes, nodes, nodes, nodes);
  } else {
    D4EST_ABORT("[D4EST_ERROR]: Apply mass matrix ref space, wrong dimension.");
  }
}

void d4est_operators_apply_invmij(d4est_operators_t* d4est_ops, double* in, int dim,
                         int deg, double* out) {
  D4EST_ASSERT(dim == 1 || dim == 2 || dim == 3);
  double* inv_mass_1d = d4est_operators_fetch_invmij_1d(d4est_ops, deg);

  int nodes = deg + 1;
  if (dim == 1)
    d4est_linalg_matvec_plus_vec(1.0, inv_mass_1d, in, 0., out, nodes, nodes);
  else if (dim == 2) {
    d4est_kron_A1A2x_nonsqr(out, inv_mass_1d, inv_mass_1d, in, nodes,
                             nodes, nodes, nodes);
  } else if (dim == 3) {
    d4est_kron_A1A2A3x_nonsqr(out, inv_mass_1d, inv_mass_1d,
                               inv_mass_1d, in, nodes, nodes, nodes, nodes,
                               nodes, nodes);
  } else {
    D4EST_ABORT("[D4EST_ERROR]: Apply mass matrix ref space, wrong dimension.");
  }
}

static double* d4est_operators_fetch_dij_1d(d4est_operators_t* d4est_ops, int deg) {
  int size = (deg + 1) * (deg + 1);
  return d4est_operators_1index_fetch
    (
     d4est_ops,
     d4est_ops->dij_1d_table,
     deg,
     size,
     d4est_operators_build_dij_1d
    );     
}



static void d4est_operators_build_hp_prolong_1d_aux(int degH, int degh, int c,
                                           double* inv_v1d_trans_degH,
                                           double* lobatto_nodes_degh,
                                           double* hp_prolong_matrix_1d) {
  D4EST_ASSERT(degH <= degh);
  int n, i;
  int nodesH = degH + 1;
  int nodesh = degh + 1;
  double* phi_modal = (double*)P4EST_ALLOC(double, nodesH);
  double r;
  for (n = 0; n < nodesh; n++) {
    r = lobatto_nodes_degh[n];
    if (c == 0 || c == 1) d4est_reference_child_to_parent_coords(c, &r);
    for (i = 0; i < nodesH; i++) {
      phi_modal[i] = d4est_lgl_jacobi(r, 0, 0, i);
    }

    d4est_linalg_matvec_plus_vec(1.0, inv_v1d_trans_degH, phi_modal, 0.0,
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
  d4est_linalg_invert_and_copy(v1d, inv_v1d, edge_nodes_degH);
  d4est_linalg_mat_transpose(inv_v1d, inv_v1d_trans_degH, edge_nodes_degH);
  double* lobatto_nodes_degh = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, degh);

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

  double* degh_nodes = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, degh);
  double* invvij = d4est_operators_fetch_invvij_1d(d4est_ops, degH);
  double* Vij_degh = P4EST_ALLOC(double, (degH+1)*(degh+1));
  
  int rows = degh + 1;
  int cols = degH + 1;
  
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      Vij_degh[i * cols + j] = d4est_lgl_jacobi(degh_nodes[i], 0., 0., j);

  d4est_linalg_mat_multiply(Vij_degh, invvij, p_prolong_1d, degh+1, degH+1, degH+1);
  P4EST_FREE(Vij_degh);
}




void d4est_operators_hp_apply_nd_prolong_transpose_with_ptr(
    double* uH, int degH, double* uh, int degh, int dim, int c,
    double* hp_prolong_transpose_matrix_1d) {
  /* sanity check */
  D4EST_ASSERT((degH <= degh) && (c < (1 << dim)));

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  if (dim == 1)
    d4est_linalg_matvec_plus_vec(1.0,
                           &hp_prolong_transpose_matrix_1d[c * nodesh * nodesH],
                           uh, 0., uH, nodesH, nodesh);
  else if (dim == 2) {
    int cx = d4est_reference_is_child_left_or_right(c, 0);
    int cy = d4est_reference_is_child_left_or_right(c, 1);


    d4est_kron_A1A2x_nonsqr(
        uH, &hp_prolong_transpose_matrix_1d[cy * nodesh * nodesH],
        &hp_prolong_transpose_matrix_1d[cx * nodesh * nodesH], uh, nodesH,
        nodesh, nodesH, nodesh);

  } else if (dim == 3) {
    int cx = d4est_reference_is_child_left_or_right(c, 0);
    int cy = d4est_reference_is_child_left_or_right(c, 1);
    int cz = d4est_reference_is_child_left_or_right(c, 2);


    d4est_kron_A1A2A3x_nonsqr(
        uH, &hp_prolong_transpose_matrix_1d[cz * nodesh * nodesH],
        &hp_prolong_transpose_matrix_1d[cy * nodesh * nodesH],
        &hp_prolong_transpose_matrix_1d[cx * nodesh * nodesH], uh, nodesH,
        nodesh, nodesH, nodesh, nodesH, nodesh);

  } else {
    D4EST_ABORT("[D4EST_ERROR]: d4est_operators_prolong_transpose_nd_U\n");
  }
}

static double* d4est_operators_fetch_hp_prolong_1d
(
 d4est_operators_t* d4est_ops, int degH,
 int degh)
{
  int size = 2 * (degh + 1) * (degH + 1);
  return d4est_operators_2index_fetch(d4est_ops,
                                      d4est_ops->hp_prolong_1d_table,
                                      degH,
                                      degh,
                                      size,
                                      d4est_operators_build_hp_prolong_1d
                                     );
  
  
}

double* d4est_operators_fetch_p_prolong_1d
(
 d4est_operators_t* d4est_ops,
 int degH,
 int degh
)
{
  int size = (degh + 1) * (degH + 1);
  return d4est_operators_2index_fetch(d4est_ops,
                                      d4est_ops->p_prolong_1d_table,
                                      degH,
                                      degh,
                                      size,
                                      d4est_operators_build_p_prolong_1d
                                     );
}

void d4est_operators_apply_hp_prolong(d4est_operators_t* d4est_ops, double* in, int degH,
                             int dim, int* degh, double* out) {
  D4EST_ASSERT(dim == 1 || dim == 2 || dim == 3);
  int children = (1 << dim);
  int c;
  double* hp_prolong_1d;

  int stride = 0;
  for (c = 0; c < children; c++) {
    hp_prolong_1d = d4est_operators_fetch_hp_prolong_1d(d4est_ops, degH, degh[c]);
    d4est_operators_hp_apply_nd_prolong_with_ptr(&out[stride], degh[c], in, degH, dim, c,
                                        hp_prolong_1d);
    stride += d4est_lgl_get_nodes(dim, degh[c]);
  }
}

void d4est_operators_apply_p_prolong(d4est_operators_t* d4est_ops, double* in, int degH,
                            int dim, int degh, double* out) {
  double* p_prolong_1d = d4est_operators_fetch_p_prolong_1d(d4est_ops, degH, degh);

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  /* if degrees are the same, p_prolong matrix should be the identity */
  if (degh == degH) {
    d4est_linalg_copy_1st_to_2nd(in, out, d4est_lgl_get_nodes(dim, degh));
    return;
  }

  if (dim == 1) {
    d4est_linalg_matvec_plus_vec(1.0, p_prolong_1d, in, 0., out, nodesh, nodesH);
  } else if (dim == 2) {
    d4est_kron_A1A2x_nonsqr(out, p_prolong_1d, p_prolong_1d, in, nodesh,
                             nodesH, nodesh, nodesH);
  } else if (dim == 3) {
    d4est_kron_A1A2A3x_nonsqr(out, p_prolong_1d, p_prolong_1d, p_prolong_1d,
                               in, nodesh, nodesH, nodesh, nodesH, nodesh,
                               nodesH);
  } else {
    D4EST_ABORT("[D4EST_ERROR]: d4est_operators_prolong_nd_U\n");
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
    d4est_linalg_matvec_plus_vec(.5, hp_prolong_matrix_1d, tmp, 0.0, &c_s[s * nodesh],
                           nodesh, nodesH);
    tmp[s] = 0.0;
  }

  d4est_linalg_mat_multiply(c_s, mass_matrix_rs_degh, c_s_x_Mh, nodesH, nodesh,
                      nodesh);

  d4est_linalg_mat_multiply(inv_mass_matrix_rs_degH, c_s_x_Mh,
                      &hp_restrict_matrix_1d[0], nodesH, nodesH, nodesh);

  P4EST_FREE(tmp);
  P4EST_FREE(c_s);
  P4EST_FREE(c_s_x_Mh);
}

static void d4est_operators_build_p_restrict_1d(d4est_operators_t* d4est_ops,
                                       double* p_restrict_1d, int degH,
                                       int degh) {
  int edge_nodes_degH = degH + 1;
  int edge_nodes_degh = degh + 1;

  double* p_prolong_1d = d4est_operators_fetch_p_prolong_1d(d4est_ops, degH, degh);
  double* inv_mij_degH =
      (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH);

  double* mij_degH = d4est_operators_fetch_mij_1d(d4est_ops, degH);
  d4est_linalg_invert_and_copy(mij_degH, inv_mij_degH, edge_nodes_degH);

  double* mij_degh = d4est_operators_fetch_mij_1d(d4est_ops, degh);
  d4est_operators_build_hp_restrict_1d_aux(degh, degH, &p_prolong_1d[0], mij_degh,
                                  inv_mij_degH, &p_restrict_1d[0]);

  d4est_linalg_vec_scale(2., p_restrict_1d, edge_nodes_degH * edge_nodes_degh);

  P4EST_FREE(inv_mij_degH);
}

double* d4est_operators_fetch_p_restrict_1d
(
 d4est_operators_t* d4est_ops,
 int degH,
 int degh
)
{
  int size = (degh + 1) * (degH + 1);
  return d4est_operators_2index_fetch(d4est_ops,
                                      d4est_ops->p_restrict_1d_table,
                                      degH,
                                      degh,
                                      size,
                                      d4est_operators_build_p_restrict_1d
                                     );

}

void d4est_operators_apply_p_restrict(d4est_operators_t* d4est_ops, double* in, int degh,
                             int dim, int degH, double* out) {
  double* p_restrict_1d = d4est_operators_fetch_p_restrict_1d(d4est_ops, degH, degh);

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  /* restrict matrix should be the identity if degrees are the same */
  if (degh == degH) {
    d4est_linalg_copy_1st_to_2nd(in, out, d4est_lgl_get_nodes(dim, degh));
    return;
  }

  if (dim == 1)
    d4est_linalg_matvec_plus_vec(1.0, p_restrict_1d, in, 0., out, nodesH, nodesh);
  else if (dim == 2) {
    d4est_kron_A1A2x_nonsqr(out, p_restrict_1d, p_restrict_1d, in, nodesH,
                             nodesh, nodesH, nodesh);
  } else if (dim == 3) {
    d4est_kron_A1A2A3x_nonsqr(out, p_restrict_1d, p_restrict_1d, p_restrict_1d,
                               in, nodesH, nodesh, nodesH, nodesh, nodesH,
                               nodesh);
  } else {
    D4EST_ABORT("[D4EST_ERROR]: d4est_operators_restrict_nd_U\n");
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

  double* inv_mij_degH =
      (double*)P4EST_ALLOC(double, edge_nodes_degH* edge_nodes_degH);

  /* TODO: replace this with a call to the database for inverse mij */
  double* mij_degH = d4est_operators_fetch_mij_1d(d4est_ops, degH);
  d4est_linalg_invert_and_copy(mij_degH, inv_mij_degH, edge_nodes_degH);

  double* mij_degh = d4est_operators_fetch_mij_1d(d4est_ops, degh);
  d4est_operators_build_hp_restrict_1d_aux(degh, degH, &hp_prolong_1d[0], mij_degh,
                                  inv_mij_degH, &hp_restrict_1d[0]);

  d4est_operators_build_hp_restrict_1d_aux(
      degh, degH, &hp_prolong_1d[edge_nodes_degH * edge_nodes_degh], mij_degh,
      inv_mij_degH, &hp_restrict_1d[edge_nodes_degH * edge_nodes_degh]);

  P4EST_FREE(inv_mij_degH);
}

static double* d4est_operators_fetch_hp_restrict_1d(d4est_operators_t* d4est_ops, int degH,
                                           int degh) {

  int size = 2 * (degh + 1) * (degH + 1);
  return d4est_operators_2index_fetch(d4est_ops,
                                      d4est_ops->hp_restrict_1d_table,
                                      degH,
                                      degh,
                                      size,
                                      d4est_operators_build_hp_restrict_1d
                                     );
  
}

void d4est_operators_apply_hp_restrict(d4est_operators_t* d4est_ops, double* in, int* degh,
                              int dim, int degH, double* out) {
  D4EST_ASSERT(dim == 1 || dim == 2 || dim == 3);
  int children = (1 << dim);
  int c;
  int nodesH = d4est_lgl_get_nodes(dim, degH);
  double* hp_restrict_1d;

  int stride = 0;
  double* tmp = P4EST_ALLOC(double, nodesH);
  d4est_linalg_fill_vec(out, 0., nodesH);

  for (c = 0; c < children; c++) {
    hp_restrict_1d = d4est_operators_fetch_hp_restrict_1d(d4est_ops, degH, degh[c]);
    d4est_operators_hp_apply_nd_restrict_with_ptr(tmp, degH, &in[stride], degh[c], dim,
                                         c, hp_restrict_1d);

    d4est_linalg_vec_axpy(1.0, tmp, out, nodesH);
    stride += d4est_lgl_get_nodes(dim, degh[c]);
  }

  P4EST_FREE(tmp);
}

static void d4est_operators_build_lobatto_rst_nd(d4est_operators_t* d4est_ops, double* ref_xyz_nd,
                                int dim, int deg) {
  int nodes = deg + 1;
  int vol_nodes = d4est_lgl_get_nodes(dim, deg);
  double* lgl = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, deg);
  int i;

  double* eye = P4EST_ALLOC_ZERO(double, nodes);
  for (i = 0; i < nodes; i++) eye[i] = 1.;

  if (dim == 2) {


    d4est_kron_AoB(eye, lgl, &ref_xyz_nd[0], (nodes), 1, (nodes), 1);
    d4est_kron_AoB(lgl, eye, &ref_xyz_nd[nodes * nodes], (nodes), 1, (nodes), 1);

  } else if (dim == 3) {


    d4est_kron_AoBoC(eye, eye, lgl, &ref_xyz_nd[0], nodes, 1, nodes, 1, nodes,
                      1);
    d4est_kron_AoBoC(eye, lgl, eye, &ref_xyz_nd[vol_nodes], nodes, 1, nodes, 1,
                      nodes, 1);
    d4est_kron_AoBoC(lgl, eye, eye, &ref_xyz_nd[2 * vol_nodes], nodes, 1,
                      nodes, 1, nodes, 1);

  } else
    D4EST_ABORT("ERROR 1: d4est_operators_build_lobatto_rst_nd");

  P4EST_FREE(eye);
}


static void d4est_operators_build_gauss_rst_nd(d4est_operators_t* d4est_ops, double* ref_xyz_nd,
                                int dim, int deg) {
  int nodes = deg + 1;
  int vol_nodes = d4est_lgl_get_nodes(dim, deg);
  double* gl = d4est_operators_fetch_gauss_nodes_1d(d4est_ops, deg);
  int i;

  double* eye = P4EST_ALLOC_ZERO(double, nodes);
  for (i = 0; i < nodes; i++) eye[i] = 1.;

  if (dim == 2) {

    d4est_kron_AoB(eye, gl, &ref_xyz_nd[0], (nodes), 1, (nodes), 1);
    d4est_kron_AoB(gl, eye, &ref_xyz_nd[nodes * nodes], (nodes), 1, (nodes), 1);

  } else if (dim == 3) {


    d4est_kron_AoBoC(eye, eye, gl, &ref_xyz_nd[0], nodes, 1, nodes, 1, nodes,
                      1);
    d4est_kron_AoBoC(eye, gl, eye, &ref_xyz_nd[vol_nodes], nodes, 1, nodes, 1,
                      nodes, 1);
    d4est_kron_AoBoC(gl, eye, eye, &ref_xyz_nd[2 * vol_nodes], nodes, 1,
                      nodes, 1, nodes, 1);

  } else
    D4EST_ABORT("ERROR 1: d4est_operators_build_lobatto_rst_nd");

  P4EST_FREE(eye);
}

double* d4est_operators_fetch_lobatto_rst_nd(d4est_operators_t* d4est_ops, int dim, int deg,
                            int dir)
{
  D4EST_ASSERT(dir < dim && dir >= 0 && dim <= 3 && dim > 0);
  
  if (dim == 1)
    return d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, deg);

  int size = dim*d4est_lgl_get_nodes(dim, deg); 
  double* op = d4est_operators_1index_2d_3d_fetch
    (
     d4est_ops,
     deg,
     dim,
     size,
     d4est_ops->lobatto_rst_2d_table,
     d4est_ops->lobatto_rst_3d_table,
     d4est_operators_build_lobatto_rst_nd
    );
   return &op[dir * d4est_lgl_get_nodes(dim, deg)];
}

void d4est_operators_apply_dij
(
 d4est_operators_t* d4est_ops,
 double* in,
 int dim,
 int deg,
 int dir,
 double* out
)
{
  D4EST_ASSERT(dir < dim && dir >= 0 && (dim == 1 || dim == 3 || dim == 2));
  double* Dr_1d = d4est_operators_fetch_dij_1d(d4est_ops, deg);
  int nodes = deg + 1;
  if (dim == 1)
    d4est_linalg_matvec_plus_vec(1.0, Dr_1d, in, 0., out, nodes, nodes);
  else if (dim == 2) {
    if (dir == 0) d4est_kron_IoMATx_SQR(out, Dr_1d, in, nodes);
    if (dir == 1) d4est_kron_MAToIx_SQR(out, Dr_1d, in, nodes);
  } else if (dim == 3) {
    if (dir == 0) d4est_kron_IoIoMATx_SQR(out, Dr_1d, in, nodes);
    if (dir == 1) d4est_kron_IoMAToIx_SQR(out, Dr_1d, in, nodes);
    if (dir == 2) d4est_kron_MAToIoIx_SQR(out, Dr_1d, in, nodes);
  } else {
    D4EST_ABORT("ERROR: d4est_operators_apply_dij");
  }
}

double* d4est_operators_fetch_gauss_rst_nd(d4est_operators_t* d4est_ops, int dim, int deg,
                            int dir) {

  D4EST_ASSERT(dir < dim && dir >= 0 && dim <= 3 && dim > 0);
  
  if (dim == 1)
    return d4est_operators_fetch_gauss_nodes_1d(d4est_ops, deg);

  int size = dim*d4est_lgl_get_nodes(dim, deg); 
  double* op = d4est_operators_1index_2d_3d_fetch
    (
     d4est_ops,
     deg,
     dim,
     size,
     d4est_ops->gauss_rst_2d_table,
     d4est_ops->gauss_rst_3d_table,
     d4est_operators_build_gauss_rst_nd
    );
   return &op[dir * d4est_lgl_get_nodes(dim, deg)];
}


static void d4est_operators_build_lift_1d(d4est_operators_t* d4est_ops, double* lift_1d, int deg) {
  memset(lift_1d, 0., sizeof(double) * 2 * (deg + 1));
  lift_1d[0] = 1.;
  lift_1d[2 * deg + 1] = 1.;
}

static double* d4est_operators_fetch_lift_1d(d4est_operators_t* d4est_ops, int deg) {

  int size = 2 * (deg + 1);
  return d4est_operators_1index_fetch
    (
     d4est_ops,
     d4est_ops->lift_1d_table,
     deg,
     size,
     d4est_operators_build_lift_1d
    );    
  
}

void d4est_operators_apply_lift(d4est_operators_t* d4est_ops, double* in, int dim, int deg,
                       int face, double* out) {
  D4EST_ASSERT(dim == 2 || dim == 3);
  D4EST_ASSERT(face < 2 * (dim));

  double* lift_1d = d4est_operators_fetch_lift_1d(d4est_ops, deg);

  /* int nodes = d4est_util_pow_int(deg+1, DIM); */

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
    D4EST_ABORT("ERROR 0: d4est_operators_lift_boundary_vec");
  }

  if (dim == 2){
    
    if (dir == 0)
      d4est_kron_IoVECx_SQR(out, &lift_1d[side * nodes], in, nodes);
    else if (dir == 1)
      d4est_kron_VECoIx_SQR(out, &lift_1d[side * nodes], in, nodes);
    else 
      D4EST_ABORT("dim == 2 so dir == 0 OR 1");
    
  }


  else if (dim == 3){

  if (dir == 0)
    d4est_kron_IoIoVECx_SQR(out, &lift_1d[side * nodes], in, nodes);
  else if (dir == 1)
    d4est_kron_IoVECoIx_SQR(out, &lift_1d[side * nodes], in, nodes);
  else if (dir == 2)
    d4est_kron_VECoIoIx_SQR(out, &lift_1d[side * nodes], in, nodes);
  else
    D4EST_ABORT("DIM = 3 so DIR = 0,1,2");
  
  }

  else {
    D4EST_ABORT("DIM not supported");
  }
}

void d4est_operators_apply_slicer(d4est_operators_t* d4est_ops, double* in, int dim,
                         int face, int deg, double* out) {
  D4EST_ASSERT(face < 2 * (dim));
  /* int nodes = d4est_util_pow_int(deg+1, DIM); */

  double* slicer_1d = d4est_operators_fetch_lift_1d(d4est_ops, deg);

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
    D4EST_ABORT("ERROR 0: d4est_operators_lift_boundary_vec");
  }

  /* d4est_util_print_matrix(slicer_1d, (deg+1), 1, "slicer_1d = ", 0); */

  if (dim == 2){

    if (dir == 0)
      d4est_kron_IoVEC_TRANSx_SQR(out, &slicer_1d[side * (deg + 1)], in,
                                 deg + 1);
    else if (dir == 1)
      d4est_kron_VEC_TRANSoIx_SQR(out, &slicer_1d[side * (deg + 1)], in,
                                 deg + 1);
    else
      D4EST_ABORT("DIM == 2, so dir == 0,1");
  }

  else if (dim == 3){

    if (dir == 0)
      d4est_kron_IoIoVEC_TRANSx_SQR(out, &slicer_1d[side * (deg + 1)], in, deg + 1);
    else if (dir == 1)
      d4est_kron_IoVEC_TRANSoIx_SQR(out, &slicer_1d[side * (deg + 1)], in, deg + 1);
    else if (dir == 2)
      d4est_kron_VEC_TRANSoIoIx_SQR(out, &slicer_1d[side * (deg + 1)], in,
                                   deg + 1);
    else
      D4EST_ABORT("DIM == 3 so DIR=0,1 or 2");
  
  }
  else {
    D4EST_ABORT("DIM must be 2 or 3");
  }
}


void d4est_operators_convert_nodal_to_modal(d4est_operators_t* d4est_ops, double* in,
                                   int dim, int deg, double* out) {
  D4EST_ASSERT(dim == 1 || dim == 2 || dim == 3);
  int nodes = deg + 1;
  double* invvij_1d = P4EST_ALLOC(double, nodes* nodes);

  /* TODO: probably could use build function for inv_Vij */
  d4est_operators_build_Vij_1d(d4est_ops, invvij_1d, deg);
  d4est_linalg_invert(invvij_1d, nodes);

  if (dim == 1)
    d4est_linalg_matvec_plus_vec(1.0, invvij_1d, in, 0., out, nodes, nodes);
  else if (dim == 2) {
    d4est_kron_A1A2x_nonsqr(out, invvij_1d, invvij_1d, in, nodes, nodes, nodes,
                             nodes);
  } else if (dim == 3) {
    d4est_kron_A1A2A3x_nonsqr(out, invvij_1d, invvij_1d, invvij_1d, in, nodes,
                               nodes, nodes, nodes, nodes, nodes);
  } else {
    D4EST_ABORT("[D4EST_ERROR]: Apply mass matrix ref space, wrong dimension.");
  }

  P4EST_FREE(invvij_1d);
}

static void d4est_operators_build_hp_prolong_transpose_1d(
    d4est_operators_t* d4est_ops, double* hp_prolong_transpose_1d, int degH,
    int degh) {
  double* hp_prolong_1d = d4est_operators_fetch_hp_prolong_1d(d4est_ops, degH, degh);
  d4est_linalg_mat_transpose_nonsqr(hp_prolong_1d, hp_prolong_transpose_1d,
                              (degh + 1), (degH + 1));
  d4est_linalg_mat_transpose_nonsqr(&hp_prolong_1d[(degh + 1) * (degH + 1)],
                              &hp_prolong_transpose_1d[(degh + 1) * (degH + 1)],
                              (degh + 1), (degH + 1));
}


static double* d4est_operators_fetch_hp_prolong_transpose_1d(d4est_operators_t* d4est_ops,
                                                    int degH, int degh) {

  int size = 2 * (degh + 1) * (degH + 1);
  return d4est_operators_2index_fetch(d4est_ops,
                                      d4est_ops->hp_prolong_transpose_1d_table,
                                      degH,
                                      degh,
                                      size,
                                      d4est_operators_build_hp_prolong_transpose_1d
                                     );
  
}





void d4est_operators_build_p_prolong_1d_inverse(d4est_operators_t* d4est_ops,
                                      double* p_prolong_1d_inverse, int degH,
                                      int degh)
{
  double* p_prolong_1d = d4est_operators_fetch_p_prolong_1d(d4est_ops, degH, degh);
  d4est_linalg_leftinverse(p_prolong_1d, p_prolong_1d_inverse, degh + 1, degH + 1);
}





static void d4est_operators_build_p_prolong_transpose_1d(d4est_operators_t* d4est_ops,
                                                double* p_prolong_transpose_1d,
                                                int degH, int degh) {
  double* p_prolong_1d = d4est_operators_fetch_p_prolong_1d(d4est_ops, degH, degh);
  d4est_linalg_mat_transpose_nonsqr(p_prolong_1d, p_prolong_transpose_1d, (degh + 1),
                              (degH + 1));
}


double* d4est_operators_fetch_p_prolong_transpose_1d(d4est_operators_t* d4est_ops,
                                                   int degH, int degh) {

  int size = (degh + 1) * (degH + 1);
  return d4est_operators_2index_fetch(d4est_ops,
                                      d4est_ops->p_prolong_transpose_1d_table,
                                      degH,
                                      degh,
                                      size,
                                      d4est_operators_build_p_prolong_transpose_1d
                                     );
 
}


void d4est_operators_build_p_prolong_transpose_1d_inverse(d4est_operators_t* d4est_ops,
                                                double* p_prolong_transpose_1d_inverse,
                                                int degH, int degh) {
  double* p_prolong_transpose_1d = d4est_operators_fetch_p_prolong_transpose_1d(d4est_ops, degH, degh);
  d4est_linalg_leftinverse(p_prolong_transpose_1d, p_prolong_transpose_1d_inverse, degH + 1, degh + 1);
}


void d4est_operators_apply_hp_prolong_transpose(d4est_operators_t* d4est_ops, double* in,
                                       int* degh, int dim, int degH,
                                       double* out) {
  D4EST_ASSERT(dim == 1 || dim == 2 || dim == 3);
  int children = (1 << dim);
  int c;
  int nodesH = d4est_lgl_get_nodes(dim, degH);
  double* hp_prolong_transpose_1d;

  int stride = 0;
  double* tmp = P4EST_ALLOC(double, nodesH);
  d4est_linalg_fill_vec(out, 0., nodesH);

  for (c = 0; c < children; c++) {
    hp_prolong_transpose_1d =
        d4est_operators_fetch_hp_prolong_transpose_1d(d4est_ops, degH, degh[c]);
    /* int nodesh = d4est_lgl_get_nodes(dim,degh[c]); */
    d4est_operators_hp_apply_nd_prolong_transpose_with_ptr(
        tmp, degH, &in[stride], degh[c], dim, c, hp_prolong_transpose_1d);
    /* printf("child = %d\n",c); */
    /* d4est_util_print_matrix(tmp, nodesH, 1, "tmp = ", 0); */
    /* d4est_util_print_matrix(&uh[stride], nodesh, 1, "uh[stride] = ", 0); */
    d4est_linalg_vec_axpy(1.0, tmp, out, nodesH);
    /* d4est_util_print_matrix(uH, nodesH, 1, "uH = ", 0); */
    stride += d4est_lgl_get_nodes(dim, degh[c]);
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
    d4est_linalg_copy_1st_to_2nd(in, out, d4est_lgl_get_nodes(dim, degh));
    return;
  }

  if (dim == 1)
    d4est_linalg_matvec_plus_vec(1.0, p_prolong_transpose_1d, in, 0., out, nodesH,
                           nodesh);
  else if (dim == 2) {
    d4est_kron_A1A2x_nonsqr(out, p_prolong_transpose_1d,
                             p_prolong_transpose_1d, in, nodesH, nodesh, nodesH,
                             nodesh);
  } else if (dim == 3) {
    d4est_kron_A1A2A3x_nonsqr(out, p_prolong_transpose_1d,
                               p_prolong_transpose_1d, p_prolong_transpose_1d,
                               in, nodesH, nodesh, nodesH, nodesh, nodesH,
                               nodesh);
  } else {
    D4EST_ABORT("[D4EST_ERROR]: d4est_operators_apply_p_prolong_transpose\n");
  }
}

static
void d4est_operators_build_hp_restrict_interp_1d_aux(
    int degh, int degH, int c, double* inv_v1d_trans_degh,
    double* lobatto_nodes_degH, double* hp_restrict_interp_matrix_1d) {
  D4EST_ASSERT(degH <= degh);
  D4EST_ASSERT(c == 0 || c == 1);

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
      phi_modal[i] = d4est_lgl_jacobi(r, 0, 0, i);
    }

    if (n == degH / 2 && degH % 2 == 0)
      d4est_linalg_matvec_plus_vec(.5, inv_v1d_trans_degh, phi_modal, 0.0,
                             &hp_restrict_interp_matrix_1d[n * nodesh], nodesh,
                             nodesh);
    else if (r < -1. || r > 1.)
      d4est_linalg_fill_vec(&hp_restrict_interp_matrix_1d[n * nodesh], 0., nodesh);
    else
      d4est_linalg_matvec_plus_vec(1., inv_v1d_trans_degh, phi_modal, 0.0,
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
  d4est_linalg_invert_and_copy(v1d, inv_v1d, edge_nodes_degh);
  d4est_linalg_mat_transpose(inv_v1d, inv_v1d_trans_degh, edge_nodes_degh);

  double* lobatto_nodes_degH = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, degH);
  d4est_operators_build_hp_restrict_interp_1d_aux(degh, degH, 0, inv_v1d_trans_degh,
                                         lobatto_nodes_degH,
                                         &hp_restrict_interp_1d[0]);
  d4est_operators_build_hp_restrict_interp_1d_aux(
      degh, degH, 1, inv_v1d_trans_degh, lobatto_nodes_degH,
      &hp_restrict_interp_1d[edge_nodes_degH * edge_nodes_degh]);

  P4EST_FREE(inv_v1d);
  P4EST_FREE(inv_v1d_trans_degh);
  P4EST_FREE(v1d);
}

static double* d4est_operators_fetch_hp_restrict_interp_1d(d4est_operators_t* d4est_ops,
                                                  int degH, int degh) {

  int size = 2 * (degh + 1) * (degH + 1);
  return d4est_operators_2index_fetch(d4est_ops,
                                      d4est_ops->hp_restrict_interp_1d_table,
                                      degH,
                                      degh,
                                      size,
                                      d4est_operators_build_hp_restrict_interp_1d
                                     );
}

void d4est_operators_apply_hp_restrict_interp(d4est_operators_t* d4est_ops, double* in,
                                     int* degh, int dim, int degH,
                                     double* out) {
  D4EST_ASSERT(dim == 1 || dim == 2 || dim == 3);
  int children = (1 << dim);
  int c;
  int nodesH = d4est_lgl_get_nodes(dim, degH);
  double* hp_restrict_1d;

  int stride = 0;
  double* tmp = P4EST_ALLOC(double, nodesH);
  d4est_linalg_fill_vec(out, 0., nodesH);

  for (c = 0; c < children; c++) {
    hp_restrict_1d = d4est_operators_fetch_hp_restrict_interp_1d(d4est_ops, degH, degh[c]);
    d4est_operators_hp_apply_nd_restrict_with_ptr(tmp, degH, &in[stride], degh[c], dim,
                                         c, hp_restrict_1d);
    d4est_linalg_vec_axpy(1.0, tmp, out, nodesH);
    stride += d4est_lgl_get_nodes(dim, degh[c]);
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
      phi_modal[i] = d4est_lgl_jacobi(r, 0, 0, i);
    }
    d4est_linalg_matvec_plus_vec(1., inv_v1d_trans_degh, phi_modal, 0.0,
                           &hp_restrict_interp_matrix_1d[n * nodesh], nodesh,
                           nodesh);
  }

  P4EST_FREE(phi_modal);
}

static void d4est_operators_build_p_restrict_interp_1d(d4est_operators_t* d4est_ops,
                                              double* hp_restrict_interp_1d,
                                              int degH, int degh) {
  D4EST_ASSERT(degH <= degh);
  /* int edge_nodes_degH = degH + 1; */
  int edge_nodes_degh = degh + 1;

  double* v1d = (double*)P4EST_ALLOC(double, edge_nodes_degh* edge_nodes_degh);
  double* inv_v1d =
      (double*)P4EST_ALLOC(double, edge_nodes_degh* edge_nodes_degh);
  double* inv_v1d_trans_degh =
      (double*)P4EST_ALLOC(double, edge_nodes_degh* edge_nodes_degh);

  d4est_operators_build_Vij_1d(d4est_ops, v1d, degh);
  d4est_linalg_invert_and_copy(v1d, inv_v1d, edge_nodes_degh);
  d4est_linalg_mat_transpose(inv_v1d, inv_v1d_trans_degh, edge_nodes_degh);

  double* lobatto_nodes_degH = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, degH);
  d4est_operators_build_p_restrict_interp_1d_aux(degh, degH, inv_v1d_trans_degh,
                                        lobatto_nodes_degH,
                                        &hp_restrict_interp_1d[0]);

  P4EST_FREE(inv_v1d);
  P4EST_FREE(inv_v1d_trans_degh);
  P4EST_FREE(v1d);
}

static double*
d4est_operators_fetch_p_restrict_interp_1d
(
 d4est_operators_t* d4est_ops,
 int degH,
 int degh
)
{
  int size = (degh + 1) * (degH + 1);
  return d4est_operators_2index_fetch(d4est_ops,
                                      d4est_ops->p_restrict_interp_1d_table,
                                      degH,
                                      degh,
                                      size,
                                      d4est_operators_build_p_restrict_interp_1d
                                     );  
}

void d4est_operators_apply_p_restrict_interp(d4est_operators_t* d4est_ops, double* in,
                                    int degh, int dim, int degH, double* out) {
  double* p_restrict_1d = d4est_operators_fetch_p_restrict_interp_1d(d4est_ops, degH, degh);

  int nodesh = degh + 1;
  int nodesH = degH + 1;

  /* restrict matrix should be the identity if degrees are the same */
  if (degh == degH) {
    d4est_linalg_copy_1st_to_2nd(in, out, d4est_lgl_get_nodes(dim, degh));
    return;
  }

  if (dim == 1)
    d4est_linalg_matvec_plus_vec(1.0, p_restrict_1d, in, 0., out, nodesH, nodesh);
  else if (dim == 2) {
    d4est_kron_A1A2x_nonsqr(out, p_restrict_1d, p_restrict_1d, in, nodesH,
                             nodesh, nodesH, nodesh);
  } else if (dim == 3) {
    d4est_kron_A1A2A3x_nonsqr(out, p_restrict_1d, p_restrict_1d, p_restrict_1d,
                               in, nodesH, nodesh, nodesH, nodesh, nodesH,
                               nodesh);
  } else {
    D4EST_ABORT("[D4EST_ERROR]: d4est_operators_restrict_nd_U\n");
  }
}

static void d4est_operators_build_flip_1d(d4est_operators_t* d4est_ops, double* flip_1d, int deg) {
  memset(flip_1d, 0., sizeof(double) * (deg + 1) * (deg + 1));
  int offset;
  int i;
  for (i = 0; i < deg + 1; i++) {
    offset = (i+1)*deg;
    flip_1d[offset] = 1.;
  }
}

static double* d4est_operators_fetch_flip_1d(d4est_operators_t* d4est_ops, int deg) {
  int size = (deg + 1) * (deg + 1);
  return d4est_operators_1index_fetch
    (
     d4est_ops,
     d4est_ops->flip_1d_table,
     deg,
     size,
     d4est_operators_build_flip_1d
    );     
}

void d4est_operators_apply_flip(d4est_operators_t* d4est_ops, double* in, int dim, int deg,
                       int dir, double* out) {
  D4EST_ASSERT(dim == 1 || dim == 2);
  double* flip_1d = d4est_operators_fetch_flip_1d(d4est_ops, deg);
  int nodes = deg + 1;
  if (dim == 1)
    d4est_linalg_matvec_plus_vec(1.0, flip_1d, in, 0., out, nodes, nodes);
  else if (dim == 2) {
    if (dir == 0) d4est_kron_IoMATx_SQR(out, flip_1d, in, nodes);
    if (dir == 1) d4est_kron_MAToIx_SQR(out, flip_1d, in, nodes);
    if (dir == 2)
      d4est_kron_A1A2x_nonsqr(out, flip_1d, flip_1d, in, nodes, nodes, nodes,
                               nodes);
  } else {
    D4EST_ABORT("[D4EST_ERROR]: flip not supported in this dimension atm.");
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
  D4EST_ASSERT((face_dim == 2 && o < 4) || (face_dim == 1 && o < 2));

  if (o == 4)
    printf("[WARNING]: o == 4 and this orientation has not been tested");
  
  int nodes = d4est_lgl_get_nodes(face_dim, deg);
  
  if (face_dim == 1){
    if (o == 1)
      d4est_operators_apply_flip
        (
         d4est_ops,
         in,
         face_dim,
         deg,
         0,
         out
        );
    else
      d4est_linalg_copy_1st_to_2nd(in, out, nodes);
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
    d4est_operators_apply_flip
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
    d4est_operators_apply_flip
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
    d4est_linalg_mat_transpose(tmp1, out, deg+1);
  }
  else {
    d4est_linalg_copy_1st_to_2nd(tmp1, out, nodes);
  }

  if (tmp0 != in)
    P4EST_FREE(tmp0);
  if(tmp1 != tmp0)
    P4EST_FREE(tmp1);
}


static void d4est_operators_build_vtk_interp_1d(d4est_operators_t* d4est_ops, double* vtk_interp_1d, int deg) {
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
  int size = (deg) * (deg + 1) * 2;
  return d4est_operators_1index_fetch
    (
     d4est_ops,
     d4est_ops->vtk_interp_1d_table,
     deg,
     size,
     d4est_operators_build_vtk_interp_1d
    );      
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
  D4EST_ASSERT((c < (1 << dim)));

  double* vtk_interp_1d = d4est_operators_fetch_vtk_interp_1d(d4est_ops, deg);
  
  int nodes_in = deg + 1;
  int nodes_out = deg;

  if (dim == 1)
    d4est_linalg_matvec_plus_vec(1.0, &vtk_interp_1d[c * nodes_in * nodes_out], in,
                           0., out, nodes_out, nodes_in);
  else if (dim == 2) {
    int cx = d4est_reference_is_child_left_or_right(c, 0);
    int cy = d4est_reference_is_child_left_or_right(c, 1);
    d4est_kron_A1A2x_nonsqr(out, &vtk_interp_1d[cx * nodes_in * nodes_out],
                             &vtk_interp_1d[cy * nodes_in * nodes_out], in,
                             nodes_out, nodes_in, nodes_out, nodes_in);
  } else if (dim == 3) {
    int cx = d4est_reference_is_child_left_or_right(c, 0);
    int cy = d4est_reference_is_child_left_or_right(c, 1);
    int cz = d4est_reference_is_child_left_or_right(c, 2);
    d4est_kron_A1A2A3x_nonsqr(out, &vtk_interp_1d[cx * nodes_in * nodes_out],
                               &vtk_interp_1d[cy * nodes_in * nodes_out],
                               &vtk_interp_1d[cz * nodes_in * nodes_out], in,
                               nodes_out, nodes_in, nodes_out, nodes_in, nodes_out, nodes_in);
  } else {
    D4EST_ABORT("Dim = 1 or 2 or 3 in vtk_interp");
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
      double* rst = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, dim, deg, dir);
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
  int children = (dim == 2) ? 4 : 8;
  int deg_dim = (dim == 2) ? deg*deg : deg*deg*deg;
  int size = 3*deg_dim*(children);

  return d4est_operators_1index_2d_3d_fetch
    (
     d4est_ops,
     deg,
     dim,
     size,
     d4est_ops->vtk_rst_2d_table,
     d4est_ops->vtk_rst_3d_table,
     d4est_operators_build_vtk_rst
    );
}


static void d4est_operators_build_dij_trans_1d(d4est_operators_t* d4est_ops, double* dij_trans_1d,
                                int deg) {

  double* Dr_1d = d4est_operators_fetch_dij_1d(d4est_ops, deg);
  /* double* Dr_1d_transpose = P4EST_ALLOC(double, (deg+1)*(deg+1)); */
  d4est_linalg_mat_transpose(Dr_1d, dij_trans_1d, (deg+1));
}


static double* d4est_operators_fetch_dij_trans_1d(d4est_operators_t* d4est_ops, int deg) {
  int size = (deg + 1) * (deg + 1);
  return d4est_operators_1index_fetch
    (
     d4est_ops,
     d4est_ops->dij_trans_1d_table,
     deg,
     size,
     d4est_operators_build_dij_trans_1d
    );     
}



void d4est_operators_apply_dij_transpose(d4est_operators_t* d4est_ops, double* in, int dim, int deg,
                      int dir, double* out) {
  D4EST_ASSERT(dir < dim && dir >= 0 && (dim == 1 || dim == 3 || dim == 2));
  double* Dr_1d_transpose = d4est_operators_fetch_dij_trans_1d(d4est_ops, deg);

  int nodes = deg + 1;
  if (dim == 1)
    d4est_linalg_matvec_plus_vec(1.0, Dr_1d_transpose, in, 0., out, nodes, nodes);
  else if (dim == 2) {

    if (dir == 0) d4est_kron_IoMATx_SQR(out, Dr_1d_transpose, in, nodes);
    if (dir == 1) d4est_kron_MAToIx_SQR(out, Dr_1d_transpose, in, nodes);
    
  } else if (dim == 3) {

    if (dir == 0) d4est_kron_IoIoMATx_SQR(out, Dr_1d_transpose, in, nodes);
    if (dir == 1) d4est_kron_IoMAToIx_SQR(out, Dr_1d_transpose, in, nodes);
    if (dir == 2) d4est_kron_MAToIoIx_SQR(out, Dr_1d_transpose, in, nodes);

  } else {
    /* P4EST_FREE(Dr_1d_transpose); */
    D4EST_ABORT("[D4EST_ERROR]: d4est_operators_apply_dij");
  }

  /* P4EST_FREE(Dr_1d_transpose); */
}
