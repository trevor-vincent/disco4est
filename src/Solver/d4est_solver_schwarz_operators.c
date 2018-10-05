#include <pXest.h>
#include <d4est_solver_schwarz_operators.h>

static
double quintic_poly
(double r)
{
  return (15.*r - 10.*r*r*r + 3.*r*r*r*r*r)/8.;
}

static
double phi_fcn
(
 double r
)
{
  if (r < -1 || r > 1){
    /* signum(r) */
    return (r > 0) - (r < 0); 
  }
  else {
    return quintic_poly(r);
  }
}

static
double poly_hat_weight_fcn
(
 double r,
 double overlap_size /* subdomain overlap size in rst space */
)
{
  double d0 = overlap_size;
  return .5*(phi_fcn((r+1)/d0) + phi_fcn((r-1)/d0));
}

void d4est_solver_schwarz_operators_build_schwarz_restrictor_1d
(
 d4est_operators_t* d4est_ops,
 double* restrict restrictor_1d,
 int deg, //degree of element before restriction
 int restricted_size //number of 1d nodes
)
{
  D4EST_ASSERT(restricted_size <= deg + 1 && restricted_size >= 1);
  int original_size = deg + 1;

  for (int i = 0; i < restricted_size; i++){
    /* face == 0 */
    restrictor_1d[(deg+1)*(i) + i] = 1.;
    /* face == 1 */
    restrictor_1d[original_size * restricted_size + (deg+1)*(i+1) - restricted_size + i] = 1.;
  }
}

double*
d4est_solver_schwarz_operators_fetch_schwarz_restrictor_1d
(
 d4est_solver_schwarz_operators_t* d4est_schwarz_ops,
 int deg,
 int restricted_size
)
{
  int size = 2 * (deg + 1) * (restricted_size);
  return d4est_operators_2index_fetch
    (
     d4est_schwarz_ops->d4est_ops,
     d4est_schwarz_ops->schwarz_restrictor_1d_table,
     deg,
     restricted_size,
     size,
     d4est_solver_schwarz_operators_build_schwarz_restrictor_1d
    );
}

void d4est_solver_schwarz_operators_build_schwarz_weights_1d
(
 d4est_operators_t* d4est_ops,
 double* restrict schwarz_weights_1d,
 int deg, //degree of element before restriction
 int restricted_size //number of 1d nodes
)
{
  double* r = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops,
                                                     deg);
  double rmax = 1.;
  double rmin = 1. - r[deg + 1 - restricted_size - 1];
  double overlap_size = rmax - rmin;

  for (int i = 0; i < restricted_size; i++){
    schwarz_weights_1d[i] = poly_hat_weight_fcn(r - 2, overlap_size);/* left subdomain element*/
    schwarz_weights_1d[restricted_size + i] = poly_hat_weight_fcn(r + 2, overlap_size);/* right subdomain element*/
  }
  for (int i = 0; i < deg + 1; i++){
    schwarz_weights_1d[2*restricted_size + i] = poly_hat_weight_fcn(r, overlap_size); /* core subdomain element*/
  }  
}

double*
d4est_solver_schwarz_operators_fetch_schwarz_weights_1d
(
 d4est_solver_schwarz_operators_t* d4est_schwarz_ops,
 int deg,
 int restricted_size
)
{
  int size = (deg + 1) + 2*restricted_size;
  return d4est_operators_2index_fetch
    (
     d4est_schwarz_ops->d4est_ops,
     d4est_schwarz_ops->schwarz_weights_1d_table,
     deg,
     restricted_size,
     size,
     d4est_solver_schwarz_operators_build_schwarz_weights_1d
    );
}


void d4est_solver_schwarz_operators_build_schwarz_restrictor_transpose_1d
(
 d4est_operators_t* d4est_ops,
 double* restrict restrictor_transpose_1d,
 int deg, //degree of element before restriction
 int res_size, //number of 1d nodes
 void* ctx
)
{

  d4est_solver_schwarz_operators_t* d4est_schwarz_ops = ctx;
  double* restrictor_1d = d4est_solver_schwarz_operators_fetch_schwarz_restrictor_1d(d4est_schwarz_ops,deg,res_size);
  d4est_linalg_mat_transpose_nonsqr(&restrictor_1d[0],
                                    &restrictor_transpose_1d[0],
                                    res_size, deg + 1);
  d4est_linalg_mat_transpose_nonsqr(&restrictor_1d[(res_size)*(deg+1)],
                                    &restrictor_transpose_1d[(res_size)*(deg+1)],
                                    res_size, deg + 1);
}

double*
d4est_solver_schwarz_operators_fetch_schwarz_restrictor_transpose_1d
(
 d4est_solver_schwarz_operators_t* d4est_schwarz_ops,
 int deg,
 int res_size
)
{
  int size = 2 * (deg + 1) * (res_size);
  return d4est_operators_2index_fetch_with_ctx
    (
     d4est_schwarz_ops->d4est_ops,
     d4est_schwarz_ops->schwarz_restrictor_transpose_1d_table,
     deg, 
     res_size,
     size,
     d4est_operators_build_schwarz_restrictor_transpose_1d,
     d4est_schwarz_ops
    );
}

void d4est_operators_apply_schwarz_restrictor
(
 d4est_solver_schwarz_operators_t* d4est_schwarz_ops,
 double *   in,
 int dim,
 int* faces,
 int deg,
 int restricted_size,
 d4est_ops_transpose_t transpose,
 double *   out
)
{  
  D4EST_ASSERT(restricted_size > 0 && restricted_size <= deg + 1); 
  double* schwarz_1d = NULL;
  int nodes_in = -1;
  int nodes_out = -1;
  
  if (transpose == D4OPS_NO_TRANSPOSE){
    nodes_in = deg + 1;
    nodes_out = restricted_size;
    schwarz_1d = d4est_operators_fetch_schwarz_restrictor_1d(d4est_schwarz_ops, deg, restricted_size);
  }
  else if (transpose == D4OPS_TRANSPOSE) {
    nodes_out = deg + 1;
    nodes_in = restricted_size;
    schwarz_1d = d4est_operators_fetch_schwarz_restrictor_transpose_1d(d4est_schwarz_ops, deg, restricted_size);
  }
  else {
    D4EST_ABORT("Transpose needs to be D4OPS_NO_TRANSPOSE or D4OPS_TRANSPOSE");
  }

  /* This should be optimized out */
  int degp1 = deg + 1;
  double* eyes = P4EST_ALLOC_ZERO(double, (degp1)*(degp1));
  for (int i = 0; i < (degp1); i++){
    eyes[i*(degp1) + i] = 1.;
  }

  int op_rows [3];
  int op_cols [3];
  double* operators [3];
  
  int dir;  /* x <-> 0; y <-> 1; z <-> 2  */
  int side; /* left <-> 0; right <-> 1 */
  
  for (int i = 0; i < dim; i++){
    operators[i] = eyes;
    op_rows[i] = degp1;
    op_cols[i] = degp1;
  }

  for (int i = 0; i < dim; i++){  
    if (faces[i] != -1) {
      d4est_reference_dir_and_side_of_face(faces[i], &dir, &side);
      operators[dir] = &schwarz_1d[side * nodes_in * nodes_out];
      op_rows[dir] = nodes_out;
      op_cols[dir] = nodes_in;
    }
  }
  
  if (dim == 2){
    d4est_kron_A1A2x_nonsqr(out,
                            operators[1],
                            operators[0],
                            in,
                            op_rows[1],
                            op_cols[1],
                            op_rows[0],
                            op_cols[0]
                           );
      
  }
  else if (dim == 3){
    d4est_kron_A1A2A3x_nonsqr(out,
                              operators[2],
                              operators[1],
                              operators[0],
                              in,
                              op_rows[2],
                              op_cols[2],
                              op_rows[1],
                              op_cols[1],
                              op_rows[0],
                              op_cols[0]
                             );
  }
  else {
    D4EST_ABORT("DIM must be 2 or 3");
  }

  P4EST_FREE(eyes);
}


/** 
 * Apply the weight function to a subdomain side element field.
 * The faces array is used to determine if the side element is
 * to the left or to the right of the core element. If faces is NULL
 * then we assume the "in" field comes from the core element.
 * 
 * @param d4est_ops 
 * @param in 
 * @param dim 
 * @param faces array of core faces that touch subdomain element side, NULL if core element
 * @param deg 
 * @param restricted_size 
 * @param out 
 */
void d4est_operators_apply_schwarz_weights
(
 d4est_solver_schwarz_operators_t* d4est_schwarz_ops,
 double *   in,
 int dim,
 int* faces,
 int deg,
 int restricted_size,
 double *   out
)
{  
  D4EST_ASSERT(restricted_size > 0 && restricted_size <= deg + 1); 

  /* This should be optimized out */
  int degp1 = deg + 1;
  double* eyes = P4EST_ALLOC_ZERO(double, (degp1));
  for (int i = 0; i < (degp1); i++){
    eyes[i] = 1.;
  }

  int vec_sizes [3];
  double* vecs [3];
  
  int dir;  /* x <-> 0; y <-> 1; z <-> 2  */
  int side; /* left <-> 0; right <-> 1 */

  double* weights_1d = d4est_solver_schwarz_operators_fetch_schwarz_weights(d4est_schwarz_ops, deg, restricted_size);

  /* populate with core weights if the core or 1's if not the core */
  for (int i = 0; i < dim; i++){
    vecs[i] = (faces == NULL) ? &weights_1d[2*restricted_size] : eyes;
    vec_sizes[i] = degp1;
  }


  /* if not the core populate vecs with side weights based on if side is left or right of core */
  if (faces != NULL){
    for (int i = 0; i < dim; i++){  
      if (faces[i] != -1) {
        d4est_reference_dir_and_side_of_face(faces[i], &dir, &side);
        vecs[dir] = &weights_1d[side * restricted_size];
        vec_rows[dir] = restricted_size;
      }
    }
  }

  if (dim == 1){
    d4est_kron_vec1_dot_x(vecs[0], in, vec_sizes[0], out);
  }
  else if (dim == 2){
    d4est_kron_vec1_o_vec2_dot_x(vecs[1], vecs[0], in, vec_sizes[1], vec_sizes[0], out);
  }
  else if (dim == 3){
    d4est_kron_vec1_o_vec2_o_vec3_dot_x(vecs[2], vecs[1], vecs[0], in, vec_sizes[2], vec_sizes[1], vec_sizes[0], out);
  }
  else {
    D4EST_ABORT("DIM must be 2 or 3");
  }

  P4EST_FREE(eyes);
}


d4est_solver_schwarz_operators_t*
d4est_solver_schwarz_operators_init
(
 d4est_operators_t* d4est_ops
)
{
  d4est_solver_schwarz_operators_t* d4est_schwarz_ops = P4EST_ALLOC(d4est_solver_schwarz_operators_t, 1);
  d4est_schwarz_ops->d4est_ops = d4est_ops;
  int max_degree = d4est_ops->max_degree;
  
  d4est_schwarz_ops->schwarz_restrictor_1d_table = P4EST_ALLOC(double**, max_degree);
  d4est_schwarz_ops->schwarz_weights_1d_table = P4EST_ALLOC(double**, max_degree);
  d4est_schwarz_ops->schwarz_restrictor_transpose_1d_table = P4EST_ALLOC(double**, max_degree);
  
  for (int i = 0; i < d4est_schwarz_ops->max_degree; i++){
    d4est_schwarz_ops->schwarz_restrictor_1d_table[i] = P4EST_ALLOC(double*, max_degree);
    d4est_schwarz_ops->schwarz_weights_1d_table[i] = P4EST_ALLOC(double*, max_degree);
    d4est_schwarz_ops->schwarz_restrictor_transpose_1d_table[i] = P4EST_ALLOC(double*, max_degree);

    for (int j = 0; j < d4est_schwarz_ops->max_degree; j++){
      d4est_schwarz_ops->schwarz_restrictor_1d_table[i][j] = NULL;
      d4est_schwarz_ops->schwarz_weights_1d_table[i][j] = NULL;
      d4est_schwarz_ops->schwarz_restrictor_transpose_1d_table[i][j] = NULL;
    }
  }
  return d4est_schwarz_ops;
}


void
d4est_solver_schwarz_operators_destroy
(
 d4est_solver_schwarz_operators_t* d4est_schwarz_ops
)
{
  d4est_schwarz_ops->d4est_ops = NULL;
  for (int i = 0; i < max_degree; i++){
    for(int j = 0; j < max_degree; j++){
      P4EST_FREE(d4est_schwarz_ops->schwarz_restrictor_1d_table[i][j]);
      P4EST_FREE(d4est_schwarz_ops->schwarz_weights_1d_table[i][j]);
      P4EST_FREE(d4est_schwarz_ops->schwarz_restrictor_transpose_1d_table[i][j]);
    }
    P4EST_FREE(d4est_schwarz_ops->schwarz_restrictor_1d_table[i]);
    P4EST_FREE(d4est_schwarz_ops->schwarz_weights_1d_table[i]);
    P4EST_FREE(d4est_schwarz_ops->schwarz_restrictor_transpose_1d_table[i]);
    P4EST_FREE(d4est_schwarz_ops->schwarz_restrictor_1d_table);
    P4EST_FREE(d4est_schwarz_ops->schwarz_weights_1d_table);
    P4EST_FREE(d4est_schwarz_ops->schwarz_restrictor_transpose_1d_table);
  }
}
