#undef NDEBUG

#include <pXest.h>
#include <d4est_operators.h>
#include <d4est_linalg.h>
#include <d4est_util.h>

#define D4EST_REAL_EPS 1e-15

static double
test_d4est_operators_poly_fcn
(
 double r,
 double s,
 double t,
 int dim,
 int deg
){
  double poly = pow(r,deg) + pow(s,deg);
  poly += (dim == 3) ? pow(t,deg) : 0;
  return poly;
}

void
test_d4est_operators_interp_lobatto_to_gauss
(
 d4est_operators_t* d4est_ops,
 int dim,
 int deg
)
{
  int deg_lobatto = deg;
  int deg_gauss = deg;

  int volume_nodes_lobatto = d4est_lgl_get_nodes((dim), deg_lobatto);
  int volume_nodes_gauss = d4est_lgl_get_nodes((dim), deg_gauss);

  int nodes_lobatto = deg_lobatto + 1;
  int nodes_gauss = deg_gauss + 1;
    
  double* rst_lobatto [] = {NULL, NULL, NULL};
  double* rst_gauss [] = {NULL, NULL, NULL};
  double* poly_lobatto = P4EST_ALLOC(double, volume_nodes_lobatto);
  double* poly_lobatto_to_gauss = P4EST_ALLOC(double, volume_nodes_gauss);
  double* poly_gauss = P4EST_ALLOC(double, volume_nodes_gauss);
 
  for (int i = 0; i < dim; i++){
    rst_lobatto[i] = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, dim, deg, i);
    rst_gauss[i] = d4est_operators_fetch_gauss_rst_nd(d4est_ops, dim, deg, i);
  }
  
  for (int i = 0; i < volume_nodes_gauss; i++){

    double xg = rst_gauss[0][i];
    double yg = rst_gauss[1][i];
    double zg = (dim == 3) ? rst_gauss[2][i] : 0.;

    double xl = rst_lobatto[0][i];
    double yl = rst_lobatto[1][i];
    double zl = (dim == 3) ? rst_lobatto[2][i] : 0.;
    
    poly_gauss[i] = test_d4est_operators_poly_fcn(xg, yg, zg, dim, deg_gauss);
    poly_lobatto[i] = test_d4est_operators_poly_fcn(xl, yl, zl, dim, deg_gauss);
  }
 
  double* interp = d4est_operators_fetch_lobatto_to_gauss_interp_1d(d4est_ops, deg_lobatto, deg_gauss);

  if (dim == 3){
    d4est_linalg_kron_A1A2A3x_nonsqr(poly_lobatto_to_gauss, interp, interp, interp, poly_lobatto,
                               nodes_gauss, nodes_lobatto, nodes_gauss, nodes_lobatto, nodes_gauss, nodes_lobatto);
  }
  else if (dim == 2){
    d4est_linalg_kron_A1A2x_nonsqr(poly_lobatto_to_gauss, interp, interp, poly_lobatto, nodes_gauss, nodes_lobatto,
                             nodes_gauss, nodes_lobatto);
  }
  else if (dim == 1){
    d4est_linalg_matvec_plus_vec(1.0,interp, poly_lobatto, 0., poly_lobatto_to_gauss, nodes_gauss, nodes_lobatto);
  }
  else {
    D4EST_ABORT("[D4EST_ERROR]: dim = 1, 2 or 3\n");
  }  

  int error = 0;
  /* double max_error = 0.; */
  for (int i = 0; i < volume_nodes_gauss; i++){
    error += !(d4est_util_compare_double(poly_lobatto_to_gauss[i],poly_gauss[i],D4EST_REAL_EPS));
    /* error_tmp = fabs(poly_lobatto_to_gauss[i] - poly_gauss[i]); */
    /* max_error = ; */
  }
 
  if (error != 0){
    printf("[D4EST_ERROR]: poly_lobatto_to_gauss[i] - poly_gauss[i] not equal\n");
    DEBUG_PRINT_2ARR_DBL(poly_gauss, poly_lobatto_to_gauss, volume_nodes_gauss);
    exit(1);
  }
  else {
    printf("test_d4est_operators_interp_lobatto_to_gauss passed for dim = %d, deg = %d\n", dim, deg);
  }

  P4EST_FREE(poly_lobatto);
  P4EST_FREE(poly_lobatto_to_gauss);
  P4EST_FREE(poly_gauss);
  
}

void
test_d4est_operators_mass_1d
(
 d4est_operators_t* d4est_ops,
 int deg
)
{
  int nodes = deg + 1;
  int N = nodes - 1;
  
  /* matrix of size nodes_H x nodes_h */
  double* mij = P4EST_ALLOC(double, nodes*nodes);

  double* w = d4est_operators_fetch_lobatto_weights_1d(d4est_ops, deg);
  double* x = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, deg);

  /* Compute mass from Teukolsky mass matrix paper */
  for (int i = 0; i < nodes; i++){
    for (int j = 0; j < nodes; j++){
      double hN = 2./(2.*(double)N + 1.);
      double gN = 2./(double)N;
      double alpha = (hN - gN)/(gN*gN);
      double M = 0.;
      double pNxi = d4est_lgl_jacobi(x[i], 0., 0., N);
      double pNxj = d4est_lgl_jacobi(x[j], 0., 0., N);
      if (i == j)
        M += w[i];
      M += alpha*w[i]*w[j]*pNxi*pNxj*hN;
      mij[i*nodes + j] = M;
    }
  }
  double* mij_2 = d4est_operators_fetch_mij_1d(d4est_ops, deg);
  if (d4est_util_compare_vecs(mij, mij_2, nodes*nodes, D4EST_REAL_EPS) == 0){
    d4est_util_print_matrices(w,x, nodes, 1, "w,x = ");
    d4est_util_print_matrices(mij, mij_2, nodes*nodes, 1, "mij, mij_2 = ");
    P4EST_FREE(mij);
    printf("[D4EST_ERROR]:mij and mij_2 not equal\n");
    exit(1);
  }
  else {
    printf("test_d4est_operators_mass_1d passed\n");
  }
  P4EST_FREE(mij);
}

void
test_d4est_operators_mass_nd
(
 d4est_operators_t* d4est_ops,
 int dim
)
{
  int deg = 2;
  int volume_nodes = d4est_lgl_get_nodes(dim, deg);

  double* mij = P4EST_ALLOC_ZERO(double, volume_nodes*volume_nodes);
  double* u = P4EST_ALLOC_ZERO(double, volume_nodes);
  double* Mu = P4EST_ALLOC_ZERO(double, volume_nodes);

  for (int i = 0; i < volume_nodes; i++){
    u[i] = 1.;
    d4est_operators_apply_mij(d4est_ops, u, dim, deg, Mu);
    d4est_linalg_set_column(mij, Mu, i, volume_nodes, volume_nodes);
    u[i] = 0.;
  }

  for (int i = 0; i < volume_nodes; i++){
    for (int j = 0; j < volume_nodes; j++){
      printf("M[%d][%d] = %.15f\n", i,j,mij[i*volume_nodes + j]);
    }
  }

  P4EST_FREE(Mu);
  P4EST_FREE(mij);
  P4EST_FREE(u);
}


void test_d4est_operators_v1d(double* v1d, double* lobatto_nodes, int degree) {
  int i, j, rows, cols;
  rows = cols = degree + 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      v1d[i * cols + j] = d4est_lgl_jacobi(lobatto_nodes[i], 0., 0., j);
}

static
double
test_d4est_operators_lagrange_defn(double x, double* lgl, int j, int deg){
  
  int N = deg;
  double l = 1.;
  for (int i = 0; i <= N; i++) {
    if (i != j)
      l *= (x - lgl[i])/(lgl[j] - lgl[i]);
  }
  return l;
}

static void
test_d4est_operators_lagrange
(
 d4est_operators_t* d4est_ops
)
{
  /* int degH = 3; */
  int deg = 5;
  /* int nodes_H = degH + 1; */
  int nodes = deg + 1;
  /* int N = nodes - 1; */
  /* P4EST_ASSERT(degh >= degH); */
  double xp = .5;
  
  /* matrix of size nodes_H x nodes_h */
  double* l = P4EST_ALLOC(double, nodes);
  double* l_2 = P4EST_ALLOC(double, nodes);

  /* double* wH = d4est_operators_fetch_LGL_weights_1d(d4est_ops, degH); */
  double* w = d4est_operators_fetch_lobatto_weights_1d(d4est_ops, deg);
  double* x = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, deg);
  
  for (int j = 0; j < nodes; j++){
    l[j] = 0.;
    for (int k = 0; k < nodes; k++){
      double hk = (2./(2.*(double)k + 1.));  
      double gk = (k == nodes - 1) ? (2./(double)k) : hk;
      double pkxj = d4est_lgl_jacobi(x[j], 0., 0., k);
      double pkx = d4est_lgl_jacobi(xp, 0., 0., k);
      l[j] += (hk/gk)*pkxj*pkx;
    }
    l[j] *= w[j];
    l_2[j] = test_d4est_operators_lagrange_defn(xp, x, j, deg);
  }

  if (d4est_util_compare_vecs(l, l_2, nodes, 100*D4EST_REAL_EPS) == 0){
    d4est_util_print_matrices(w,x, nodes, 1, "w,x = ");
    d4est_util_print_matrices(l, l_2, nodes, 1, "l, l_2 = ");
    P4EST_FREE(l);
    P4EST_FREE(l_2);
    printf("[D4EST_ERROR]:l and l_2 not equal\n");
    exit(1);
  }
  else {
    printf("test_d4est_operators_lagrange passed\n");
  }

  
  P4EST_FREE(l);
  P4EST_FREE(l_2);
}

void
test_d4est_operators_inv_vandermonde
(
 d4est_operators_t* d4est_ops
)
{
  /* int degH = 3; */
  int deg = 5;
  /* int nodes_H = degH + 1; */
  int nodes = deg + 1;
  /* int N = nodes - 1; */
  /* P4EST_ASSERT(degh >= degH); */
  
  /* matrix of size nodes_H x nodes_h */
  double* invVij = P4EST_ALLOC(double, nodes*nodes);
  double* invVij_2 = P4EST_ALLOC(double, nodes*nodes);

  /* double* wH = d4est_operators_fetch_lobatto_weights_1d(d4est_ops, degH); */
  double* w = d4est_operators_fetch_lobatto_weights_1d(d4est_ops, deg);
  double* x = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, deg);

  for (int k = 0; k < nodes; k++){
    double hk = (2./(2.*(double)k + 1.));  
    double gk = (k == nodes - 1) ? (2./(double)k) : hk;
    for (int j = 0; j < nodes; j++){
      double hi = (2./(2.*(double)k + 1.));
      double pkxj = d4est_lgl_jacobi(x[j], 0., 0., k)*sqrt(hk)*sqrt(hi);
      invVij[k*nodes + j] = (1./gk)*w[j]*pkxj;
    }
  }

  test_d4est_operators_v1d(invVij_2, x, deg);
  d4est_linalg_invert(invVij_2, nodes);

  if (d4est_util_compare_vecs(invVij, invVij_2, nodes*nodes, 100*D4EST_REAL_EPS) == 0){
    d4est_util_print_matrices(w,x, nodes, 1, "w,x = ");
    d4est_util_print_matrices(invVij, invVij_2, nodes*nodes, 1, "invVij, invVij_2 = ");
    P4EST_FREE(invVij_2);
    P4EST_FREE(invVij);
    printf("[D4EST_ERROR]:invVij and invVij_2 not equal\n");
    exit(1);
  }
  else {
    printf("test_d4est_operators_inv_vandermonde passed\n");
  }
  
  P4EST_FREE(invVij_2);
  P4EST_FREE(invVij);
}


static void
test_d4est_operators_p_projection
(
 d4est_operators_t* d4est_ops
)
{
  int deg_H = 3;
  int deg_h = 5;
  int nodes_H = deg_H + 1;
  int nodes_h = deg_h + 1;
  P4EST_ASSERT(deg_h >= deg_H);
  
  /* matrix of size nodes_H x nodes_h */
  double* Pij = P4EST_ALLOC(double, nodes_H*nodes_h);

  double* wH = d4est_operators_fetch_lobatto_weights_1d(d4est_ops, deg_H);
  double* wh = d4est_operators_fetch_lobatto_weights_1d(d4est_ops, deg_h);

  double* xH = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, deg_H);
  double* xh = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, deg_h);

  int Ns = nodes_H - 1;

  for (int i = 0; i < nodes_H; i++){
    for (int j = 0; j < nodes_h; j++){
      double wmort_j = wh[j];
      double w_i = wH[i];
      double xmort_j = xh[j];
      double li_xmort_j = test_d4est_operators_lagrange_defn(xmort_j, xH, i, deg_H);

      double hNs = 2./(2.*(double)Ns + 1.);
      double gNs = 2./(double)Ns;
      double beta = -(hNs - gNs)/(gNs*hNs);
      
      double p_Ns_xi = d4est_lgl_jacobi(xH[i], 0., 0., Ns);
      double p_Ns_xmort_j = d4est_lgl_jacobi(xmort_j, 0., 0., Ns);

      Pij[i*nodes_h + j] = (wmort_j/w_i)*li_xmort_j + beta*p_Ns_xi*p_Ns_xmort_j*wmort_j*hNs;
    }
  }

  double* Pij_2 = d4est_operators_fetch_p_restrict_1d(d4est_ops, deg_H, deg_h);  
  double* Pij_3 = P4EST_ALLOC(double, nodes_H*nodes_h);
  double* invVij_h = P4EST_ALLOC(double, nodes_h*nodes_h);
  double* Vij_H = P4EST_ALLOC(double, nodes_H*nodes_H);
  double* invVij_h_trunc = P4EST_ALLOC(double, nodes_h*nodes_H);
 
  test_d4est_operators_v1d(invVij_h, xh, deg_h);
  test_d4est_operators_v1d(Vij_H, xH, deg_H);
  d4est_linalg_invert(invVij_h, nodes_h);
  for(int i = 0; i < nodes_H; i++)
    for (int j = 0; j < nodes_h; j++)
      invVij_h_trunc[i*nodes_h + j] = invVij_h[i*nodes_h + j];
  d4est_linalg_mat_multiply(Vij_H, invVij_h_trunc, Pij_3, nodes_H, nodes_H, nodes_h);
  /* d4est_util_print_matrices(Pij, Pij_2, nodes_H*nodes_h, 1, "Pij, Pij_2 = "); */
  /* d4est_util_print_matrices(Pij, Pij_3, nodes_H*nodes_h, 1, "Pij, Pij_3 = "); */
    
  if (d4est_util_compare_vecs(Pij, Pij_2, nodes_H*nodes_h, 100*D4EST_REAL_EPS) == 0){
    /* d4est_d4est_util_print_matrices(w,x, nodes, 1, "w,x = "); */
    d4est_util_print_matrices(Pij, Pij_2, nodes_H*nodes_h, 1, "Pij, Pij_2 = ");
    d4est_util_print_matrices(Pij, Pij_3, nodes_H*nodes_h, 1, "Pij, Pij_3 = ");
    P4EST_FREE(Pij);
    P4EST_FREE(invVij_h_trunc);
    P4EST_FREE(Vij_H);
    P4EST_FREE(invVij_h);
    P4EST_FREE(Pij_3);  
    /* P4EST_FREE(Pij_2); */
    printf("[D4EST_ERROR]:Pij and Pij_2 not equal\n");
    exit(1);
  }
  else {
    printf("test_d4est_operators_p_projection passed\n");
  }
  
  P4EST_FREE(invVij_h_trunc);
  P4EST_FREE(Vij_H);
  P4EST_FREE(invVij_h);
  P4EST_FREE(Pij_3);  
  P4EST_FREE(Pij);
}



int main(int argc, char *argv[])
{
  d4est_operators_t* d4est_ops = d4est_ops_init(20);  
  test_d4est_operators_mass_1d(d4est_ops, 3);
  test_d4est_operators_inv_vandermonde(d4est_ops);
  test_d4est_operators_p_projection(d4est_ops);
  test_d4est_operators_lagrange(d4est_ops);
  for (int dim = 3; dim < 4; dim++)
    for (int deg = 1; deg < 20; deg++)
      test_d4est_operators_interp_lobatto_to_gauss(d4est_ops,dim,deg);
  d4est_ops_destroy(d4est_ops);
  return 0;
}
