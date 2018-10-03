#undef NDEBUG

#include <pXest.h>
#include <d4est_operators.h>
#include <d4est_linalg.h>
#include <d4est_kron.h>
#include <d4est_util.h>

#define D4EST_REAL_EPS 100*1e-15

static double
d4est_test_operators_poly_fcn
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
d4est_test_operators_interpolate
(
 d4est_operators_t* d4est_ops,
 int dim,
 int deg
)
{    
  int volume_nodes = d4est_lgl_get_nodes((dim), deg);

  double* rst_lobatto [] = {NULL, NULL, NULL};
  double* rst_gauss [] = {NULL, NULL, NULL};
  double* poly_lobatto = P4EST_ALLOC(double, volume_nodes);
 
  for (int i = 0; i < dim; i++){
    rst_lobatto[i] = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, dim, deg, i);
  }
  
  for (int i = 0; i < volume_nodes; i++){
    double xl = rst_lobatto[0][i];
    double yl = rst_lobatto[1][i];
    double zl = (dim == 3) ? rst_lobatto[2][i] : 0.;    
    poly_lobatto[i] = d4est_test_operators_poly_fcn(xl, yl, zl, dim, deg);
  }

  double rst_new [3];
  rst_new[0] = .5;
  rst_new[1] = 0.1;
  rst_new[2] = -0.1;
  double ans = d4est_test_operators_poly_fcn(rst_new[0], rst_new[1], rst_new[2], dim, deg);
  double ans_check = d4est_operators_interpolate(d4est_ops, &rst_new[0], poly_lobatto, dim, deg);

  P4EST_FREE(poly_lobatto);
  D4EST_ASSERT(fabs(ans -ans_check) < 1e-10);
}


void
d4est_test_operators_interp_lobatto_to_gauss
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
    
    poly_gauss[i] = d4est_test_operators_poly_fcn(xg, yg, zg, dim, deg_gauss);
    poly_lobatto[i] = d4est_test_operators_poly_fcn(xl, yl, zl, dim, deg_gauss);
  }
 
  double* interp = d4est_operators_fetch_lobatto_to_gauss_interp_1d(d4est_ops, deg_lobatto, deg_gauss);

  if (dim == 3){
    d4est_kron_A1A2A3x_nonsqr(poly_lobatto_to_gauss, interp, interp, interp, poly_lobatto,
                               nodes_gauss, nodes_lobatto, nodes_gauss, nodes_lobatto, nodes_gauss, nodes_lobatto);
  }
  else if (dim == 2){
    d4est_kron_A1A2x_nonsqr(poly_lobatto_to_gauss, interp, interp, poly_lobatto, nodes_gauss, nodes_lobatto,
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
    printf("d4est_test_operators_interp_lobatto_to_gauss passed for dim = %d, deg = %d\n", dim, deg);
  }

  P4EST_FREE(poly_lobatto);
  P4EST_FREE(poly_lobatto_to_gauss);
  P4EST_FREE(poly_gauss);
}

void
d4est_test_operators_mass_1d
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
    printf("d4est_test_operators_mass_1d passed\n");
  }
  P4EST_FREE(mij);
}

void
d4est_test_operators_mass_nd
(
 d4est_operators_t* d4est_ops,
 int dim,
 int deg
)
{
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


void d4est_test_operators_v1d(double* v1d, double* lobatto_nodes, int degree) {
  int i, j, rows, cols;
  rows = cols = degree + 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      v1d[i * cols + j] = d4est_lgl_jacobi(lobatto_nodes[i], 0., 0., j);
}

static
double
d4est_test_operators_lagrange_defn(double x, double* lgl, int j, int deg){
  
  int N = deg;
  double l = 1.;
  for (int i = 0; i <= N; i++) {
    if (i != j)
      l *= (x - lgl[i])/(lgl[j] - lgl[i]);
  }
  return l;
}

void
d4est_test_operators_p_prolong_1d
(
 d4est_operators_t* d4est_ops,
 int degH,
 int degh
)
{

  int nodesH = degH + 1;
  int nodesh = degh + 1;
  double* P = d4est_operators_fetch_p_prolong_1d(d4est_ops, degH, degh);
  double* lgl_h = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, degh);
  double* lgl_H = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, degH);

  for (int i = 0; i < nodesh; i++){
    for (int j = 0; j < nodesH; j++){
      double x = lgl_h[i];
      double p_ops = P[i*nodesH + j];
      double p_lag = d4est_test_operators_lagrange_defn(x, lgl_H, j, degH);
      double error = fabs(p_ops - p_lag);
      if (error > D4EST_REAL_EPS){
        printf("[D4EST_ERROR]: ERROR > D4EST_REAL_EPS -> p_ops, p_lag, error = %.15f, %.15f, %.15f\n", p_ops, p_lag, error);
        exit(1);
      }   
    }
  }
  printf("d4est_test_operators_p_prolong_1d degH = %d, degh = %d passed\n", degH, degh);
}

static void
d4est_test_operators_lagrange
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
    l_2[j] = d4est_test_operators_lagrange_defn(xp, x, j, deg);
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
    printf("d4est_test_operators_lagrange passed\n");
  }

  
  P4EST_FREE(l);
  P4EST_FREE(l_2);
}

void
d4est_test_operators_inv_vandermonde
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

  d4est_test_operators_v1d(invVij_2, x, deg);
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
    printf("d4est_test_operators_inv_vandermonde passed\n");
  }
  
  P4EST_FREE(invVij_2);
  P4EST_FREE(invVij);
}

static void
d4est_test_nodal_to_modal
(
 d4est_operators_t* d4est_ops,
 int dim,
 int deg
)
{
  int nodes = d4est_lgl_get_nodes(dim,deg + 1);
  double* r [3];
  for (int d = 0; d < dim; d++){
    r[d] = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, dim, deg, d);
  }
  /* double* vij = P4EST_ALLOC_ZERO(double, nodes*nodes); */
  double* u_nodal = P4EST_ALLOC_ZERO(double, nodes);
  double* u_nodal_2 = P4EST_ALLOC_ZERO(double, nodes);
  double* u_modal = P4EST_ALLOC_ZERO(double, nodes);
  double* u_modal_2 = P4EST_ALLOC_ZERO(double, nodes);
  double* error = P4EST_ALLOC_ZERO(double, nodes);
  double* error_modal = P4EST_ALLOC_ZERO(double, nodes);
  double* error_modal_check = P4EST_ALLOC_ZERO(double, nodes);


  for (int i = 0; i < deg + 1; i++){
    for (int d = 0; d < dim; d++){
      u_nodal[i] += pow(r[d][i],deg-1);
      u_nodal_2[i] += pow(r[d][i],deg-1) + .001*r[d][i] + .001;
    }
    error[i] = u_nodal[i] - u_nodal_2[i];
  }
  
  /* d4est_test_operators_v1d(vij, r, deg); */
  d4est_operators_convert_nodal_to_modal(d4est_ops, u_nodal, dim, deg, u_modal);
  d4est_operators_convert_nodal_to_modal(d4est_ops, u_nodal_2, dim, deg, u_modal_2);
  d4est_operators_convert_nodal_to_modal(d4est_ops, error, dim, deg, error_modal);

  for (int i = 0; i < nodes; i++){
    error_modal_check[i] = u_modal[i] - u_modal_2[i];
  }

  if (d4est_util_compare_vecs(error_modal, error_modal_check, nodes, 100*D4EST_REAL_EPS) == 0){
    printf("[D4EST_ERROR]: error_modal and errorl_modal_check are not equal\n");

    for (int j = 0; j < nodes; j++){
      printf("u_modal u_modal_2 error_modal = %.15f %.15f %.15f\n",  u_modal[j], u_modal_2[j], error_modal[j]);
    } 
    P4EST_FREE(u_nodal);
    P4EST_FREE(u_nodal_2);
    P4EST_FREE(u_modal);
    P4EST_FREE(u_modal_2);
    P4EST_FREE(error);
    P4EST_FREE(error_modal);
    P4EST_FREE(error_modal_check);
    exit(1);
  }
  else {
    printf("d4est_test_nodal_to_modal passed\n");
  }
  

  
  
  P4EST_FREE(u_nodal);
  P4EST_FREE(u_nodal_2);
  P4EST_FREE(u_modal);
  P4EST_FREE(u_modal_2);
  P4EST_FREE(error);
  P4EST_FREE(error_modal);
  P4EST_FREE(error_modal_check);
}



static void
d4est_test_operators_p_projection
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
      double li_xmort_j = d4est_test_operators_lagrange_defn(xmort_j, xH, i, deg_H);

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
 
  d4est_test_operators_v1d(invVij_h, xh, deg_h);
  d4est_test_operators_v1d(Vij_H, xH, deg_H);
  d4est_linalg_invert(invVij_h, nodes_h);
  for(int i = 0; i < nodes_H; i++)
    for (int j = 0; j < nodes_h; j++)
      invVij_h_trunc[i*nodes_h + j] = invVij_h[i*nodes_h + j];
  d4est_linalg_mat_multiply(Vij_H, invVij_h_trunc, Pij_3, nodes_H, nodes_H, nodes_h);
  /* d4est_util_print_matrices(Pij, Pij_2, nodes_H*nodes_h, 1, "Pij, Pij_2 = "); */
  /* d4est_util_print_matrices(Pij, Pij_3, nodes_H*nodes_h, 1, "Pij, Pij_3 = "); */
    
  if (d4est_util_compare_vecs(Pij, Pij_2, nodes_H*nodes_h, 100*D4EST_REAL_EPS) == 0){
    /* d4est_util_print_matrices(w,x, nodes, 1, "w,x = "); */
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
    printf("d4est_test_operators_p_projection passed\n");
  }
  
  P4EST_FREE(invVij_h_trunc);
  P4EST_FREE(Vij_H);
  P4EST_FREE(invVij_h);
  P4EST_FREE(Pij_3);  
  P4EST_FREE(Pij);
}

static void
d4est_test_operators_schwarz_nd
(
 d4est_operators_t* d4est_ops,
 int dim
)
{
  int res_deg = 1;
  int deg = 1;
  
  double* rst_lobatto [3];
  for (int i = 0; i < dim; i++){
    rst_lobatto[i] = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, dim, deg, i);  
  }

  int volume_nodes = d4est_lgl_get_nodes((dim), deg);
  double* frst = P4EST_ALLOC(double, volume_nodes);
  double* res_frst = P4EST_ALLOC(double, volume_nodes);
  double* res_frst_trans = P4EST_ALLOC(double, volume_nodes);
  double* res_rst_lobatto [3];

  for (int i = 0; i < volume_nodes; i++){
    frst[i] = 1.;
    for (int d = 0; d < (dim); d++){
      frst[i] *= exp(rst_lobatto[d][i]);
    }
  }

  int res_volume_nodes = 1;
  int faces [3];
  if (dim == 2){
    faces [0] = 1;
    faces [1] = 2;
    faces [2] = -1;
  }
  else {
    faces [0] = 0;
    faces [1] = 3;
    faces [2] = 4;
  }
  for (int i = 0; i < dim; i++){
    if (faces[i] != -1){
      res_volume_nodes *= (res_deg + 1);
    }
    else {
      res_volume_nodes *= (deg + 1);
    }
    res_rst_lobatto[i] = P4EST_ALLOC(double, volume_nodes);

    d4est_operators_apply_schwarz_restrictor
      (
       d4est_ops,
       rst_lobatto[i],
       dim,
       &faces[0],
       deg,
       res_deg+1,
       0,
       res_rst_lobatto[i]
      );
  }
  
  d4est_operators_apply_schwarz_restrictor
    (
     d4est_ops,
     frst,
     dim,
     &faces[0],
     deg,
     res_deg+1,
     D4OPS_NO_TRANSPOSE,
     res_frst
    );

  for (int i = 0; i < res_volume_nodes; i++){
    double exprst = 1.;
    for (int d = 0; d < (dim); d++){
      exprst *= exp(res_rst_lobatto[d][i]);
    }
    if (fabs(res_frst[i] - exprst) > 100*1e-14){
      D4EST_ABORT("fabs(res_frst[i] - exp(r)*exp(s)*exp(t)) > 100*1e-14");
    }
  }

  d4est_operators_apply_schwarz_restrictor
    (
     d4est_ops,
     res_frst,
     dim,
     &faces[0],
     deg,
     res_deg+1,
     D4OPS_TRANSPOSE,
     res_frst_trans
    );
  
  DEBUG_PRINT_ARR_DBL(frst, volume_nodes);
  DEBUG_PRINT_ARR_DBL(res_frst, res_volume_nodes);
  DEBUG_PRINT_ARR_DBL(res_frst_trans, volume_nodes);
  if (dim == 3){
    DEBUG_PRINT_3ARR_DBL(res_rst_lobatto[0],res_rst_lobatto[1],res_rst_lobatto[2], res_volume_nodes);
  }
  else {
    DEBUG_PRINT_2ARR_DBL(res_rst_lobatto[0],res_rst_lobatto[1], res_volume_nodes);
  }
  P4EST_FREE(frst);
  P4EST_FREE(res_frst);
  P4EST_FREE(res_frst_trans);

  for (int i = 0; i < dim; i++){
    P4EST_FREE(res_rst_lobatto[i]);
  }
}

static void
d4est_test_operators_schwarz
(
 d4est_operators_t* d4est_ops
){


  int deg = 3;
  int res_deg = 3;

  int nodes_res = res_deg + 1;
  int nodes_deg = deg + 1;
  
  double* x = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, deg);
  double* fx = P4EST_ALLOC(double, deg + 1);

  double* res_x_side0 = P4EST_ALLOC(double, res_deg + 1);
  double* res_fx_side0 = P4EST_ALLOC(double, res_deg + 1);

  double* res_x_side1 = P4EST_ALLOC(double, res_deg + 1);
  double* res_fx_side1 = P4EST_ALLOC(double, res_deg + 1);

  double* res_x_side0_trans = P4EST_ALLOC(double, deg + 1);
  double* res_fx_side0_trans = P4EST_ALLOC(double, deg + 1);
  double* res_x_side1_trans = P4EST_ALLOC(double, deg + 1);
  double* res_fx_side1_trans = P4EST_ALLOC(double, deg + 1);

  for (int i = 0; i < deg + 1; i++){
    fx[i] = exp(x[i]);
  }
  
  double* schwarz_1d = d4est_operators_fetch_schwarz_restrictor_1d
                       (
                        d4est_ops,
                        deg,
                        res_deg
                       );
  
  double* schwarz_1d_trans = d4est_operators_fetch_schwarz_restrictor_transpose_1d
                             (
                              d4est_ops,
                              deg,
                              res_deg
                             );


  
  d4est_linalg_matvec_plus_vec(1.0, schwarz_1d, x, 0., res_x_side0, nodes_res, nodes_deg);       
  d4est_linalg_matvec_plus_vec(1.0, schwarz_1d, fx, 0., res_fx_side0, nodes_res, nodes_deg);
  d4est_linalg_matvec_plus_vec(1.0, &schwarz_1d[nodes_deg*nodes_res], x, 0., res_x_side1, nodes_res, nodes_deg);       
  d4est_linalg_matvec_plus_vec(1.0, &schwarz_1d[nodes_deg*nodes_res], fx, 0., res_fx_side1, nodes_res, nodes_deg);

  d4est_linalg_matvec_plus_vec(1.0, schwarz_1d_trans, res_x_side0, 0., res_x_side0_trans, nodes_deg, nodes_res);       
  d4est_linalg_matvec_plus_vec(1.0, schwarz_1d_trans, res_fx_side0, 0., res_fx_side0_trans, nodes_deg, nodes_res);

  
  d4est_linalg_matvec_plus_vec(1.0, &schwarz_1d_trans[nodes_deg*nodes_res], res_x_side1, 0., res_x_side1_trans, nodes_deg, nodes_res);       
  d4est_linalg_matvec_plus_vec(1.0, &schwarz_1d_trans[nodes_deg*nodes_res], res_fx_side1, 0., res_fx_side1_trans, nodes_deg, nodes_res);

  /* d4est_util_print_matrix(&schwarz_1d[0], */
  /*                         res_deg+1, deg + 1, "schwarz_1d_side0 = ", 0); */
  
  /* d4est_util_print_matrix(&schwarz_1d[(deg + 1)*(res_deg+1)], */
  /*                         res_deg+1, deg + 1, "schwarz_1d_side1 = ", 0); */


  /* d4est_util_print_matrix(&schwarz_1d_trans[0], */
  /*                         deg+1, res_deg + 1, "schwarz_1d_trans_side0 = ", 0); */
  
  /* d4est_util_print_matrix(&schwarz_1d_trans[(deg + 1)*(res_deg+1)], */
  /*                         deg+1, res_deg + 1, "schwarz_1d_trans_side1 = ", 0); */
  
  /* DEBUG_PRINT_2ARR_DBL(x,fx,deg+1); */
  /* DEBUG_PRINT_2ARR_DBL(res_x_side0,res_fx_side0,res_deg+1); */
  /* DEBUG_PRINT_2ARR_DBL(res_x_side1,res_fx_side1,res_deg+1); */
  /* DEBUG_PRINT_2ARR_DBL(res_x_side0_trans,res_fx_side0_trans,deg+1); */
  /* DEBUG_PRINT_2ARR_DBL(res_x_side1_trans,res_fx_side1_trans,deg+1); */

  for (int i = 0; i < res_deg + 1; i++){
    D4EST_ASSERT(fabs(exp(res_x_side0[i]) - res_fx_side0[i]) < 1e-14);
    D4EST_ASSERT(fabs(exp(res_x_side1[i]) - res_fx_side1[i]) < 1e-14);
    D4EST_ASSERT(res_fx_side0_trans[i] == 0. || fabs(exp(res_x_side0_trans[i]) - res_fx_side0_trans[i]) < 1e-14);
    D4EST_ASSERT(res_fx_side1_trans[i] == 0. || fabs(exp(res_x_side1_trans[i]) - res_fx_side1_trans[i]) < 1e-14);
  }
  
  P4EST_FREE(fx);
  P4EST_FREE(res_x_side0);
  P4EST_FREE(res_fx_side0);
  P4EST_FREE(res_x_side1);
  P4EST_FREE(res_fx_side1);

  P4EST_FREE(res_x_side0_trans);
  P4EST_FREE(res_fx_side0_trans);
  P4EST_FREE(res_x_side1_trans);
  P4EST_FREE(res_fx_side1_trans);
  
}

int main(int argc, char *argv[])
{
  d4est_operators_t* d4est_ops = d4est_ops_init(20);  

  d4est_test_nodal_to_modal
    (
     d4est_ops,
     3,
     3
    );
  
  /* printf("mass 2d = \n"); */
  d4est_test_operators_mass_nd
    (
     d4est_ops,
     2,
     1
    );
  
  /* printf("mass 3d = \n"); */

  d4est_test_operators_mass_nd
    (
     d4est_ops,
     3,
     1
    );

  /* d4est_test_operators_schwarz_nd */
  /*   ( */
  /*    d4est_ops, */
  /*    3 */
  /*   ); */

  /* d4est_test_operators_schwarz */
  /*   ( */
  /*    d4est_ops */
  /*   ); */
  
  d4est_ops_destroy(d4est_ops);
  return 0;
}


