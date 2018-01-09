#include <pXest.h>
#include <d4est_lgl.h>
#include <math.h>
#include <d4est_util.h>
#include <d4est_kron.h>

int d4est_lgl_get_nodes(int dim, int deg) {
  D4EST_ASSERT(dim > 0 && dim < 4 && deg > 0);
  int volume_nodes = 1;
  for (int i = 0; i < dim; i++) volume_nodes *= deg + 1;
  return volume_nodes;
}

double d4est_lgl_jacobi(double r, double alpha, double beta, int N) {
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

double d4est_lgl_gradjacobi(double r, double alpha, double beta, int N) {
  if (N == 0)
    return 0;
  else
    return sqrt(N * (N + alpha + beta + 1.0)) *
           d4est_lgl_jacobi(r, alpha + 1.0, beta + 1.0, N - 1);
}

double
d4est_lgl_lagrange_1d(double x, double* lgl, int j, int deg){
  
  int N = deg;
  double l = 1.;
  for (int i = 0; i <= N; i++) {
    if (i != j)
      l *= (x - lgl[i])/(lgl[j] - lgl[i]);
  }
  return l;
}

/* double */
/* d4est_lgl_lagrange_2d(double x, double y, double* lgl, int j, int deg){ */
/* ind2sub2drowbeginat0[i_, cols_] := IntegerPart[i/cols] */
/* ind2sub2dcolbeginat0[i_, cols_] := IntegerPart[Mod[i, cols]] */
/* ind2sub2drowbeginat1[i_, cols_] := IntegerPart[1 + (i/cols)] */
/* ind2sub2dcolbeginat1[i_, cols_] := IntegerPart[1 + Mod[i, cols]] */

/* L2D[r_, s_, ind_, deg_] :=  */
/*  L[r, ind2sub2dcolbeginat1[ind - 1, deg + 1], deg]* */
/*   L[s, ind2sub2drowbeginat1[ind - 1, deg + 1], deg] */

/* Do[Print[StringJoin["mass[", ToString[i - 1], "][", ToString[j - 1],  */
/*     "] = ", ToString[ */
/*      NIntegrate[ */
/*       L2D[r, s, i, degree]*L2D[r, s, j, degree], {r, -1, 1}, {s, -1,  */
/*        1}]]]], {i, 1, (degree + 1)*(degree + 1)}, {j,  */
/*    1, (degree + 1)*(degree + 1)}]; */
/* } */



/* double */
/* d4est_lgl_lagrange_3d(double x, double y, double */
/* (*index starts at 1*) */

/* ind2sub3dbeginat1i[index_, nx_, ny_] :=  */
/*   IntegerPart[Mod[index, nx]] + 1; */
/* ind2sub3dbeginat1j[index_, nx_, ny_] :=  */
/*   IntegerPart[Mod[(index - Mod[index, nx])/nx, ny]] + 1; */
/* ind2sub3dbeginat1k[index_, nx_, ny_] :=  */
/*   IntegerPart[(((index - Mod[index, nx])/nx) -  */
/*        IntegerPart[Mod[(index - Mod[index, nx])/nx, ny]] + 1)/ny] + 1 ; */

/* (*ind must be 1 to (degree+1)^3 inclusive*) */

/* L3D[r_, s_, t_, ind_, deg_] :=  */
/*   L[r, ind2sub3dbeginat1i[ind - 1, deg + 1, deg + 1], deg]* */
/*    L[s, ind2sub3dbeginat1j[ind - 1, deg + 1, deg + 1], deg]* */
/*    L[t, ind2sub3dbeginat1k[ind - 1, deg + 1, deg + 1], deg]; */

/* Do[Print[StringJoin["mass[", ToString[i - 1], "][", ToString[j - 1],  */
/*     "] = ", ToString[ */
/*      NIntegrate[ */
/*       L3D[r, s, t, i, degree]*L3D[r, s, t, j, degree], {r, -1,  */
/*        1}, {s, -1, 1}, {t, -1, 1}]]]], {i,  */
/*    1, (degree + 1)*(degree + 1)*(degree + 1)}, {j,  */
/*    1, (degree + 1)*(degree + 1)*(degree + 1)}]; */
