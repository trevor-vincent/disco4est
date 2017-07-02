#define _GNU_SOURCE
/* #define TEST_SYM  /\* Test if matrix is symmetric *\/ */
/* #define NDEBUG /\* TURN DEBUG OFF/ON *\/ */

/* TURN DEBUG OFF/ON */
#define NDEBUG

#include <sc_reduce.h>
#include <pXest.h>
#include <util.h>
#include <d4est_linalg.h>
#include <element_data.h>
#include <sipg_flux_vector_fcns.h>
#include <sipg_flux_scalar_fcns.h>
#include <problem.h>
#include <problem_data.h>
#include <d4est_elliptic_eqns.h>
#include <poisson_operator.h>
#include <hp_amr_smooth_pred.h>
#include <hp_amr_uniform.h>
#include <bi_estimator.h>
#include <bi_estimator_flux_fcns.h>
#include <hp_amr.h>
#include <ini.h>
#include <newton_petsc.h>
#include <hacked_p4est_vtk.h>
#include <dg_norm.h>
#include <grad.h>
#include "time.h"

#define NUM_PUNCTURES 2
static const double pi = 3.1415926535897932384626433832795;
static const double puncture_eps = 0.00000000001;// .000000000001;
static double xyz_bh [NUM_PUNCTURES][3];
static double P_bh [NUM_PUNCTURES][3];
static double S_bh [NUM_PUNCTURES][3];
static double M_bh [NUM_PUNCTURES];

static
double levi_civita(int a, int b, int c)
{
  double eps;
  if( ( ( a == 0 )&&( b == 1 )&&( c == 2 ) ) ||
      ( ( a == 1 )&&( b == 2 )&&( c == 0 ) ) ||
      ( ( a == 2 )&&( b == 0 )&&( c == 1 ) ) ) {
    eps = 1.;
  } else
    if( ( ( a == 1 )&&( b == 0 )&&( c == 2 ) ) ||
        ( ( a == 0 )&&( b == 2 )&&( c == 1 ) ) ||
        ( ( a == 2 )&&( b == 1 )&&( c == 0 ) ) ) {
      eps = -1.;
    } else {
      eps = 0.;
    }
  return eps;
}

typedef struct {

  int deg_offset_for_gauss_quad;
  
} problem_ctx_t;

static
double kronecker(int a, int b)
{
  /* printf("a = %d\n", a); */
  /* printf("b = %d\n", b); */
  /* printf("a == b = %d\n", (a==b)); */
  
  if (a != b)
    return 0.0;
  else
    return 1.0;

  /* return (double)(a==b); */
}


static
double boundary_fcn
(
 double x,
 double y,
 double z
)
{
  return zero_fcn(x,y,z);
}

static
double psi_fcn
(
 double x,
 double y,
 double z,
 double u
)
{
  double sumn_mn_o_2rn = 0.;
  double dxn, dyn, dzn, r;
  int n;
  for (n = 0; n < NUM_PUNCTURES; n++){
    dxn = x - xyz_bh[n][0];
    dyn = y - xyz_bh[n][1];
    dzn = z - xyz_bh[n][2];
    r = sqrt(dxn*dxn + dyn*dyn + dzn*dzn);
    if (r == 0.)
      r += puncture_eps;
    sumn_mn_o_2rn += M_bh[n]/(2.*r);
  }

  return 1. + u + sumn_mn_o_2rn;
}


double two_punctures_adm_quadral_face_contribution
(
 element_data_t* elem_data,
 int face,
 void* ctx,
 d4est_operators_t* d4est_ops
)
{
  double* u = (double*)ctx;
  int face_nodes = d4est_lgl_get_nodes((P4EST_DIM)-1, elem_data->deg);
  int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), elem_data->deg);
  double* restrict u_elem = &u[elem_data->stride];

  
  double* psi_elem = P4EST_ALLOC(double, volume_nodes);

  double h = elem_data->h;
  double* x = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM), elem_data->deg, 0);
  double xl = elem_data->xyz_corner[0];
  double* y = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM), elem_data->deg, 1);
  double yl = elem_data->xyz_corner[1];  
  double* z = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM), elem_data->deg, 2);
  double zl = elem_data->xyz_corner[2];

  int i;
  for (i = 0; i < volume_nodes; i++){
    double xp = d4est_reference_rtox(x[i], xl, h);
    double yp = d4est_reference_rtox(y[i], yl, h);
    double zp = d4est_reference_rtox(z[i], zl, h);
    psi_elem[i] = psi_fcn(xp,yp,zp,u_elem[i]);
  }
                                 
  double* grad_psi_elem [(P4EST_DIM)];

  for (int d = 0; d < (P4EST_DIM); d++){
    grad_psi_elem[d] = P4EST_ALLOC(double, volume_nodes);
  }
  
  grad(psi_elem, grad_psi_elem, elem_data->h, elem_data->deg, d4est_ops);
  double n [(P4EST_DIM)];
  d4est_reference_get_normal(face, (P4EST_DIM), &n[0]);
  
  double* grad_psi_face = P4EST_ALLOC(double, face_nodes);
  double* grad_psi_face_dot_n = P4EST_ALLOC_ZERO(double, face_nodes);
  double* M_grad_psi_face_dot_n = P4EST_ALLOC(double, face_nodes);
  
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_slicer(d4est_ops, &grad_psi_elem[d][0], (P4EST_DIM), face, elem_data->deg, grad_psi_face);
    for (int fn = 0; fn < face_nodes; fn++){
      grad_psi_face_dot_n[fn] += grad_psi_face[fn]*n[d];
    }
  }

  d4est_operators_apply_mij(d4est_ops, grad_psi_face_dot_n, (P4EST_DIM)-1, elem_data->deg, M_grad_psi_face_dot_n);
  d4est_linalg_vec_scale(elem_data->surface_jacobian, M_grad_psi_face_dot_n, face_nodes);
  
  double face_contrib = 0.;
  for (int fn = 0; fn < face_nodes; fn++){
    face_contrib += M_grad_psi_face_dot_n[fn];
  }
   
  face_contrib *= -(1./(2.*M_PI));

  for (int d = 0; d < (P4EST_DIM); d++){
    P4EST_FREE(grad_psi_elem[d]);
  }
  P4EST_FREE(grad_psi_face);
  P4EST_FREE(grad_psi_face_dot_n);
  P4EST_FREE(M_grad_psi_face_dot_n);
  P4EST_FREE(psi_elem);
  return face_contrib;

}


double two_punctures_adm_quadral_volume
(
 double* u,
 double m1,
 double m2,
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_elliptic_problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops
)
{
  double* tmp = prob_vecs->u;

  poisson_apply_aij
    (
     p4est,
     ghost,
     ghost_data,
     prob_vecs,
     d4est_ops
    );

  double* M_Au = P4EST_ALLOC(double, prob_vecs->local_nodes);

  element_data_apply_mij_on_vec(p4est, prob_vecs->Au, M_Au, d4est_ops);
 
  double adm_quadral = 0.;
  for (int i = 0; i < prob_vecs->local_nodes; i++){
    adm_quadral += M_Au[i];
  }

  /* don't need negative here because Au is -\nabla^2 */
  adm_quadral *= (1./(2.*M_PI));
  adm_quadral += m1 + m2;
  prob_vecs->u = tmp;
  P4EST_FREE(M_Au);
  return adm_quadral;
}

static
double compute_confAij
(
 int i,
 int j,
 double x,
 double y,
 double z
)
{
  double nvec[NUM_PUNCTURES][3];
  double confAij = 0.;
  
  /* initialize and check if (x,y,z) is a puncture point */
  int n;
  for (n = 0; n < NUM_PUNCTURES; n++){
    double dxn = (x - xyz_bh[n][0]);
    double dyn = (y - xyz_bh[n][1]);
    double dzn = (z - xyz_bh[n][2]);
    double rn = sqrt(dxn*dxn + dyn*dyn + dzn*dzn);
    if (rn == 0.)
      rn += puncture_eps;
    nvec[n][0] = (x - xyz_bh[n][0])/rn;
    nvec[n][1] = (y - xyz_bh[n][1])/rn;
    nvec[n][2] = (z - xyz_bh[n][2])/rn;

    double lcikl_dot_Sk_dot_nln_times_njn = 0.;
    double lcjkl_dot_Sk_dot_nln_times_nin = 0.;
    double Pn_dot_nvec = 0.;
    
    int k,l;
    for (k = 0; k < 3; k++){
      Pn_dot_nvec += P_bh[n][k]*nvec[n][k];
      for (l = 0; l < 3; l++){
        lcikl_dot_Sk_dot_nln_times_njn += (levi_civita(i,k,l))*S_bh[n][k]*nvec[n][l]*nvec[n][j];
        lcjkl_dot_Sk_dot_nln_times_nin += (levi_civita(j,k,l))*S_bh[n][k]*nvec[n][l]*nvec[n][i];
      }
    }

    double a = nvec[n][i]*P_bh[n][j] + nvec[n][j]*P_bh[n][i] - (kronecker(i,j) - nvec[n][i]*nvec[n][j])*Pn_dot_nvec;
    double b = lcikl_dot_Sk_dot_nln_times_njn + lcjkl_dot_Sk_dot_nln_times_nin;

    confAij += (1./(rn*rn))*(a + (2./rn)*b);
  }

  confAij *= 3./2.;
  return confAij;
}

static
double compute_confAij_sqr
(
 double x,
 double y,
 double z
)
{
  double confAij;
  double confAij_sqr = 0.;
  int i,j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++){
      confAij = compute_confAij(i,j,x,y,z);
      confAij_sqr += confAij*confAij;
    }
  return confAij_sqr;
}

static
double neg_1o8_K2_psi_neg7
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  double sumn_mn_o_2rn = 0.;
  double dxn, dyn, dzn, r;
  int n;
  for (n = 0; n < NUM_PUNCTURES; n++){
    dxn = x - xyz_bh[n][0];
    dyn = y - xyz_bh[n][1];
    dzn = z - xyz_bh[n][2];
    r = sqrt(dxn*dxn + dyn*dyn + dzn*dzn);
    sumn_mn_o_2rn += M_bh[n]/(2.*r);
  }

  double psi_0 = 1. + u + sumn_mn_o_2rn;
  double confAij_sqr = compute_confAij_sqr(x,y,z);

  if (r > puncture_eps)
    return (-1./8.)*confAij_sqr/(psi_0*psi_0*psi_0*psi_0*psi_0*psi_0*psi_0);
  else
    return 0.;
  /* return (1000./8.)*confAij_sqr/(psi_0*psi_0*psi_0*psi_0*psi_0*psi_0*psi_0); */
}



static
double plus_7o8_K2_psi_neg8
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  double sumn_mn_o_2rn = 0.;
  double dxn, dyn, dzn, r;
  int n;
  for (n = 0; n < NUM_PUNCTURES; n++){
    dxn = x - xyz_bh[n][0];
    dyn = y - xyz_bh[n][1];
    dzn = z - xyz_bh[n][2];
    r = sqrt(dxn*dxn + dyn*dyn + dzn*dzn);
    sumn_mn_o_2rn += M_bh[n]/(2.*r);
  }

  double psi_0 = 1. + u + sumn_mn_o_2rn;
  double confAij_sqr = compute_confAij_sqr(x,y,z);

  if (r > puncture_eps)
    return (7./8.)*confAij_sqr/(psi_0*psi_0*psi_0*psi_0*psi_0*psi_0*psi_0*psi_0);
  else
    return 0.;
}

/* static */
/* double f_fcn */
/* ( */
/*  double x, */
/*  double y */
/* #if (P4EST_DIM)==3 */
/*  , */
/*  double z */
/* #endif */
/* ) */
/* { */
/*   return 0.; */
/* } */

/** 
 * ASSUMES F(u0) = f - u0'' - u0^2
 * with Mass*f already stored in rhs
 *
 * @param p4est 
 * @param ghost 
 * @param ghost_data 
 * @param prob_vecs 
 */
static
void
build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_elliptic_problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops
)
{
  /* printf("HI FROM RESIDUAL\n"); */
  /* printf("xbh,P,S,M_bh = %f,%f,%f,%f\n",x_bh,P_bh,S_bh,M_bh); */
  /* -u'' */
  
  poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops);  

  /* printf("OLD VERSION\n"); */
  /* /\* integral(f(u)) =+ Au *\/ */
  /* element_data_quadrate_fofuv_andaddto(p4est, */
  /*                                       prob_vecs->u, */
  /*                                       NULL, */
  /*                                       prob_vecs->Au, */
  /*                                       /\* -1, *\/ */
  /*                                       neg_1o8_K2_psi_neg7, */
  /*                                       d4est_ops */
  /*                                      ); */


  /* printf("NEW VERSION\n"); */
  /* double* neg_1o8_K2_psi_neg7_vec = P4EST_ALLOC(double, prob_vecs->local_nodes); */

  double* M_neg_1o8_K2_psi_neg7_vec = P4EST_ALLOC(double, prob_vecs->local_nodes);
  
  /* element_data_compute_f_of_uxyz */
  /*   ( */
  /*    p4est, */
  /*    prob_vecs->u, */
  /*    neg_1o8_K2_psi_neg7_vec, */
  /*    neg_1o8_K2_psi_neg7, */
  /*    d4est_ops */
  /*   ); */

  /* element_data_apply_mij_on_vec */
  /*   ( */
  /*    p4est, */
  /*    neg_1o8_K2_psi_neg7_vec, */
  /*    M_neg_1o8_K2_psi_neg7_vec, */
  /*    d4est_ops */
  /*   ); */


  element_data_apply_mij_on_f_of_vec
    (
     p4est,
     prob_vecs->u,
     M_neg_1o8_K2_psi_neg7_vec,
     d4est_ops,
     neg_1o8_K2_psi_neg7,
     1
    );

  d4est_linalg_vec_axpy(1.0, M_neg_1o8_K2_psi_neg7_vec, prob_vecs->Au, prob_vecs->local_nodes);

  /* P4EST_FREE(neg_1o8_K2_psi_neg7_vec); */
  P4EST_FREE(M_neg_1o8_K2_psi_neg7_vec);
  
}

static
void
build_residual_gauss
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_elliptic_problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops
)
{
  prob_vecs->scalar_flux_fcn_data.bndry_fcn = boundary_fcn;
  prob_vecs->vector_flux_fcn_data.bndry_fcn = boundary_fcn;
  
  poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops);  

  double* M_neg_1o8_K2_psi_neg7_vec= P4EST_ALLOC(double, prob_vecs->local_nodes);
 
  problem_ctx_t* ctx = (problem_ctx_t*)prob_vecs->user;

  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        element_data_t* ed = quad->p.user_data;        
        int deg_gauss = ed->deg + ctx->deg_offset_for_gauss_quad;

        int volume_nodes_gauss = d4est_lgl_get_nodes((P4EST_DIM), deg_gauss);
        double* r_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 0);
        double* s_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 1);
        double* t_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 2);
        double *jac_gauss = P4EST_ALLOC(double, volume_nodes_gauss);

        double* x_GL = P4EST_ALLOC(double, volume_nodes_gauss);
        double* y_GL = P4EST_ALLOC(double, volume_nodes_gauss);
        double* z_GL = P4EST_ALLOC(double, volume_nodes_gauss);

        d4est_linalg_fill_vec(jac_gauss, ed->jacobian, volume_nodes_gauss);  
        d4est_reference_rtox_array(r_GL, ed->xyz_corner[0], ed->h, x_GL, volume_nodes_gauss);
        d4est_reference_rtox_array(s_GL, ed->xyz_corner[1], ed->h, y_GL, volume_nodes_gauss);
        d4est_reference_rtox_array(t_GL, ed->xyz_corner[2], ed->h, z_GL, volume_nodes_gauss);
        double* xyz_gauss [3] = {x_GL, y_GL, z_GL};
  
        d4est_operators_apply_fofufofvlj_gaussnodes
          (
           d4est_ops,
           &prob_vecs->u[ed->stride],
           NULL,
           ed->deg,
           jac_gauss,
           xyz_gauss,
           deg_gauss,
           (P4EST_DIM),
           &M_neg_1o8_K2_psi_neg7_vec[ed->stride],
           neg_1o8_K2_psi_neg7,
           NULL,
           NULL,
           NULL
          );

        P4EST_FREE(z_GL);
        P4EST_FREE(y_GL);
        P4EST_FREE(x_GL);
        P4EST_FREE(jac_gauss);
      }
    }


  d4est_linalg_vec_axpy(1.0, M_neg_1o8_K2_psi_neg7_vec, prob_vecs->Au, prob_vecs->local_nodes);

  P4EST_FREE(M_neg_1o8_K2_psi_neg7_vec); 
}





static
void apply_jac
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_elliptic_problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops
)
{

  /* printf("HI FROM JAC\n"); */
  /* -u'' */
  poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops);

  /* integral(f(u)) =+ Au */
  /* element_data_quadrate_fofuv_andaddto(p4est, */
  /*                                       prob_vecs->u0, */
  /*                                       prob_vecs->u, */
  /*                                       prob_vecs->Au, */
  /*                                       /\* -1, *\/ */
  /*                                       plus_7o8_K2_psi_neg8, */
  /*                                      d4est_ops); */


  /* double* plus_7o8_K2_psi_neg8_of_u0_vec = P4EST_ALLOC(double, prob_vecs->local_nodes); */
  /* double* plus_7o8_K2_psi_neg8_of_u0_u_vec = P4EST_ALLOC(double, prob_vecs->local_nodes); */
  double* M_plus_7o8_K2_psi_neg8_of_u0_u_vec = P4EST_ALLOC(double, prob_vecs->local_nodes);
  
  /* element_data_compute_f_of_uxyz */
  /*   ( */
  /*    p4est, */
  /*    prob_vecs->u0, */
  /*    plus_7o8_K2_psi_neg8_of_u0_vec, */
  /*    plus_7o8_K2_psi_neg8, */
  /*    d4est_ops */
  /*   ); */

  /* d4est_linalg_component_mult */
  /*   ( */
  /*    plus_7o8_K2_psi_neg8_of_u0_vec, */
  /*    prob_vecs->u, */
  /*    plus_7o8_K2_psi_neg8_of_u0_u_vec, */
  /*    prob_vecs->local_nodes */
  /*   ); */

  /* element_data_apply_mij_on_vec */
  /*   ( */
  /*    p4est, */
  /*    plus_7o8_K2_psi_neg8_of_u0_u_vec, */
  /*    M_plus_7o8_K2_psi_neg8_of_u0_u_vec, */
  /*    d4est_ops */
  /*   ); */


  element_data_apply_mij_on_f_of_vec1_x_vec2
    (
     p4est,
     prob_vecs->u0,
     prob_vecs->u,
     M_plus_7o8_K2_psi_neg8_of_u0_u_vec,
     d4est_ops,
     plus_7o8_K2_psi_neg8,
     1
    );
  
  d4est_linalg_vec_axpy(1.0, M_plus_7o8_K2_psi_neg8_of_u0_u_vec, prob_vecs->Au, prob_vecs->local_nodes);

  /* P4EST_FREE(plus_7o8_K2_psi_neg8_of_u0_u_vec); */
  /* P4EST_FREE(plus_7o8_K2_psi_neg8_of_u0_vec); */
  P4EST_FREE(M_plus_7o8_K2_psi_neg8_of_u0_u_vec);
}


static
void apply_jac_gauss
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_elliptic_problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops
)
{
  prob_vecs->scalar_flux_fcn_data.bndry_fcn = zero_fcn;
  prob_vecs->vector_flux_fcn_data.bndry_fcn = zero_fcn;
  
  poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops);

  /* double* neg_10pi_rho_up1_neg4_of_u0_vec = P4EST_ALLOC(double, prob_vecs->local_nodes); */
  /* double* neg_10pi_rho_up1_neg4_of_u0_u_vec = P4EST_ALLOC(double, prob_vecs->local_nodes); */
  double* M_plus_7o8_K2_psi_neg8_of_u0_u_vec = P4EST_ALLOC(double, prob_vecs->local_nodes);

  problem_ctx_t* ctx = (problem_ctx_t*)prob_vecs->user;

  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        element_data_t* ed = quad->p.user_data;        
        int deg_gauss = ed->deg + ctx->deg_offset_for_gauss_quad;

        int volume_nodes_gauss = d4est_lgl_get_nodes((P4EST_DIM), deg_gauss);
        double* r_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 0);
        double* s_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 1);
        double* t_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 2);
        double *jac_gauss = P4EST_ALLOC(double, volume_nodes_gauss);

        double* x_GL = P4EST_ALLOC(double, volume_nodes_gauss);
        double* y_GL = P4EST_ALLOC(double, volume_nodes_gauss);
        double* z_GL = P4EST_ALLOC(double, volume_nodes_gauss);

        d4est_linalg_fill_vec(jac_gauss, ed->jacobian, volume_nodes_gauss);  
        d4est_reference_rtox_array(r_GL, ed->xyz_corner[0], ed->h, x_GL, volume_nodes_gauss);
        d4est_reference_rtox_array(s_GL, ed->xyz_corner[1], ed->h, y_GL, volume_nodes_gauss);
        d4est_reference_rtox_array(t_GL, ed->xyz_corner[2], ed->h, z_GL, volume_nodes_gauss);
        double* xyz_gauss [3] = {x_GL, y_GL, z_GL};
  
        d4est_operators_apply_fofufofvlilj_gaussnodes
          (
           d4est_ops,
           &prob_vecs->u[ed->stride],
           &prob_vecs->u0[ed->stride],
           NULL,
           ed->deg,
           jac_gauss,
           xyz_gauss,
           deg_gauss,
           (P4EST_DIM),
           &M_plus_7o8_K2_psi_neg8_of_u0_u_vec[ed->stride],
           plus_7o8_K2_psi_neg8,
           NULL,
           NULL,
           NULL
          );
        P4EST_FREE(z_GL);
        P4EST_FREE(y_GL);
        P4EST_FREE(x_GL);
        P4EST_FREE(jac_gauss);
      }
    }

  d4est_linalg_vec_axpy(1.0, M_plus_7o8_K2_psi_neg8_of_u0_u_vec, prob_vecs->Au, prob_vecs->local_nodes);

  P4EST_FREE(M_plus_7o8_K2_psi_neg8_of_u0_u_vec);

  prob_vecs->scalar_flux_fcn_data.bndry_fcn = boundary_fcn;
  prob_vecs->vector_flux_fcn_data.bndry_fcn = boundary_fcn;
}

/* static */
/* void problem_build_rhs */
/* ( */
/*  p4est_t* p4est, */
/*  d4est_elliptic_problem_data_t* prob_vecs, */
/*  d4est_elliptic_eqns_t* prob_fcns, */
/*  p4est_ghost_t* ghost, */
/*  element_data_t* ghost_data, */
/*  penalty_calc_t peanalty_fcn, */
/*  d4est_operators_t* d4est_ops, */
/*  double sipg_flux_prefactor */
/* ) */
/* { */
/*   d4est_linalg_fill_vec(prob_vecs->rhs,0.0,prob_vecs->local_nodes);   */
/* } */

p4est_connectivity_t*
problem_build_conn()
{
  return p8est_connectivity_new_unitcube();
}

p4est_geometry_t*
problem_build_geom
(
 p4est_connectivity_t* conn
)
{
  return NULL;
}

p4est_t*
problem_build_p4est
(
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t* conn,
 p4est_locidx_t min_quadrants,
 int min_level, 
 int fill_uniform
)
{
  return p4est_new_ext
    (
     mpicomm,
     conn,
     min_quadrants,
     min_level,
     fill_uniform,
     sizeof(element_data_t),
     NULL,
     NULL
    );
}

p4est_t*
problem_load_p4est_from_checkpoint
(
 const char* filename,
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t** conn
){
  int autopartition = 1;
  int load_data = 1;
  int broadcasthead = 0;
  
  return p4est_load_ext (filename,
                         mpicomm,
                         sizeof(element_data_t),
                         load_data,
                         autopartition,
                         broadcasthead,
                         NULL,
                         conn);
}


typedef struct {

  int endlevel;
  int degree;
  double gamma_h;
  double gamma_p;
  double domain_size;
  double rho0_div_rhoc;
  double ip_flux_penalty;
  int hrefine_til_inview;
  int percentile;
  int use_gauss_quad;
  int deg_offset_for_gauss_quad;
  int degmax;
  KSPType krylov_type;
  int amr_inflation_size;
  
  int count;
  
} problem_input_t;


static
int problem_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  problem_input_t* pconfig = (problem_input_t*)user;
  if (util_match_couple(section,"amr",name,"amr_levels")) {
    mpi_assert(pconfig->endlevel == -1);
    pconfig->endlevel = atoi(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"amr",name, "initial_degree")) {
    mpi_assert(pconfig->degree == -1);
    pconfig->degree = atoi(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"amr",name,"percentile")) {
    mpi_assert(pconfig->percentile == -1);
    pconfig->percentile = atoi(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"amr",name,"gamma_h")) {
    mpi_assert(pconfig->gamma_h == -1);
    pconfig->gamma_h = atof(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"amr",name,"gamma_p")) {
    mpi_assert(pconfig->gamma_p == -1);
    pconfig->gamma_p = atof(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"amr",name,"degmax")) {
    mpi_assert(pconfig->degmax == -1);
    pconfig->degmax = atoi(value);
    pconfig->count += 1;
  }else if (util_match_couple(section,"flux",name,"ip_flux_penalty")) {
    mpi_assert(pconfig->ip_flux_penalty == -1);
    pconfig->ip_flux_penalty = atof(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"problem",name,"hrefine_til_inview")) {
    mpi_assert(pconfig->hrefine_til_inview == -1);
    pconfig->hrefine_til_inview = atoi(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"problem",name,"domain_size")) {
    mpi_assert(pconfig->domain_size == -1);
    pconfig->domain_size = atof(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"use_gauss_quad")) {
    mpi_assert(pconfig->use_gauss_quad == -1);
    pconfig->use_gauss_quad = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_offset_for_gauss_quad")) {
    mpi_assert(pconfig->deg_offset_for_gauss_quad == -1);
    pconfig->deg_offset_for_gauss_quad = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"amr",name,"amr_inflation_size")) {
    mpi_assert(pconfig->amr_inflation_size == -1);
    pconfig->amr_inflation_size = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"solver",name,"krylov_type")) {

    if (strcmp("cg", value) == 0){
      pconfig->krylov_type = KSPCG;
    }
    else if (strcmp("gmres", value) == 0){
      pconfig->krylov_type = KSPGMRES;
    }
    else {
      mpi_abort("not a support KSP solver type");
    }
    pconfig->count += 1;
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

static
problem_input_t
problem_input
(
 const char* input_file
)
{
  int num_of_options = 13;
  
  problem_input_t input;
  input.degree = -1;
  input.degmax = -1;
  input.domain_size = -1;
  input.endlevel = -1;
  input.gamma_h = -1;
  input.gamma_p = -1;
  input.ip_flux_penalty = -1;
  input.percentile = -1;
  input.use_gauss_quad = -1;
  input.deg_offset_for_gauss_quad = -1;
  input.amr_inflation_size = -1;
  input.hrefine_til_inview = -1;
  input.krylov_type = "";
  input.count = 0;
  
  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  if (input.count != num_of_options){
    printf("[D4EST_INFO]: input.count = %d\n", input.count);
    printf("[D4EST_INFO]: num_of_options = %d\n", num_of_options);
    mpi_abort("num_of_options != input.count\n");
  }
  return input;
}


void
problem_init
(
 int argc,
 char* argv [],
 p4est_t* p4est,
 p4est_geometry_t* p4est_geom,
 d4est_operators_t* d4est_ops,
 int proc_size,
 sc_MPI_Comm mpicomm,
 int load_from_checkpoint
)
{

  mpi_assert((P4EST_DIM) == 3);
  /* mpi_assert(argc == 12); */
  
  int world_rank,world_size;
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &world_rank);
  sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &world_size);

  problem_input_t input = problem_input("options.input");
  
  int level;
  int endlevel = input.endlevel;
  int degree = input.degree;         
  double gamma_h = input.gamma_h;
  double gamma_p = input.gamma_p;
  double domain_size = input.domain_size;
  int percentile = input.percentile;
  int inflation_percentile = 25;
  int hrefine_til_inview = input.hrefine_til_inview;
  
  if (world_rank == 0){
    printf("\n");
    printf("HP_LEVS = %d\n", endlevel);
    printf("INIT_DEG = %d\n", degree);
    printf("PENALTY = %f\n", input.ip_flux_penalty);
    printf("SMOOTH_PRED_GAMMA_H = %f\n", gamma_h);
    printf("SMOOTH_PRED_GAMMA_P = %f\n", gamma_p);
    printf("DOMAIN_SIZE = %f\n", domain_size);
    printf("INITIAL NUMBER OF ELEMENTS PER PROC = %d\n", p4est->local_num_quadrants);
  }
  
  /* int save_hp_img = 1; */
  /* int save_hp_img_period = 1; */
  /* int save_u_puncture = 4; */
  
  double M = 1./domain_size;
  
  M_bh[0] = .5*M;
  M_bh[1] = .5*M;

  xyz_bh[0][0] = .5 - 3*M;
  xyz_bh[0][1] = 0.5;
  xyz_bh[0][2] = 0.5;

  xyz_bh[1][0] = .5 + 3*M;
  xyz_bh[1][1] = 0.5;
  xyz_bh[1][2] = 0.5;

  P_bh[0][0] = 0.;
  P_bh[0][1] = -0.2*M;
  P_bh[0][2] = 0.;

  P_bh[1][0] = 0.;
  P_bh[1][1] = 0.2*M;
  P_bh[1][2] = 0.;

  S_bh[0][0] = 0.;
  S_bh[0][1] = 0.;
  S_bh[0][2] = 0.;
  
  S_bh[1][0] = 0.;
  S_bh[1][1] = 0.;
  S_bh[1][2] = 0.;

  double sep_dx = xyz_bh[0][0] - xyz_bh[1][0];
  double sep_dy = xyz_bh[0][1] - xyz_bh[1][1];
  double sep_dz = xyz_bh[0][2] - xyz_bh[1][2];
  double separation = sqrt(sep_dx*sep_dx + sep_dy*sep_dy + sep_dz*sep_dz);
  double global_h_min = 1.;
  
  double adm_energy_analytic = M_bh[0] + M_bh[1] + (5./8.)*P_bh[0][1]*P_bh[0][1]/M_bh[0] + (5./8.)*P_bh[1][1]*P_bh[1][1]/M_bh[1];
  printf("ADM energy theoretical estimate = %f\n", adm_energy_analytic);

  p4est_partition_ext(p4est, 0, NULL);
  p4est_balance_ext(p4est, P4EST_CONNECT_FACE, NULL, NULL);
  
  int local_nodes;
  p4est_topidx_t      t, flt, llt;
  p4est_tree_t       *tree;
  int init_max_level = 0;

  /* compute the timestep by finding the smallest quadrant */
  flt = p4est->first_local_tree;
  llt = p4est->last_local_tree;

  for (t = flt; t <= llt; t++) {
    tree = p4est_tree_array_index (p4est->trees, t);
    init_max_level = SC_MAX (init_max_level, tree->maxlevel);
  }
  
  p4est_reset_data(p4est, sizeof(element_data_t), NULL, NULL);
  element_data_init(p4est, degree);

  local_nodes = element_data_get_local_nodes(p4est);
  double* Au = P4EST_ALLOC_ZERO(double, local_nodes);
  /* double* rhs = P4EST_ALLOC_ZERO(double, local_nodes); */
  double* u = P4EST_ALLOC_ZERO(double, local_nodes);
  /* double* f = P4EST_ALLOC_ZERO(double, local_nodes); */
  double* u_prev = P4EST_ALLOC_ZERO(double, local_nodes);

  d4est_linalg_fill_vec(u, 0., local_nodes);

  double local_eta2 = -1.;

  penalty_calc_t bi_u_penalty_fcn = bi_u_prefactor_conforming_maxp_minh;
  penalty_calc_t bi_u_dirichlet_penalty_fcn = bi_u_prefactor_conforming_maxp_minh;
  penalty_calc_t bi_gradu_penalty_fcn = bi_gradu_prefactor_maxp_minh;

  ip_flux_params_t ip_flux_params;
  ip_flux_params.ip_flux_penalty_prefactor = input.ip_flux_penalty;
  ip_flux_params.ip_flux_penalty_calculate_fcn = sipg_flux_vector_calc_penalty_maxp2_over_minh;
  
  d4est_elliptic_problem_data_t prob_vecs;
  /* prob_vecs.rhs = rhs; */
  prob_vecs.Au = Au;
  prob_vecs.u = u;
  /* prob_vecs.f = f; */
  prob_vecs.local_nodes = local_nodes;
  prob_vecs.vector_flux_fcn_data = sipg_flux_vector_dirichlet_fetch_fcns
                                   (
                                    zero_fcn,
                                    &ip_flux_params
                                   );
  prob_vecs.scalar_flux_fcn_data = sipg_flux_scalar_dirichlet_fetch_fcns(zero_fcn);


  /* set weak equations */
  d4est_elliptic_eqns_t prob_fcns;
  problem_ctx_t prob_ctx;
  if (input.use_gauss_quad){
    mpi_assert(input.deg_offset_for_gauss_quad > -1);
    prob_ctx.deg_offset_for_gauss_quad = input.deg_offset_for_gauss_quad;
    prob_fcns.apply_lhs = apply_jac_gauss;
    prob_fcns.build_residual = build_residual_gauss;
    prob_vecs.user = (void*)&prob_ctx;
  }
  else{
    prob_fcns.apply_lhs = apply_jac;
    prob_fcns.build_residual = build_residual;
  }

  hp_amr_scheme_t* scheme =
    hp_amr_smooth_pred_init
    (
     p4est,
     gamma_h,
     gamma_p,
     1.,
     input.degmax,
     hp_amr_smooth_pred_get_NULL_marker()
    );

  /* force a balance if we're loading from checkpoint with autopartition */
  if(load_from_checkpoint){
    mpi_abort("Not supported yet");
  }
  
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  /* create space for storing the ghost data */
  element_data_t* ghost_data = P4EST_ALLOC (element_data_t,
                                            ghost->ghosts.elem_count);
  

  /* element_data_init_node_vec(p4est, f, f_fcn, d4est_ops); */
  /* problem_build_rhs */
  /*   ( */
  /*    p4est, */
  /*    &prob_vecs, */
  /*    &prob_fcns, */
  /*    ghost, */
  /*    ghost_data, */
  /*    penalty_fcn, */
  /*    d4est_ops, */
  /*    sipg_flux_prefactor */
  /*   );    */

  double* dof_data_for_fit = P4EST_ALLOC(double, endlevel);
  double* dgerr_data_for_fit = P4EST_ALLOC(double, endlevel);
  
  for (level = 0; level < endlevel; ++level){

    /* int local_check = 1; */
    /* int global_check = 0; */
    /* sc_reduce */
    /*   ( */
    /*    &local_check, */
    /*    &global_check, */
    /*    1, */
    /*    sc_MPI_INT, */
    /*    sc_MPI_SUM, */
    /*    0, */
    /*    sc_MPI_COMM_WORLD */
    /*   ); */

    /* if (world_rank == 0){ */
    /*   printf(" We have %d processes at this point \n", global_check); */
    /* } */


    bi_estimator_compute
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       bi_u_penalty_fcn,
       bi_u_dirichlet_penalty_fcn,
       bi_gradu_penalty_fcn,
       zero_fcn,
       ip_flux_params.ip_flux_penalty_prefactor,
       ghost,
       ghost_data,
       d4est_ops
      );

    
    estimator_stats_t stats;
    estimator_stats_compute(p4est, &stats,0);

    if(world_rank == 0)
      estimator_stats_print(&stats, 0);

    local_eta2 = stats.total;
       
    /* if(save_hp_img && (level && !(level % save_hp_img_period))){ */
    /*   char save_as [500]; */
    /*   sprintf(save_as, "%s_hp_amr_level_%d", P4EST_STRING,level-1); */
    /*   hp_amr_save_to_vtk */
    /*     ( */
    /*      p4est, */
    /*      save_as, */
    /*      1 */
    /*     ); */
    /* } */

    double* u_vertex = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN));
    double* eta2_vertex = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN));
    element_data_store_local_estimator_in_corner_array
      (
       p4est,
       eta2_vertex
      );

    element_data_store_nodal_vec_in_vertex_array
      (
       p4est,
       u,
       u_vertex
      );
    
    char sol_save_as [500];
    sprintf(sol_save_as, "%s_hp_amr_level_%d_sols", P4EST_STRING, level);
    hacked_p4est_vtk_write_all
      (p4est,
       NULL,
       1.,
       0,   
       1,   
       1,   
       0,
       2,
       0,  
       sol_save_as,
       "eta2",
       eta2_vertex,
       "psi",
       u_vertex
      );

    P4EST_FREE(u_vertex);   
    P4EST_FREE(eta2_vertex);   
    

    if (world_rank == 0)
      printf("[INFO]: global_h_min = %f, separation = %f\n", global_h_min, separation);

    printf("p4est->local_num_quadrants*p4est->mpisize < input.amr_inflation_size -> %d < %d = %d\n", p4est->local_num_quadrants*p4est->mpisize,input.amr_inflation_size, p4est->local_num_quadrants*p4est->mpisize < input.amr_inflation_size);
    if (p4est->local_num_quadrants*p4est->mpisize < input.amr_inflation_size){
      hp_amr_smooth_pred_set_marker
        (
         scheme,
         hp_amr_smooth_pred_get_percentile_marker(&inflation_percentile)
        );
    }
    else {
      hp_amr_smooth_pred_set_marker
        (
         scheme,
         hp_amr_smooth_pred_get_percentile_marker(&percentile)
        );
    }

    if (global_h_min > separation && hrefine_til_inview){
      printf("[INFO]: Still only h-refining\n");
      hp_amr_smooth_pred_set_gammah(scheme, 0.);
    }
    else {
      printf("[INFO]: Starting to p-refine\n");
      hp_amr_smooth_pred_set_gammah(scheme, gamma_h);
    }
    
    hp_amr(p4est,
           d4est_ops,
           &u,
           &stats,
           scheme
          );
        
    
    /* if (level < endlevel - 1) */
    /* hp_amr_smooth_pred_save_predictor(p4est, smooth_pred_data); */
        
    p4est_ghost_destroy(ghost);
    P4EST_FREE(ghost_data);

    ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
    ghost_data = P4EST_ALLOC(element_data_t, ghost->ghosts.elem_count);
    
    element_data_init(p4est, -1);
    local_nodes = element_data_get_local_nodes(p4est);

    Au = P4EST_REALLOC(Au, double, local_nodes);
    /* f = P4EST_REALLOC(f, double, local_nodes); */
    /* rhs = P4EST_REALLOC(rhs, double, local_nodes); */
    u_prev = P4EST_REALLOC(u_prev, double, local_nodes);

    /* initialize vectors */
    /* element_data_init_node_vec(p4est, f, f_fcn, d4est_ops); */
  
    /* prob_vecs.rhs = rhs; */
    prob_vecs.Au = Au;
    prob_vecs.u = u;
    prob_vecs.u0 = u;
    /* prob_vecs.f = f; */
    prob_vecs.local_nodes = local_nodes;
    prob_vecs.vector_flux_fcn_data = sipg_flux_vector_dirichlet_fetch_fcns
                                     (
                                      zero_fcn,
                                      &ip_flux_params
                                     );
  
    prob_vecs.scalar_flux_fcn_data = sipg_flux_scalar_dirichlet_fetch_fcns(zero_fcn);


    /* problem_build_rhs */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs, */
    /*    &prob_fcns, */
    /*    ghost, */
    /*    ghost_data, */
    /*    penalty_fcn, */
    /*    d4est_ops, */
    /*    sipg_flux_prefactor */
    /*   ); */
    
    clock_t begin = 0;
    clock_t end = -1;

    if (world_rank == 0){
      begin = clock();
      /* util_print_matrix(u, local_nodes, 1, "u = ", 0); */
    }

    flt = p4est->first_local_tree;
    llt = p4est->last_local_tree;

    int cur_max_level = 0;
    for (t = flt; t <= llt; t++) {
      tree = p4est_tree_array_index (p4est->trees, t);
      cur_max_level = SC_MAX (cur_max_level, tree->maxlevel);
    }

    int local_num_levels = cur_max_level - init_max_level + 1;
    int global_num_levels = -1.;

    sc_allreduce
      (
       &local_num_levels,
       &global_num_levels,
       1,
       sc_MPI_INT,
       sc_MPI_MIN,
       sc_MPI_COMM_WORLD
      );

    /* printf("global_num_levels = %d\n", global_num_levels); */
    
    /* multigrid_nr_data_t mg_data; */
    /* mg_data.num_of_levels = global_num_levels; */
    /* mg_data.lmax_lmin_rat = 30.; */
    /* mg_data.save_vtk_snapshot = 0; */
    /* mg_data.perform_checksum = 0; */
    /* mg_data.fine_nodes = local_nodes; */
    /* mg_data.coarse_nodes = -1; */
    /* mg_data.mpi_rank = world_rank; */
    /* mg_data.smooth_iter = 8; */
    /* mg_data.coarse_iter = 1000; */
    /* mg_data.coarse_rtol = 1e-10; */
    /* mg_data.vcycle_rtol = 1e-9; */
    /* mg_data.vcycle_imax = 5; */
    /* mg_data.cg_eigs_iter = 10; */
    /* mg_data.d4est_ops = d4est_ops; */

    /* newton_multigrid_solver_params_t params; */
    /* params.atol = 1e-8; */
    /* params.rtol = 1e-5; */
    /* params.max_iter = 5; */
    /* params.mpi_rank = world_rank; */
    /* params.eta_max = -.01; */
    /* params.gamma = .9; */
    /* params.mg_data = &mg_data; */

    /* newton_multigrid_solve */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs, */
    /*    &prob_fcns, */
    /*    &ghost, */
    /*    &ghost_data, */
    /*    (void*)&params, */
    /*    d4est_ops */
    /*   ); */

    /* newton_petsc_solve */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs, */
    /*    &prob_fcns, */
    /*    NULL, */
    /*    d4est_ops */
    /*   ); */

    /* multigrid_nr_data_t mg_data; */
    /* mg_data.num_of_levels = 2 + level;//global_num_levels; */
    /* mg_data.lmax_lmin_rat = 30.; */
    /* mg_data.save_vtk_snapshot = 0; */
    /* mg_data.perform_checksum = 0; */
    /* mg_data.fine_nodes = local_nodes; */
    /* mg_data.coarse_nodes = -1; */
    /* mg_data.mpi_rank = world_rank; */
    /* mg_data.smooth_iter = 15; */
    /* mg_data.coarse_iter = 100; */
    /* mg_data.coarse_rtol = 1e-10; */
    /* mg_data.vcycle_rtol = 1e-9; */
    /* mg_data.vcycle_imax = 5 + level; */
    /* mg_data.cg_eigs_iter = 20; */
    /* mg_data.d4est_ops = d4est_ops; */
    
    /* newton_multigrid_nols_solver_params_t params; */
    /* params.atol = 1e-18; */
    /* params.rtol = 1e-4; */
    /* params.max_iter = 3 + level; */
    /* /\* params.max_iter_linear = util_int_pow_int(10, level+1); *\/ */
    /* params.eta_max = .01; */
    /* params.mg_data = &mg_data; */
    /* params.mpi_rank = world_rank; */
    
    /* newton_multigrid_nols_solve */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs, */
    /*    &prob_fcns, */
    /*    &ghost, */
    /*    &ghost_data, */
    /*    (void*)&params, */
    /*    d4est_ops */
    /*   );    */

    d4est_linalg_copy_1st_to_2nd(u, u_prev, local_nodes);
    
    /* newton_krylov_nols_solver_params_t params; */
    /* params.atol = 1e-18; */
    /* params.rtol = newton_rtol; */
    /* params.max_iter = 3 + level; */
    /* params.max_iter_linear = util_int_pow_int(10, level+1); */
    /* params.eta_max = newton_eta_max; */
    /* /\* params.gamma = .9; *\/ */
    /* params.mpi_rank = world_rank; */
    
    /* newton_krylov_nols_solve */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs, */
    /*    &prob_fcns, */
    /*    (void*)&params, */
    /*    ghost, */
    /*    ghost_data, */
    /*    d4est_ops */
    /*   ); */

    /* newton_krylov_nols_solver_params_t params; */
    /* params.final_fnrm = -1.; */
    /* params.final_iter = -1; */
    
    /* newton_petsc_solve */
    /* ( */
    /*    p4est, */
    /*    &prob_vecs, */
    /*    &prob_fcns, */
    /*    ghost, */
    /*    ghost_data, */
    /*    d4est_ops, */
    /*    0 */
    /*   ); */


    krylov_petsc_params_t params;
    params.user_defined_pc = 0;
    params.ksp_monitor = 0;
    params.krylov_type = input.krylov_type;

    newton_petsc_solve
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       &ghost,
       &ghost_data,
       d4est_ops,
       &params
      );
    
    
    
    d4est_linalg_vec_axpy(-1., u, u_prev, local_nodes);

    
    double local_adm_energy = element_data_compute_boundary_quadral
                              (
                               p4est,
                               two_punctures_adm_quadral_face_contribution,
                               u,
                               d4est_ops
                              );

    double local_adm_energy2 = two_punctures_adm_quadral_volume
                               (
                                u,
                                M_bh[0],
                                M_bh[1],
                                p4est,
                                ghost,
                                ghost_data,
                                &prob_vecs,
                                d4est_ops
                               );
    
    /* dg norm should always have the boundary fcn set to zero */
    double local_dg_norm_sqr = element_data_compute_DG_norm_sqr
                               (
                                p4est,
                                u_prev,
                                zero_fcn,
                                &ip_flux_params,
                                d4est_ops,
                                ghost,
                                ghost_data
                               );
    
    double local_l2_norm_sqr =  element_data_compute_l2_norm_sqr_no_local
                                (
                                 p4est,
                                 u_prev,
                                 d4est_ops
                                );


    double local_h_min = element_data_get_local_h_min
                         (
                          p4est
                         );

    double local_nodes_dbl = (double)local_nodes;
    double local_reduce [6];
    local_reduce[0] = local_nodes_dbl;
    local_reduce[1] = local_l2_norm_sqr;
    local_reduce[2] = local_dg_norm_sqr;
    local_reduce[3] = local_eta2;
    local_reduce[4] = local_adm_energy;
    local_reduce[5] = local_adm_energy2;

    /* printf("local_eta2 = %f\n",local_eta2); */
    
    double global_reduce [6];

    sc_reduce
      (
       &local_reduce[0],
       &global_reduce[0],
       6,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       0,
       sc_MPI_COMM_WORLD
      );


    sc_allreduce
      (
       &local_h_min,
       &global_h_min,
       1,
       sc_MPI_DOUBLE,
       sc_MPI_MIN,
       sc_MPI_COMM_WORLD
      );
    

    double global_nodes_dbl = global_reduce[0];
    double global_l2_norm_sqr = global_reduce[1];
    double global_dg_norm_sqr = global_reduce[2];
    double global_eta2 = global_reduce[3];
    double global_adm_energy = global_reduce[4];
    double global_adm_energy2 = global_reduce[5];
    
    if (world_rank == 0){
      end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
      printf
        (
         "[HP_AMR]: %d %d %d %d %.25f .25f %.25f %.25f %.25f %f\n",
         level,
         degree,
         (int)p4est->global_num_quadrants,
         (int)global_nodes_dbl,
         sqrt(global_eta2),
         sqrt(global_l2_norm_sqr),
         sqrt(global_dg_norm_sqr),
         global_adm_energy,
         global_adm_energy2,
         /* params.final_fnrm, */
         /* params.final_iter, */
         time_spent
        );

      dgerr_data_for_fit[level] = log(sqrt(global_dg_norm_sqr));
      dof_data_for_fit[level] = pow(global_nodes_dbl, 1./(2.*(P4EST_DIM)-1.));

      if (level > 0){
        double slope;
        double intercept;
        int num_of_hpamr_levels = level + 1;
        util_linear_regression
          (
           dgerr_data_for_fit,
           dof_data_for_fit,
           &slope,
           &intercept,
           num_of_hpamr_levels
          );
        printf("[HP_AMR_FIT](1): ||err||DG = C1*exp(-C2*DOF^(1/%d))\n",2*(P4EST_DIM)-1);
        printf("[HP_AMR_FIT](2): LEV SLOPE DG_ERR\n");
        printf("[HP_AMR_FIT](3): %d %.25f %.25f, \n\n", level, slope, sqrt(global_dg_norm_sqr));
      }      


      
    }


    int checkpoint_period = 3;
    if (level > 8)
      checkpoint_period = 1;
    
    if (level && level % checkpoint_period == 0){
      /* char* p4est_filename = NULL; */
      /* D4EST_ASPRINTF(p4est_filename, "%s_%s_%d", "checkpoint","level", level); */
      /* element_data_copy_from_vec_to_storage(p4est, u); */
      /* p4est_save_ext(p4est_filename, p4est, 1, 1); */
      /* free(p4est_filename); */

      sc_MPI_Barrier(sc_MPI_COMM_WORLD);
      char save_as1 [500];
      sprintf(save_as1, "two_punctures_u_3d_hp_amr_level_%d", level+1);
      element_data_print_node_vec
        (
         p4est,
         u,
         &save_as1[0],
         world_rank,
         1 /* save to file */,
         d4est_ops
        );      
    }
  }

  double* u_vertex = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN));

  element_data_store_nodal_vec_in_vertex_array
    (
     p4est,
     u,
     u_vertex
    );
  
  char sol_save_as [500];
  sprintf(sol_save_as, "%s_hp_amr_level_%d_sols_noeta2", P4EST_STRING, level);

  hacked_p4est_vtk_write_all
    (p4est,
     NULL,
     1.,
     0,   
     1,   
     1,   
     0,
     1,
     0,  
     sol_save_as,
     "u",
     u_vertex
    );

  P4EST_FREE(u_vertex);   
  
  hp_amr_smooth_pred_destroy(scheme);
  
  /* if (ghost) { */
  p4est_ghost_destroy (ghost);
  P4EST_FREE (ghost_data);
  ghost = NULL;
  ghost_data = NULL;
  /* } */
  
  /* P4EST_FREE(f); */
  P4EST_FREE(dof_data_for_fit);
  P4EST_FREE(dgerr_data_for_fit);  
  P4EST_FREE(Au);
  /* P4EST_FREE(rhs); */
  P4EST_FREE(u_prev);
  P4EST_FREE(u);
}
