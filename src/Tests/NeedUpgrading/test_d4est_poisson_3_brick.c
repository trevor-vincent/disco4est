#include <pXest.h>
#include <d4est_quadrature.h>
#include <d4est_element_data.h>
#include <d4est_quadrature_legendre.h>
#include <d4est_poisson_flux_sipg.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_geometry_brick.h>
#include <petscsnes.h>
#include <d4est_linalg.h>
#include <d4est_mortars.h>
#include <d4est_amr.h>
#include <d4est_poisson.h>
#include <d4est_poisson_flux.h>
#include <d4est_util.h>
#include <limits.h>

#define D4EST_REAL_EPS 100*1e-15
#define TEST_DEG_INIT 1

typedef struct {

  double boundary_term_3_err;
  double boundary_term_2_err;
  double boundary_term_1_err;
  double interface_term_1_err;
  d4est_poisson_flux_sipg_params_t* sipg_params;
  
} testd4est_poisson_2_brick_data_t;

static void
problem_set_degrees_init
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  /* d4est_element_data_t* elem_data = elem_data_tmp; */
  elem_data->deg = TEST_DEG_INIT;
  elem_data->deg_quad = TEST_DEG_INIT;
}

double
poly_vec_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
){
  int deg = TEST_DEG_INIT;
  double poly = pow(x,deg) + pow(y,deg);
#if (P4EST_DIM)==3
  poly += ((P4EST_DIM)==3) ? pow(z,deg) : 0.;
#endif
  return poly;
}

double
laplacian_poly_vec_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
){
  
  int deg = TEST_DEG_INIT;
  if (deg-2 < 0){
    return 0.;
  }
  double poly = pow(x,deg-2) + pow(y,deg-2);
#if (P4EST_DIM)==3
  poly += ((P4EST_DIM)==3) ? pow(z,deg-2) : 0.;
#endif
  double factor = deg;
  while (factor != 0){
    poly *= factor;
    factor -= 1;
  }
  return poly;
}

static void
problem_set_degrees_amr
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  elem_data->deg_quad = elem_data->deg;

}

static p4est_t*
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
     sizeof(d4est_element_data_t),
     NULL,
     NULL
    );
}

static double
testd4est_poisson_2_brick_dirichlet_term3
(
 int deg,
 int face,
 int i,
 double penalty,
 d4est_geometry_t* d4est_geom,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq
)
{
  d4est_geometry_brick_attr_t* brick_attrs = d4est_geom->user;
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;
  
  double x0 = brick_attrs->X0;
  double x1 = brick_attrs->X1;
  double y0 = brick_attrs->Y0;
  double y1 = brick_attrs->Y1;

  /* transform element corners to [X0,X1] X [YO,Y1] X [Z0,Z1] topological space */
  amin = (x1-x0)*amin + x0;
  amax = (x1-x0)*amax + x0;
  bmin = (y1-y0)*bmin + y0;
  bmax = (y1-y0)*bmax + y0;
  
  if (deg == 1){
    double term3_boundary [4][4];
    term3_boundary[0][0] = (sqrt(pow(amax - amin,-2))*sqrt(pow(bmax - bmin,2))*(6*amin + 2*bmax + 4*bmin)*penalty)/3.;
    term3_boundary[0][1] = 0;
    term3_boundary[0][2] = (sqrt(pow(amax - amin,-2))*sqrt(pow(bmax - bmin,2))*(6*amin + 4*bmax + 2*bmin)*penalty)/3.;
    term3_boundary[0][3] = 0;
    term3_boundary[1][0] = 0;
    term3_boundary[1][1] = sqrt(pow(amax - amin,-2))*sqrt(pow(bmax - bmin,2))*(2*amax + (2*(bmax + 2*bmin))/3.)*penalty;
    term3_boundary[1][2] = 0;
    term3_boundary[1][3] = sqrt(pow(amax - amin,-2))*sqrt(pow(bmax - bmin,2))*(2*amax + (2*(2*bmax + bmin))/3.)*penalty;
    term3_boundary[2][0] = (sqrt(pow(amax - amin,2))*sqrt(pow(bmax - bmin,-2))*(2*amax + 4*amin + 6*bmin)*penalty)/3.;
    term3_boundary[2][1] = (sqrt(pow(amax - amin,2))*sqrt(pow(bmax - bmin,-2))*(4*amax + 2*amin + 6*bmin)*penalty)/3.;
    term3_boundary[2][2] = 0;
    term3_boundary[2][3] = 0;
    term3_boundary[3][0] = 0;
    term3_boundary[3][1] = 0;
    term3_boundary[3][2] = (sqrt(pow(amax - amin,2))*(2*amax + 4*amin + 6*bmax)*sqrt(pow(bmax - bmin,-2))*penalty)/3.;
    term3_boundary[3][3] = (sqrt(pow(amax - amin,2))*(4*amax + 2*amin + 6*bmax)*sqrt(pow(bmax - bmin,-2))*penalty)/3.;
    return term3_boundary[face][i];
  }
  else {
    D4EST_ABORT("no other degrees supported atm\n");
  }
}

static double
testd4est_poisson_2_brick_dirichlet_term1
(
 int deg,
 int face,
 int i,
 double penalty,
 d4est_geometry_t* d4est_geom,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq
)
{
  d4est_geometry_brick_attr_t* brick_attrs = d4est_geom->user;
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;
  
  double x0 = brick_attrs->X0;
  double x1 = brick_attrs->X1;
  double y0 = brick_attrs->Y0;
  double y1 = brick_attrs->Y1;

  /* transform element corners to [X0,X1] X [YO,Y1] X [Z0,Z1] topological space */
  amin = (x1-x0)*amin + x0;
  amax = (x1-x0)*amax + x0;
  bmin = (y1-y0)*bmin + y0;
  bmax = (y1-y0)*bmax + y0;
  
  if (deg == 1){
    double term1_boundary [4][4];
    term1_boundary[0][0] = (bmax - bmin)/2.;
    term1_boundary[0][1] = 0;
    term1_boundary[0][2] = (bmax - bmin)/2.;
    term1_boundary[0][3] = 0;
    term1_boundary[1][0] = 0;
    term1_boundary[1][1] = (-bmax + bmin)/2.;
    term1_boundary[1][2] = 0;
    term1_boundary[1][3] = (-bmax + bmin)/2.;
    term1_boundary[2][0] = (amax - amin)/2.;
    term1_boundary[2][1] = (amax - amin)/2.;
    term1_boundary[2][2] = 0;
    term1_boundary[2][3] = 0;
    term1_boundary[3][0] = 0;
    term1_boundary[3][1] = 0;
    term1_boundary[3][2] = (-amax + amin)/2.;
    term1_boundary[3][3] = (-amax + amin)/2.;
    return term1_boundary[face][i];
  }
  else {
    D4EST_ABORT("no other degrees supported atm\n");
  }
}



static double
testd4est_poisson_2_brick_dirichlet_term2
(
 int deg,
 int face,
 int i,
 double penalty,
 d4est_geometry_t* d4est_geom,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq
)
{
  d4est_geometry_brick_attr_t* brick_attrs = d4est_geom->user;
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;
  
  double x0 = brick_attrs->X0;
  double x1 = brick_attrs->X1;
  double y0 = brick_attrs->Y0;
  double y1 = brick_attrs->Y1;

  /* transform element corners to [X0,X1] X [YO,Y1] X [Z0,Z1] topological space */
  amin = (x1-x0)*amin + x0;
  amax = (x1-x0)*amax + x0;
  bmin = (y1-y0)*bmin + y0;
  bmax = (y1-y0)*bmax + y0;
  
  if (deg == 1){
    double term2_boundary [4][4];
    term2_boundary[0][0] = -((bmax - bmin)*(6*amin + 2*bmax + 4*bmin))/(12.*(amax - amin));
    term2_boundary[0][1] = ((bmax - bmin)*(2*amin + (2*(bmax + 2*bmin))/3.))/(4.*(amax - amin));
    term2_boundary[0][2] = -((bmax - bmin)*(6*amin + 4*bmax + 2*bmin))/(12.*(amax - amin));
    term2_boundary[0][3] = ((bmax - bmin)*(2*amin + (2*(2*bmax + bmin))/3.))/(4.*(amax - amin));
    term2_boundary[1][0] = ((bmax - bmin)*(6*amax + 2*bmax + 4*bmin))/(12.*(amax - amin));
    term2_boundary[1][1] = -((bmax - bmin)*(2*amax + (2*(bmax + 2*bmin))/3.))/(4.*(amax - amin));
    term2_boundary[1][2] = ((bmax - bmin)*(6*amax + 4*bmax + 2*bmin))/(12.*(amax - amin));
    term2_boundary[1][3] = -((bmax - bmin)*(2*amax + (2*(2*bmax + bmin))/3.))/(4.*(amax - amin));
    term2_boundary[2][0] = -((amax - amin)*(2*amax + 4*amin + 6*bmin))/(12.*(bmax - bmin));
    term2_boundary[2][1] = -((amax - amin)*(4*amax + 2*amin + 6*bmin))/(12.*(bmax - bmin));
    term2_boundary[2][2] = ((amax - amin)*((2*amax)/3. + (4*amin)/3. + 2*bmin))/(4.*(bmax - bmin));
    term2_boundary[2][3] = ((amax - amin)*((4*amax)/3. + (2*amin)/3. + 2*bmin))/(4.*(bmax - bmin));
    term2_boundary[3][0] = ((amax - amin)*(2*amax + 4*amin + 6*bmax))/(12.*(bmax - bmin));
    term2_boundary[3][1] = ((amax - amin)*(4*amax + 2*amin + 6*bmax))/(12.*(bmax - bmin));
    term2_boundary[3][2] = -((amax - amin)*((2*amax)/3. + (4*amin)/3. + 2*bmax))/(4.*(bmax - bmin));
    term2_boundary[3][3] = -((amax - amin)*((4*amax)/3. + (2*amin)/3. + 2*bmax))/(4.*(bmax - bmin));
    return term2_boundary[face][i];
  }
  else {
    D4EST_ABORT("no other degrees supported atm\n");
  }
}

static void
testd4est_poisson_2_brick_interface_old_style
(
 d4est_element_data_t** e_m,
 int faces_m,
 int f_m,
 int mortar_side_id_m,
 d4est_element_data_t** e_p,
 int faces_p,
 int f_p,
 int mortar_side_id_p,
 int* e_m_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_poisson_flux_interface_data_t* mortar_data,
 void* params
)
{
  d4est_quadrature_mortar_t* mortar_face_object = mortar_data->mortar_face_object;
  int total_volume_nodes_m_lobatto = 0;
  for (int i = 0; i < faces_m; i++){
    total_volume_nodes_m_lobatto += d4est_lgl_get_nodes((P4EST_DIM), e_m[i]->deg);
  }
  
  double* u_m_on_f_m_mortar_quad = mortar_data->u_m_on_f_m_mortar_quad;
  double* sj_on_f_m_mortar_quad = mortar_data->sj_on_f_m_mortar_quad;
  double* j_div_sj_on_f_m_mortar_quad = mortar_data->j_div_sj_on_f_m_mortar_quad;
  double* u_p_on_f_p_mortar_quad = mortar_data->u_p_on_f_p_mortar_quad;
  double* j_div_sj_on_f_p_mortar_quad = mortar_data->j_div_sj_on_f_p_mortar_quad;
  
  double* drst_dxyz_m_on_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_quad [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_quad [(P4EST_DIM)];
  double* n_sj_on_f_m_mortar_quad [(P4EST_DIM)];

  D4EST_COPY_DBYD_MAT(mortar_data->drst_dxyz_m_on_mortar_quad, drst_dxyz_m_on_mortar_quad);
  D4EST_COPY_DIM_VEC(mortar_data->dudx_m_on_f_m_mortar_quad, dudx_m_on_f_m_mortar_quad);
  D4EST_COPY_DIM_VEC(mortar_data->dudx_p_on_f_p_mortar_quad, dudx_p_on_f_p_mortar_quad);
  D4EST_COPY_DIM_VEC(mortar_data->n_sj_on_f_m_mortar_quad, n_sj_on_f_m_mortar_quad);
  
  int faces_mortar = mortar_data->faces_mortar;
  int total_side_nodes_m_lobatto = mortar_data->total_side_nodes_m_lobatto;
  int total_side_nodes_p_lobatto = mortar_data->total_side_nodes_p_lobatto;
  int total_nodes_mortar_lobatto = mortar_data->total_nodes_mortar_lobatto;
  int total_nodes_mortar_quad = mortar_data->total_nodes_mortar_quad;
  
  int* deg_mortar_quad = mortar_data->deg_mortar_quad;
  int* nodes_mortar_quad = mortar_data->nodes_mortar_quad;
  int* nodes_mortar_lobatto = mortar_data->nodes_mortar_lobatto;
  int* deg_mortar_lobatto = mortar_data->deg_mortar_lobatto;
  int* face_nodes_m_lobatto = mortar_data->face_nodes_m_lobatto;
  int* face_nodes_p_lobatto = mortar_data->face_nodes_p_lobatto;
  int* deg_m_lobatto = mortar_data->deg_m_lobatto;
  int* deg_p_lobatto = mortar_data->deg_p_lobatto;
  
  d4est_poisson_flux_sipg_params_t* ip_flux_params = (d4est_poisson_flux_sipg_params_t*) params;
  double sipg_kronbichler_flux_penalty_prefactor = ip_flux_params->sipg_penalty_prefactor;
  penalty_calc_t sipg_kronbichler_flux_penalty_calculate_fcn = ip_flux_params->sipg_penalty_fcn;

  double* du_m_on_f_m [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(du_m_on_f_m, total_side_nodes_m_lobatto);
  double* du_m_on_f_m_mortar [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(du_m_on_f_m_mortar, total_nodes_mortar_lobatto);
  double* du_p_on_f_p [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(du_p_on_f_p, total_side_nodes_p_lobatto);
  double* du_p_on_f_p_mortar [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(du_p_on_f_p_mortar, total_nodes_mortar_lobatto);
  double* term1_old_style_mortar [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(term1_old_style_mortar, total_nodes_mortar_lobatto);
  double* term1_old_style_side_m [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(term1_old_style_side_m, total_side_nodes_m_lobatto);
  double* M_term1_old_style_side_m [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(M_term1_old_style_side_m, total_side_nodes_m_lobatto);
  double* M_term1_old_style_mortar [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(M_term1_old_style_mortar, total_nodes_mortar_lobatto);
  double* proj_M_term1_old_style_mortar [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(proj_M_term1_old_style_mortar, total_nodes_mortar_lobatto);
  double* lifted_proj_M_term1_old_style_mortar [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(lifted_proj_M_term1_old_style_mortar, total_volume_nodes_m_lobatto);
  double* lifted_M_term1_old_style_side_m [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(lifted_M_term1_old_style_side_m, total_volume_nodes_m_lobatto);
  
  for (int dir = 0; dir < (P4EST_DIM); dir++){
    /* compute the (-)-u-derivative and project on the (-)-side faces and project q onto the (-)-side faces */
    int stride = 0;
    for (int f = 0; f < faces_m; f++){
      d4est_operators_apply_slicer
        (
         d4est_ops,
         e_m[f]->dudr_elem[dir],
         (P4EST_DIM),
         f_m,
         e_m[f]->deg,
         &du_m_on_f_m[dir][stride]
        );
      double h = e_m[f]->dq/(double)P4EST_ROOT_LEN;
      d4est_linalg_vec_scale(2./h, &du_m_on_f_m[dir][stride], face_nodes_m_lobatto[f]);
      stride += face_nodes_m_lobatto[f];
    }

    /* compute the (+)-u-derivative and project on the (+)-side faces and project q onto the (+)-side faces */
    stride = 0;
    for (int f = 0; f < faces_p; f++){
      d4est_operators_apply_slicer
        (
         d4est_ops,
         e_p[f]->dudr_elem[dir],
         (P4EST_DIM),
         f_p,
         e_p[f]->deg,
         &du_p_on_f_p[dir][stride]
        );
      double h = e_p[f]->dq/(double)P4EST_ROOT_LEN;
      d4est_linalg_vec_scale(2./h, &du_p_on_f_p[dir][stride], face_nodes_p_lobatto[f]);
      stride += face_nodes_p_lobatto[f];
    }

    /* project the derivatives from (-) and (+) sides onto the mortar space */
    d4est_mortars_project_side_onto_mortar_space
      (
       d4est_ops,
       du_m_on_f_m[dir],
       faces_m,
       deg_m_lobatto,
       du_m_on_f_m_mortar[dir],
       faces_mortar,
       deg_mortar_lobatto
      );

    d4est_mortars_project_side_onto_mortar_space
      (
       d4est_ops,
       du_p_on_f_p[dir],
       faces_p,
       deg_p_lobatto,
       du_p_on_f_p_mortar[dir],
       faces_mortar,
       deg_mortar_lobatto
      );

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){
      for (int k = 0; k < nodes_mortar_lobatto[f]; k++){
        int ks = k + stride;
        /* term1_old_style_mortar[dir][ks] = .5*(du_p_on_f_p_mortar[dir][ks] + du_m_on_f_m_mortar[dir][ks]); */
        term1_old_style_mortar[dir][ks] = 5.;
        /* printf (".5*(du_p_on_f_p_mortar[%d][%d] + du_m_on_f_m_mortar[][] = %.25f\n", dir, ks, */
                /* .5*(du_p_on_f_p_mortar[dir][ks] + du_m_on_f_m_mortar[dir][ks]) */
               /* ); */
      }
      stride += nodes_mortar_lobatto[f];
    }

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){
      for (int k = 0; k < nodes_mortar_quad[f]; k++){
        int ks = k + stride;
        double temp = .5*(mortar_data->dudx_m_on_f_m_mortar_quad[dir][ks] + mortar_data->dudx_p_on_f_p_mortar_quad[dir][ks]);
        /* printf ("FROM MORTAR DATA .5*(du_p_on_f_p_mortar[%d][%d] + du_m_on_f_m_mortar[][] = %.25f\n", dir, ks, temp); */
      }
      stride += nodes_mortar_lobatto[f];
    }

    double* term1_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
    double* ones_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
    double* VT_w_term1_mortar_lobatto = P4EST_ALLOC(double, total_nodes_mortar_lobatto);
    d4est_util_fill_array(ones_mortar_quad, 1., total_nodes_mortar_quad);
  stride = 0;
  int stride_lobatto = 0;
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_quad[f]; k++){
      int ks = k + stride;

      term1_mortar_quad[ks] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
      term1_mortar_quad[ks] += 5.*n_sj_on_f_m_mortar_quad[d][ks];
      }
      /* term1_mortar_quad[ks] = 0.; */
      /* for (int d = 0; d < (P4EST_DIM); d++){ */
        /* term1_mortar_quad[ks] += -1.*n_sj_on_f_m_mortar_quad[d][ks] */
                                 /* *.5*(dudx_p_on_f_p_mortar_quad[d][ks] + dudx_m_on_f_m_mortar_quad[d][ks]); */
      /* } */
    }


    d4est_quadrature_apply_galerkin_integral
      (
       d4est_ops,
       d4est_geom,
       d4est_quad,
       mortar_face_object,
       QUAD_OBJECT_MORTAR,
       QUAD_INTEGRAND_UNKNOWN,
       &term1_mortar_quad[stride],
       deg_mortar_lobatto[f],
       &ones_mortar_quad[stride],
       deg_mortar_quad[f],
       &VT_w_term1_mortar_lobatto[stride_lobatto]
      );

   for (int k = 0; k < nodes_mortar_lobatto[f]; k++){
      int ks = k + stride;
        printf ("USING QUADRATURE VT_w_term1_mortar_lobatto[%d] = %.25f\n", ks, VT_w_term1_mortar_lobatto[ks]);
    }
    
    
    stride += nodes_mortar_quad[f];
    stride_lobatto += nodes_mortar_lobatto[f];
  }
    
    P4EST_FREE(term1_mortar_quad);
    P4EST_FREE(ones_mortar_quad);
    P4EST_FREE(VT_w_term1_mortar_lobatto);
    
    stride = 0;
    for (int f = 0; f < faces_m; f++){
      d4est_operators_apply_mij(
                                d4est_ops,
                                &term1_old_style_mortar[dir][stride],
                                (P4EST_DIM) - 1,
                                deg_mortar_lobatto[f],
                                &M_term1_old_style_mortar[dir][stride]
                              );

      double h = (double)e_m[f]->dq/(double)P4EST_ROOT_LEN;
      double surface_jacobian = d4est_util_dbl_pow_int(.5*h, (P4EST_DIM) - 1);
      double n [(P4EST_DIM)];
      d4est_reference_get_normal(f_m, (P4EST_DIM), &n[0]);
      d4est_linalg_vec_scale(-n[dir]*surface_jacobian, &M_term1_old_style_mortar[dir][stride], nodes_mortar_lobatto[f]);


      for (int k = 0; k < nodes_mortar_lobatto[f]; k++){
        int ks = k + stride;
        double temp = 0.;
        for (int d = 0; d < (P4EST_DIM); d++)
          temp += M_term1_old_style_mortar[d][ks];
        printf ("USING MIJ VT_w_term1_mortar_lobatto[%d] = %.25f\n", ks, temp);
      }
      stride += nodes_mortar_lobatto[f];
    }
    
   d4est_mortars_project_mass_mortar_onto_side
      (
       d4est_ops,
       M_term1_old_style_mortar[dir],
       faces_mortar,
       deg_mortar_lobatto,
       proj_M_term1_old_style_mortar[dir],
       faces_m,
       deg_m_lobatto
      );

   
   d4est_mortars_project_mortar_onto_side
      (
       d4est_ops,
       term1_old_style_mortar[dir],
       faces_mortar,
       deg_mortar_lobatto,
       term1_old_style_side_m[dir],
       faces_m,
       deg_m_lobatto
      );

   stride = 0;
   int stride_volume = 0;

   for (int f = 0; f < faces_m; f++){
     d4est_operators_apply_mij(
                               d4est_ops,
                               &term1_old_style_side_m[dir][stride],
                               (P4EST_DIM) - 1,
                               e_m[f]->deg,
                               &M_term1_old_style_side_m[dir][stride]
                              );

     double h = (double)e_m[f]->dq/(double)P4EST_ROOT_LEN;
     double surface_jacobian = d4est_util_dbl_pow_int(.5*h, (P4EST_DIM) - 1);
     double n [(P4EST_DIM)];
     d4est_reference_get_normal(f_m, (P4EST_DIM), &n[0]);
     d4est_linalg_vec_scale(-n[dir]*surface_jacobian,
                            &M_term1_old_style_side_m[dir][stride],
                            face_nodes_m_lobatto[f]);
      
     d4est_operators_apply_lift(d4est_ops,
                                &M_term1_old_style_side_m[dir][stride],
                                (P4EST_DIM),
                                e_m[f]->deg,
                                f_m, &lifted_M_term1_old_style_side_m[dir][stride_volume]);

     d4est_operators_apply_lift(d4est_ops, &proj_M_term1_old_style_mortar[dir][stride], (P4EST_DIM), e_m[f]->deg, f_m, &lifted_proj_M_term1_old_style_mortar[dir][stride_volume]);
     
     stride += face_nodes_m_lobatto[f];
     stride_volume += d4est_lgl_get_nodes((P4EST_DIM), e_m[f]->deg);
   }
  }

  testd4est_poisson_2_brick_data_t* data = params;
  int volume_stride = 0;
  for (int f = 0; f < faces_m; f++){
    int volume_nodes_m = d4est_lgl_get_nodes((P4EST_DIM), e_m[f]->deg);
    for (int i = 0; i < volume_nodes_m; i++){
      double term1_final = 0.;
      for (int d = 0; d < (P4EST_DIM); d++)
        term1_final += lifted_M_term1_old_style_side_m[d][i+volume_stride];

      double term1_final_2 = 0.;
      for (int d = 0; d < (P4EST_DIM); d++)
        term1_final_2 += lifted_proj_M_term1_old_style_mortar[d][i+volume_stride];
      
      
      double term1_check = testd4est_poisson_2_brick_dirichlet_term1
                           (
                            e_m[f]->deg,
                            f_m,
                            i,
                            data->sipg_params->sipg_penalty_prefactor,
                            d4est_geom,
                            e_m[f]->q,
                            e_m[f]->dq
                           );

      /* printf("i, term1_final, term1_final_2, term1_check = %d, %f, %f, %f\n", i, term1_final, term1_final_2, term1_check); */
    }
    volume_stride += volume_nodes_m;
  }


  D4EST_FREE_DIM_VEC(du_m_on_f_m);
  D4EST_FREE_DIM_VEC(du_m_on_f_m_mortar);
  D4EST_FREE_DIM_VEC(du_p_on_f_p);
  D4EST_FREE_DIM_VEC(du_p_on_f_p_mortar);
  
  D4EST_FREE_DIM_VEC(term1_old_style_mortar);
  D4EST_FREE_DIM_VEC(term1_old_style_side_m);
  D4EST_FREE_DIM_VEC(M_term1_old_style_side_m);
  D4EST_FREE_DIM_VEC(lifted_M_term1_old_style_side_m);

  D4EST_FREE_DIM_VEC(M_term1_old_style_mortar);
  D4EST_FREE_DIM_VEC(proj_M_term1_old_style_mortar);
  D4EST_FREE_DIM_VEC(lifted_proj_M_term1_old_style_mortar);
}

static void
testd4est_poisson_2_brick_interface
(
 d4est_element_data_t** e_m,
 int faces_m,
 int f_m,
 int mortar_side_id_m,
 d4est_element_data_t** e_p,
 int faces_p,
 int f_p,
 int mortar_side_id_p,
 int* e_m_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_poisson_flux_interface_data_t* mortar_data,
 void* params
)
{
  testd4est_poisson_2_brick_data_t* data = params;
  int faces_mortar = mortar_data->faces_mortar;
  int total_side_nodes_m_lobatto = mortar_data->total_side_nodes_m_lobatto;
  int total_nodes_mortar_lobatto = mortar_data->total_nodes_mortar_lobatto;
  int total_nodes_mortar_quad = mortar_data->total_nodes_mortar_quad;
  
  int total_volume_nodes_m_lobatto = 0;
  for (int i = 0; i < faces_m; i++){
    total_volume_nodes_m_lobatto += d4est_lgl_get_nodes((P4EST_DIM), e_m[i]->deg);
  }

  double* lifted_proj_VT_w_term1_mortar_lobatto = P4EST_ALLOC(double, total_volume_nodes_m_lobatto);
  double* proj_VT_w_term1_mortar_lobatto = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
  double* VT_w_term1_mortar_lobatto = P4EST_ALLOC(double, total_nodes_mortar_lobatto);
  double* term1_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* DT_lifted_proj_VT_w_term2_mortar_lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(DT_lifted_proj_VT_w_term2_mortar_lobatto, total_volume_nodes_m_lobatto);
  double* lifted_proj_VT_w_term2_mortar_lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(lifted_proj_VT_w_term2_mortar_lobatto, total_volume_nodes_m_lobatto);
  double* proj_VT_w_term2_mortar_lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(proj_VT_w_term2_mortar_lobatto, total_side_nodes_m_lobatto);
  double* VT_w_term2_mortar_lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(VT_w_term2_mortar_lobatto, total_nodes_mortar_lobatto);
  double* term2_mortar_quad [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(term2_mortar_quad, total_nodes_mortar_quad);
  double* lifted_proj_VT_w_term3_mortar_lobatto = P4EST_ALLOC(double, total_volume_nodes_m_lobatto);
  double* proj_VT_w_term3_mortar_lobatto = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
  double* VT_w_term3_mortar_lobatto = P4EST_ALLOC(double, total_nodes_mortar_lobatto);
  double* term3_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  
  d4est_poisson_flux_sipg_interface_aux
    (
     e_m,
     faces_m,
     f_m,
     mortar_side_id_m,
     e_p,
     faces_p,
     f_p,
     mortar_side_id_p,
     e_m_is_ghost,
     orientation,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     mortar_data,
     data->sipg_params,
     lifted_proj_VT_w_term1_mortar_lobatto,
     proj_VT_w_term1_mortar_lobatto,
     VT_w_term1_mortar_lobatto,
     term1_mortar_quad,
     DT_lifted_proj_VT_w_term2_mortar_lobatto,
     lifted_proj_VT_w_term2_mortar_lobatto,
     proj_VT_w_term2_mortar_lobatto,
     VT_w_term2_mortar_lobatto,
     term2_mortar_quad,
     lifted_proj_VT_w_term3_mortar_lobatto,
     proj_VT_w_term3_mortar_lobatto,
     VT_w_term3_mortar_lobatto,
     term3_mortar_quad
    );

  double* lifted_proj_VT_w_term1_mortar_lobatto_check = P4EST_ALLOC(double, total_volume_nodes_m_lobatto);
  
  int stride = 0;
  for (int f = 0; f < faces_m; f++){
    int volume_nodes_m = d4est_lgl_get_nodes((P4EST_DIM), e_m[f]->deg);
    if (e_m_is_ghost[f] == 0){
      for (int i = 0; i < volume_nodes_m; i++){

        lifted_proj_VT_w_term1_mortar_lobatto_check[i] = testd4est_poisson_2_brick_dirichlet_term1
                                             (
                                              e_m[f]->deg,
                                              f_m,
                                              i,
                                              data->sipg_params->sipg_penalty_prefactor,
                                              d4est_geom,
                                              e_m[f]->q,
                                              e_m[f]->dq
                                             );

        printf("lifted_proj_VT_w_term1_mortar_lobatto[%d] = %.25f, lifted_proj_VT_w_term1_mortar_lobatto[%d] = %.25f\n",i, lifted_proj_VT_w_term1_mortar_lobatto[i + stride], i, lifted_proj_VT_w_term1_mortar_lobatto_check[i]);
        
        data->interface_term_1_err += fabs(lifted_proj_VT_w_term1_mortar_lobatto_check[i] - lifted_proj_VT_w_term1_mortar_lobatto[i + stride]);
      }
    }
    stride += volume_nodes_m;
  }

  P4EST_FREE(lifted_proj_VT_w_term1_mortar_lobatto_check);

  P4EST_FREE(lifted_proj_VT_w_term1_mortar_lobatto);
  P4EST_FREE(proj_VT_w_term1_mortar_lobatto);
  P4EST_FREE(VT_w_term1_mortar_lobatto);
  P4EST_FREE(term1_mortar_quad);

  D4EST_FREE_DIM_VEC(DT_lifted_proj_VT_w_term2_mortar_lobatto);
  D4EST_FREE_DIM_VEC(lifted_proj_VT_w_term2_mortar_lobatto);
  D4EST_FREE_DIM_VEC(proj_VT_w_term2_mortar_lobatto);
  D4EST_FREE_DIM_VEC(VT_w_term2_mortar_lobatto);
  D4EST_FREE_DIM_VEC(term2_mortar_quad);

  P4EST_FREE(lifted_proj_VT_w_term3_mortar_lobatto);
  P4EST_FREE(proj_VT_w_term3_mortar_lobatto);
  P4EST_FREE(VT_w_term3_mortar_lobatto);
  P4EST_FREE(term3_mortar_quad);
}


static void
testd4est_poisson_2_brick_dirichlet
(
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_xyz_fcn_t boundary_condition,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_poisson_flux_boundary_data_t* boundary_data,
 void* params
)
{
  testd4est_poisson_2_brick_data_t* data = params;
  
  int face_nodes_m_lobatto = boundary_data->face_nodes_m_lobatto;
  int face_nodes_m_quad = boundary_data->face_nodes_m_quad;
  int volume_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM), e_m->deg);
  
  double* term1_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  double* VT_w_term1_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* lifted_VT_w_term1_lobatto = P4EST_ALLOC(double, volume_nodes_m_lobatto);

  double* term2_quad [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(term2_quad, face_nodes_m_quad);
  double* VT_w_term2_lobatto [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(VT_w_term2_lobatto, face_nodes_m_lobatto);
  double* lifted_VT_w_term2_lobatto [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(lifted_VT_w_term2_lobatto, volume_nodes_m_lobatto);
  double* DT_lifted_VT_w_term2_lobatto [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(DT_lifted_VT_w_term2_lobatto, volume_nodes_m_lobatto);
 
  double* term3_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  double* VT_w_term3_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* lifted_VT_w_term3_lobatto = P4EST_ALLOC(double, volume_nodes_m_lobatto);
  
  double* sigma = P4EST_ALLOC(double, face_nodes_m_quad);

  double* u_at_bndry_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* u_at_bndry_lobatto_to_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  
  d4est_poisson_flux_sipg_dirichlet_aux
    (
     e_m,
     f_m,
     mortar_side_id_m,
     boundary_condition,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     boundary_data,
     data->sipg_params,
     term1_quad,
     VT_w_term1_lobatto,
     lifted_VT_w_term1_lobatto,
     term2_quad,
     VT_w_term2_lobatto,
     lifted_VT_w_term2_lobatto,
     DT_lifted_VT_w_term2_lobatto,
     term3_quad,
     VT_w_term3_lobatto,
     lifted_VT_w_term3_lobatto,
     sigma,
     u_at_bndry_lobatto,
     u_at_bndry_lobatto_to_quad
    );

  double* lifted_VT_w_term3_lobatto_check = P4EST_ALLOC(double, volume_nodes_m_lobatto);
  double* lifted_VT_w_term1_lobatto_check = P4EST_ALLOC(double, volume_nodes_m_lobatto);
  double* DT_lifted_VT_w_term2_lobatto_check = P4EST_ALLOC(double, volume_nodes_m_lobatto);
 
  d4est_poisson_flux_sipg_params_t* ip_flux_params = (d4est_poisson_flux_sipg_params_t*) data->sipg_params;
  double sipg_kronbichler_flux_penalty_prefactor = ip_flux_params->sipg_penalty_prefactor;

  /* int volume_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM), e_m->deg); */
  for (int i = 0; i < volume_nodes_m_lobatto; i++){
    lifted_VT_w_term3_lobatto_check[i] = testd4est_poisson_2_brick_dirichlet_term3
                                         (
                                          e_m->deg,
                                          f_m,
                                          i,
                                          sipg_kronbichler_flux_penalty_prefactor,
                                          d4est_geom,
                                          e_m->q,
                                          e_m->dq
                                         );

    data->boundary_term_3_err += fabs(lifted_VT_w_term3_lobatto_check[i] - lifted_VT_w_term3_lobatto[i]);
  }


  /* printf("e_m->id, f_m = %d, %d\n", e_m->id, f_m);   */
  for (int i = 0; i < volume_nodes_m_lobatto; i++){
    lifted_VT_w_term1_lobatto_check[i] = testd4est_poisson_2_brick_dirichlet_term1
                                         (
                                          e_m->deg,
                                          f_m,
                                          i,
                                          sipg_kronbichler_flux_penalty_prefactor,
                                          d4est_geom,
                                          e_m->q,
                                          e_m->dq
                                         );



    data->boundary_term_1_err += fabs(lifted_VT_w_term1_lobatto_check[i] - lifted_VT_w_term1_lobatto[i]);
  }

  /* printf(" DT_lifted_VT_w_term2_lobatto_check, DT_lifted_VT_w_term2_lobatto = \n"); */
  for (int i = 0; i < volume_nodes_m_lobatto; i++){
    DT_lifted_VT_w_term2_lobatto_check[i] = testd4est_poisson_2_brick_dirichlet_term2
                                         (
                                          e_m->deg,
                                          f_m,
                                          i,
                                          sipg_kronbichler_flux_penalty_prefactor,
                                          d4est_geom,
                                          e_m->q,
                                          e_m->dq
                                         );


    double term2_real = 0.;
    for (int d = 0; d < (P4EST_DIM); d++){
      term2_real += DT_lifted_VT_w_term2_lobatto[d][i];
    }
    

    data->boundary_term_2_err += fabs(DT_lifted_VT_w_term2_lobatto_check[i] - term2_real);
    /* if (data->boundary_term_2_err > D4EST_REAL_EPS){ */
      /* printf("%.25f, %.25f\n",DT_lifted_VT_w_term2_lobatto_check[i], term2_real); */
    /* } */
  }

  /* DEBUG_PRINT_2ARR_DBL(lifted_VT_w_term3_lobatto, lifted_VT_w_term3_lobatto_check, volume_nodes_m_lobatto); */

  /* DEBUG_PRINT_2ARR_DBL(lifted_VT_w_term1_lobatto, lifted_VT_w_term1_lobatto_check, volume_nodes_m_lobatto); */

  
  P4EST_FREE(lifted_VT_w_term3_lobatto_check);
  P4EST_FREE(DT_lifted_VT_w_term2_lobatto_check);
  P4EST_FREE(lifted_VT_w_term1_lobatto_check);

  P4EST_FREE(u_at_bndry_lobatto);
  P4EST_FREE(u_at_bndry_lobatto_to_quad);
  P4EST_FREE(sigma);
  P4EST_FREE(term1_quad);
  P4EST_FREE(VT_w_term1_lobatto);
  P4EST_FREE(lifted_VT_w_term1_lobatto);
  D4EST_FREE_DIM_VEC(term2_quad);
  D4EST_FREE_DIM_VEC(VT_w_term2_lobatto);
  D4EST_FREE_DIM_VEC(lifted_VT_w_term2_lobatto);
  D4EST_FREE_DIM_VEC(DT_lifted_VT_w_term2_lobatto);
  P4EST_FREE(term3_quad);
  P4EST_FREE(VT_w_term3_lobatto);
  P4EST_FREE(lifted_VT_w_term3_lobatto);
}

void
testd4est_poisson_2_brick_on_interfaces
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_elliptic_data_t* prob_vecs,
 double err_tol
)
{
  d4est_poisson_flux_data_t* d4est_poisson_flux_data = P4EST_ALLOC(d4est_poisson_flux_data_t,1);
  testd4est_poisson_2_brick_data_t* data = P4EST_ALLOC(testd4est_poisson_2_brick_data_t, 1);
  data->boundary_term_3_err = 0.;
  data->boundary_term_2_err = 0.;
  data->boundary_term_1_err = 0.;
  data->interface_term_1_err = 0.;

  d4est_poisson_flux_sipg_params_t* d4est_poisson_flux_sipg_params = P4EST_ALLOC(d4est_poisson_flux_sipg_params_t, 1);
  d4est_poisson_flux_sipg_params_input(p4est, "flux", "testd4est_poisson_2_brick.input", d4est_poisson_flux_sipg_params);

  data->sipg_params = d4est_poisson_flux_sipg_params;
  
  d4est_poisson_flux_data->user = data;
  d4est_poisson_flux_data->interface_fcn = testd4est_poisson_2_brick_interface_old_style;
  d4est_poisson_flux_data->boundary_fcn = testd4est_poisson_2_brick_dirichlet;
  d4est_poisson_flux_data->boundary_condition = zero_fcn;
  d4est_poisson_flux_data->destroy = NULL;

  d4est_poisson_flux_init_element_data(p4est, d4est_ops, prob_vecs->u, prob_vecs->Au);
  
  d4est_mortars_fcn_ptrs_t flux_fcns = d4est_poisson_flux_fetch_fcns(d4est_poisson_flux_data);
  d4est_mortars_compute_flux_on_local_elements
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &flux_fcns,
     EXCHANGE_GHOST_DATA
    );

  printf("data->boundary_term_3_err = %.25f\n", data->boundary_term_3_err);
  printf("data->boundary_term_2_err = %.25f\n", data->boundary_term_2_err);
  printf("data->boundary_term_1_err = %.25f\n", data->boundary_term_1_err);
  printf("data->interface_term_1_err = %.25f\n", data->interface_term_1_err);
  int err = (data->boundary_term_3_err < err_tol);
  
  P4EST_FREE(d4est_poisson_flux_data);
  P4EST_FREE(d4est_poisson_flux_sipg_params);
  P4EST_FREE(data);
  D4EST_ASSERT(err);
}
  
int main(int argc, char *argv[])
{
  sc_MPI_Comm mpicomm;
  PetscInitialize(&argc,&argv,(char*)0,NULL);
  mpicomm = PETSC_COMM_WORLD;
  
  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  p4est_init(NULL, SC_LP_ERROR);
  
  d4est_geometry_t* d4est_geom = P4EST_ALLOC(d4est_geometry_t, 1);
  d4est_geom->geom_type = GEOM_BRICK;
  d4est_geom->DX_compute_method = GEOM_COMPUTE_ANALYTIC;
  d4est_geom->JAC_compute_method = GEOM_COMPUTE_NUMERICAL;

  d4est_geometry_brick_new(proc_rank, "testd4est_poisson_2_brick.input", "geometry", "[Geometry]:", d4est_geom);
    
  p4est_t* p4est = problem_build_p4est
                   (
                    mpicomm,
                    d4est_geom->p4est_conn,
                    -1,
                    1,
                    1
                   );
  
  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_data_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t, 1);
  d4est_quad->quad_type = QUAD_TYPE_GAUSS_LEGENDRE;
  d4est_quadrature_legendre_new(d4est_quad, d4est_geom,"");

  
  d4est_poisson_flux_data_t* flux_data = d4est_poisson_flux_new(p4est, "testd4est_poisson_2_brick.input", zero_fcn, NULL);
  d4est_poisson_flux_data_t* flux_data_with_bc = d4est_poisson_flux_new(p4est, "testd4est_poisson_2_brick.input", poly_vec_fcn, NULL);
      
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);

  d4est_amr_t* d4est_amr = d4est_amr_init(
                                          p4est,
                                          "testd4est_poisson_2_brick.input",
                                          NULL
  );

  int local_nodes = d4est_mesh_update
                    (
                     p4est,
                     ghost,
                     ghost_data,
                     d4est_ops,
                     d4est_geom,
                     d4est_quad,
                     geometric_factors,
                     INITIALIZE_QUADRATURE_DATA,
                     INITIALIZE_GEOMETRY_DATA,
                     INITIALIZE_GEOMETRY_ALIASES,
                     problem_set_degrees_init,
                     NULL
                    );
  
  double* poly_vec = P4EST_ALLOC(double, local_nodes);
  int same = 1;
  int same2 = 1;
  
  for (int level = 0; level < d4est_amr->num_of_amr_steps; ++level){

    local_nodes = d4est_mesh_update
                  (
                   p4est,
                   ghost,
                   ghost_data,
                   d4est_ops,
                   d4est_geom,
                   d4est_quad,
                   geometric_factors,
                   INITIALIZE_QUADRATURE_DATA,
                   INITIALIZE_GEOMETRY_DATA,
                   INITIALIZE_GEOMETRY_ALIASES,
                   problem_set_degrees_amr,
                   NULL
                  );

    
    printf("level = %d, elements = %d, nodes = %d\n", level, p4est->local_num_quadrants, local_nodes);

    if (level == 0){
      d4est_mesh_init_field(p4est, poly_vec, poly_vec_fcn, d4est_ops, d4est_geom, NULL);
    }
    else {

      double* poly_vec_compare = P4EST_ALLOC(double, local_nodes);
      d4est_mesh_init_field(p4est, poly_vec_compare, poly_vec_fcn, d4est_ops, d4est_geom, NULL);
      same = d4est_util_compare_vecs(poly_vec, poly_vec_compare, local_nodes, D4EST_REAL_EPS);
      if (!same){
        DEBUG_PRINT_2ARR_DBL(poly_vec, poly_vec_compare, local_nodes);
        D4EST_ABORT("poly_vec and poly_vec_compare are not the same");
      }
      P4EST_FREE(poly_vec_compare);
      
      double* Apoly_vec = P4EST_ALLOC(double, local_nodes);
      d4est_elliptic_data_t elliptic_data;
      elliptic_data.u = poly_vec;
      elliptic_data.Au = Apoly_vec;
      elliptic_data.local_nodes = local_nodes;

      testd4est_poisson_2_brick_on_interfaces
        (
         p4est,
         ghost,
         ghost_data,
         d4est_ops,
         d4est_geom,
         d4est_quad,
         &elliptic_data,
         (D4EST_REAL_EPS)*100
        );
      
      P4EST_FREE(Apoly_vec);
    }

    
    d4est_amr_step
      (
       p4est,
       &ghost,
       &ghost_data,
       d4est_ops,
       d4est_amr,
       &poly_vec,
       NULL
      );
    
  }

  P4EST_FREE(poly_vec);

  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }
    
  d4est_poisson_flux_destroy(flux_data);
  d4est_poisson_flux_destroy(flux_data_with_bc);
  d4est_mesh_geometry_storage_destroy(geometric_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_amr_destroy(d4est_amr);
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);
  
  PetscFinalize();

}
