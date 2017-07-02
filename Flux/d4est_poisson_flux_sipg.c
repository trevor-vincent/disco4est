#include <util.h>
#include <d4est_operators.h>
#include <d4est_element_data.h>
#include <d4est_linalg.h>
#include <d4est_poisson_flux_sipg.h>
#include <d4est_mortars.h>
#include <d4est_quadrature_compactified.h>
#include <d4est_quadrature_lobatto.h>

/* #define NASTY_DEBUG */
#define NASTY_DEBUG2
#define NASTY_TRUMP_DEBUG

static void
d4est_poisson_flux_sipg_dirichlet
(
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_grid_fcn_t bndry_fcn,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* params
)
{
  ip_flux_params_t* ip_flux_params = (ip_flux_params_t*) params;
  double sipg_kronbichler_flux_penalty_prefactor = ip_flux_params->ip_flux_penalty_prefactor;
  penalty_calc_t sipg_kronbichler_flux_penalty_calculate_fcn = ip_flux_params->ip_flux_penalty_calculate_fcn;

  
  d4est_grid_fcn_t u_at_bndry = bndry_fcn;
  int face_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_quad = d4est_lgl_get_nodes((P4EST_DIM) - 1, e_m->deg_quad);

  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* u_m_on_f_m_quad = P4EST_ALLOC(double, face_nodes_m_quad);

  double* dudr_m_on_f_m [(P4EST_DIM)];
  double* dudr_m_on_f_m_quad [(P4EST_DIM)];
  double* dudx_m_on_f_m_quad [(P4EST_DIM)];
  
  double* u_m_on_f_m_min_u_at_bndry_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  double* u_m_on_f_m_min_u_at_bndry_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* sj_on_f_m_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  double* n_on_f_m_quad [(P4EST_DIM)];
  double* n_sj_on_f_m_quad [(P4EST_DIM)];

  double* drst_dxyz_quad [(P4EST_DIM)][(P4EST_DIM)]; D4EST_ALLOC_DBYD_MAT(drst_dxyz_quad, face_nodes_m_quad);
  
  for (int d = 0; d < (P4EST_DIM); d++) {
    dudr_m_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m_lobatto);
    dudr_m_on_f_m_quad[d] = P4EST_ALLOC(double, face_nodes_m_quad);
    dudx_m_on_f_m_quad[d] = P4EST_ALLOC(double, face_nodes_m_quad);
    n_on_f_m_quad[d] = P4EST_ALLOC(double, face_nodes_m_quad);
    n_sj_on_f_m_quad[d] = P4EST_ALLOC(double, face_nodes_m_quad);
  }

  double* u_at_bndry_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* u_at_bndry_quad = P4EST_ALLOC(double, face_nodes_m_quad);


  int volume_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM), e_m->deg);
  
  double* ones_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  d4est_linalg_fill_vec(ones_quad, 1., face_nodes_m_quad);
  
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

  double* J_div_SJ_quad = NULL;
  if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ || ip_flux_params->ip_flux_h_calc == H_EQ_J_DIV_SJ_MIN){
    J_div_SJ_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  }

  d4est_mortars_compute_geometric_data_on_mortar
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     QUAD_INTEGRAND_UNKNOWN,
     e_m->tree,
     e_m->q,
     e_m->dq,
     mortar_side_id_m,
     1,
     1,
     &e_m->deg_quad,
     f_m,
     drst_dxyz_quad,
     sj_on_f_m_quad,
     n_on_f_m_quad,
     n_sj_on_f_m_quad,
     J_div_SJ_quad,
     COMPUTE_NORMAL_USING_JACOBIAN
    );
  
  
  double* sigma = P4EST_ALLOC(double, face_nodes_m_quad);
  double h, h_min;
  if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ_MIN){
    h_min = util_min_dbl_array(J_div_SJ_quad, face_nodes_m_quad);
  }
    
  for (int i = 0; i < face_nodes_m_quad; i++){
    if (ip_flux_params->ip_flux_h_calc == H_EQ_VOLUME_DIV_AREA){
      mpi_abort("H_EQ_VOLUME_DIV_AREA no longer supported, will be completely deprecated soon\n");
      /* h = (e_m->volume/e_m->surface_area[f_m]); */
    }
    else if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ){
      h = J_div_SJ_quad[i];
    }
    else if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ_MIN){
      h = h_min;
    }
    sigma[i] = sipg_kronbichler_flux_penalty_calculate_fcn
               (
                e_m->deg,
                h,
                e_m->deg,
                h,
                sipg_kronbichler_flux_penalty_prefactor
               );

#ifdef NASTY_TRUMP_DEBUG
    printf("sigma[%d] = %f, h = %f, e_m->deg = %d\n", i, sigma[i], h, e_m->deg);
#endif

    
  }

  if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ || ip_flux_params->ip_flux_h_calc == H_EQ_J_DIV_SJ_MIN){
    P4EST_FREE(J_div_SJ_quad);
  }

  
  d4est_quadrature_mortar_t face_object;
  face_object.dq = e_m->dq;
  face_object.tree = e_m->tree;
  face_object.face = f_m;
  face_object.mortar_side_id = mortar_side_id_m;
  face_object.mortar_subface_id = 0;
  
  face_object.q[0] = e_m->q[0];
  face_object.q[1] = e_m->q[1];
#if (P4EST_DIM)==3
  face_object
    .q[2] = e_m->q[2];
#endif

  
  for (int d = 0; d < (P4EST_DIM); d++){
    
    d4est_operators_apply_slicer
      (
       d4est_ops,
       e_m->dudr_elem[d],
       (P4EST_DIM),
       f_m,
       e_m->deg,
       dudr_m_on_f_m[d]
      );

    d4est_quadrature_interpolate
      (
       d4est_ops,
       d4est_quad,
       d4est_geom,
       &face_object,
       QUAD_OBJECT_MORTAR,
       QUAD_INTEGRAND_UNKNOWN,
       dudr_m_on_f_m[d],
       e_m->deg,
       dudr_m_on_f_m_quad[d],
       e_m->deg_quad
      );
    
   
  }


  
  d4est_operators_apply_slicer(d4est_ops, e_m->u_elem, (P4EST_DIM), f_m, e_m->deg, u_m_on_f_m);
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_linalg_fill_vec
      (
       dudx_m_on_f_m_quad[d],
       0.0,
       face_nodes_m_quad
      );
  }
  
  for (int j = 0; j < (P4EST_DIM); j++){
    for (int i = 0; i < (P4EST_DIM); i++){   
      for (int k = 0; k < face_nodes_m_quad; k++){
        dudx_m_on_f_m_quad[j][k] += drst_dxyz_quad[i][j][k]*dudr_m_on_f_m_quad[i][k];
      }
    }
  }


#ifdef NASTY_TRUMP_DEBUG
  mpi_assert((P4EST_DIM)==2);
  d4est_rst_t rst_lobatto;
  d4est_rst_t rst_compactified;

  rst_lobatto.r = d4est_quadrature_lobatto_get_rst(d4est_ops, d4est_geom, d4est_quad, &face_object, QUAD_OBJECT_MORTAR, QUAD_INTEGRAND_UNKNOWN, e_m->deg, 0);

 
  if (d4est_quadrature_compactified_check_type(d4est_quad)){

 
  rst_compactified.r = d4est_quadrature_compactified_get_rst(d4est_ops,
                                                             d4est_geom,
                                                             d4est_quad,
                                                             &face_object,
                                                             QUAD_OBJECT_MORTAR,
                                                             QUAD_INTEGRAND_UNKNOWN,
                                                             e_m->deg_quad,
                                                             0);

    
  printf("rst_compactified.r[0] = %f\n", rst_compactified.r[0]);
  
  double* xyz_lobatto [(P4EST_DIM)];
  double* xyz_compactified [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(xyz_lobatto, face_nodes_m_lobatto);
  D4EST_ALLOC_DIM_VEC(xyz_compactified, face_nodes_m_quad);
  
  d4est_geometry_compute_xyz_face_analytic
    (
     d4est_ops,
     d4est_geom,
     rst_lobatto,
     e_m->q,
     e_m->dq,
     e_m->tree,
     f_m,
     e_m->deg,
     xyz_lobatto
    );


  d4est_geometry_compute_xyz_face_analytic
    (
     d4est_ops,
     d4est_geom,
     rst_compactified,
     e_m->q,
     e_m->dq,
     e_m->tree,
     f_m,
     e_m->deg,
     xyz_compactified
    );

  for (int i = 0; i < face_nodes_m_lobatto; i++){
    printf("rst_lobatto[%d], x, y, = %f, %f, %f\n",i, rst_lobatto.r[i], xyz_lobatto[0][i], xyz_lobatto[1][i]);
  }

  for (int i = 0; i < face_nodes_m_quad; i++){
    printf("rst_compactified[%d], x, y, = %f, %f, %f\n",i,rst_compactified.r[i], xyz_compactified[0][i], xyz_compactified[1][i]);
  }

  D4EST_FREE_DIM_VEC(xyz_lobatto);
  D4EST_FREE_DIM_VEC(xyz_compactified);
  }  

#endif
  

  if (ip_flux_params->ip_flux_bc_eval == BC_EVAL_ON_LOBATTO_POINTS){

    double* xyz_on_f_m [(P4EST_DIM)];
    D4EST_ALLOC_DIM_VEC(xyz_on_f_m, face_nodes_m_lobatto);
    
    for (int d = 0; d < (P4EST_DIM); d++){

      d4est_operators_apply_slicer(d4est_ops,
                                   e_m->xyz[d],
                                   (P4EST_DIM),
                                   f_m,
                                   e_m->deg,
                                   xyz_on_f_m[d]);

    }

    
    for (int i = 0; i < face_nodes_m_lobatto; i++){
      u_at_bndry_lobatto[i] = u_at_bndry
                              (
                               xyz_on_f_m[0][i],
                               xyz_on_f_m[1][i]
#if (P4EST_DIM)==3
                               ,
                               xyz_on_f_m[2][i]
#endif
                              );
      u_m_on_f_m_min_u_at_bndry_lobatto[i] = u_m_on_f_m[i]
                                             - u_at_bndry_lobatto[i];


#ifdef NASTY_TRUMP_DEBUG
      printf(" u_m_on_f_m = %f, u_at_bndry = %f\n", u_m_on_f_m[i],
             u_at_bndry_lobatto[i]);
#endif
      
    }
    

    d4est_quadrature_interpolate
      (
       d4est_ops,
       d4est_quad,
       d4est_geom,
       &face_object,
       QUAD_OBJECT_MORTAR,
       QUAD_INTEGRAND_UNKNOWN,
       u_m_on_f_m_min_u_at_bndry_lobatto,
       e_m->deg,
       u_m_on_f_m_min_u_at_bndry_quad,
       e_m->deg_quad
      );

    d4est_quadrature_interpolate
      (
       d4est_ops,
       d4est_quad,
       d4est_geom,
       &face_object,
       QUAD_OBJECT_MORTAR,
       QUAD_INTEGRAND_UNKNOWN,
       u_m_on_f_m_min_u_at_bndry_lobatto,
       e_m->deg,
       u_m_on_f_m_min_u_at_bndry_quad,
       e_m->deg_quad
      );  

    DEBUG_PRINT_2ARR_DBL(u_m_on_f_m_min_u_at_bndry_lobatto,
                         u_m_on_f_m_min_u_at_bndry_quad,
                         face_nodes_m_quad
                        );
    
    D4EST_FREE_DIM_VEC(xyz_on_f_m);
  }
  else {
    mpi_abort("At this time we only support BC_EVAL_ON_LOBATTO_POINTS");
  }
  
  for(int i = 0; i < face_nodes_m_quad; i++){
     
    term1_quad[i] = 0.;
    for (int d = 0; d < (P4EST_DIM); d++){
      /* term1_quad[i] += -1.*sj_on_f_m_quad[i] */
      /* *n_on_f_m_quad[d][i] */
      /* *(dudx_m_on_f_m_quad[d][i]); */
      term1_quad[i] += -1.*n_sj_on_f_m_quad[d][i]
                       *(dudx_m_on_f_m_quad[d][i]);

      
    }
    
    for (int l = 0; l < (P4EST_DIM); l++){
      term2_quad[l][i] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        /* term2_quad[l][i] += -.5*sj_on_f_m_quad[i] */
        /* *drst_dxyz_quad[l][d][i] */
        /* *n_on_f_m_quad[d][i] */
        /* *2.*u_m_on_f_m_min_u_at_bndry_quad[i]; */
        term2_quad[l][i] += -.5*drst_dxyz_quad[l][d][i]
                            *n_sj_on_f_m_quad[d][i]
                            *2.*u_m_on_f_m_min_u_at_bndry_quad[i];
      }
    }

    term3_quad[i] = sj_on_f_m_quad[i]
                    *sigma[i]
                    *2.*u_m_on_f_m_min_u_at_bndry_quad[i];
#ifdef NASTY_TRUMP_DEBUG
    /* printf("sj_on_f_m_quad[%d] = %f, sigma = %f, u_m_on_f_m_min_u_term3= %f, term1_quad[i] = %f, term2_quad[0][i] = %f, term3_quad[i] = %f\n", i, sj_on_f_m_quad[i], sigma[i], u_m_on_f_m_min_u_at_bndry_quad[i], term1_quad[i], term2_quad[0][i], term3_quad[i]);     */
    printf("sj, sig, u_m-u, term3 = %.15f, %.15f, %.15f, %.15f\n",
           sj_on_f_m_quad[i],
           sigma[i],
           u_m_on_f_m_min_u_at_bndry_quad[i],
           term3_quad[i]
          );
           

    
#endif
  }
  
  d4est_quadrature_apply_galerkin_integral
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &face_object,
     QUAD_OBJECT_MORTAR,
     QUAD_INTEGRAND_UNKNOWN,
     term1_quad,
     e_m->deg,
     ones_quad,
     e_m->deg_quad,
     VT_w_term1_lobatto
    );

  
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_quadrature_apply_galerkin_integral
      (
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &face_object,
       QUAD_OBJECT_MORTAR,
       QUAD_INTEGRAND_UNKNOWN,
       term2_quad[d],
       e_m->deg,
       ones_quad,
       e_m->deg_quad,
       VT_w_term2_lobatto[d]
      );
  }

  d4est_quadrature_apply_galerkin_integral
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &face_object,
     QUAD_OBJECT_MORTAR,
     QUAD_INTEGRAND_UNKNOWN,
     term3_quad,
     e_m->deg,
     ones_quad,
     e_m->deg_quad,
     VT_w_term3_lobatto
    );
  
  d4est_operators_apply_lift(
                             d4est_ops,
                             VT_w_term1_lobatto,
                             (P4EST_DIM),
                             e_m->deg,
                             f_m,
                             lifted_VT_w_term1_lobatto);


  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_lift(
                               d4est_ops,
                               VT_w_term2_lobatto[d],
                               (P4EST_DIM),
                               e_m->deg,
                               f_m,
                               lifted_VT_w_term2_lobatto[d]);

    d4est_operators_apply_dij_transpose(d4est_ops,
                                        lifted_VT_w_term2_lobatto[d],
                                        (P4EST_DIM),
                                        e_m->deg,
                                        d,
                                        DT_lifted_VT_w_term2_lobatto[d]
                                       );

  }
        

  d4est_operators_apply_lift(
                             d4est_ops,
                             VT_w_term3_lobatto,
                             (P4EST_DIM),
                             e_m->deg,
                             f_m,
                             lifted_VT_w_term3_lobatto);

#ifdef NASTY_DEBUG2
  printf("Element m id = %d, f_m = %d, BOUNDARY, t1 = %f, t20 = %f, t21 = %f, t22 = %f, t3 = %f\n",
         e_m->id,
         f_m,
         util_sum_array_dbl(VT_w_term1_lobatto, face_nodes_m_lobatto),
         util_sum_array_dbl(VT_w_term2_lobatto[0], face_nodes_m_lobatto),
         util_sum_array_dbl(VT_w_term2_lobatto[1], face_nodes_m_lobatto),
#if (P4EST_DIM)==3
         util_sum_array_dbl(VT_w_term2_lobatto[2], face_nodes_m_lobatto),
#else
         -1.,
#endif
         util_sum_array_dbl(VT_w_term3_lobatto, face_nodes_m_lobatto)
        );


  DEBUG_PRINT_ARR_DBL(VT_w_term3_lobatto, face_nodes_m_lobatto);
#endif
  

  
  int volume_nodes_m = d4est_lgl_get_nodes((P4EST_DIM), e_m->deg);
  for (int i = 0; i < volume_nodes_m; i++){
    for (int d = 0; d < (P4EST_DIM); d++){
      e_m->Au_elem[i] += DT_lifted_VT_w_term2_lobatto[d][i];
    }
    e_m->Au_elem[i] += lifted_VT_w_term3_lobatto[i];
    e_m->Au_elem[i] += lifted_VT_w_term1_lobatto[i];
  }
  
  P4EST_FREE(u_at_bndry_lobatto);
  P4EST_FREE(u_at_bndry_quad);  
  P4EST_FREE(sigma);
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_m_on_f_m_quad);
  for (int i = 0; i < (P4EST_DIM); i++) {
    P4EST_FREE(dudr_m_on_f_m[i]);
    P4EST_FREE(dudr_m_on_f_m_quad[i]);
    P4EST_FREE(dudx_m_on_f_m_quad[i]);
    P4EST_FREE(n_on_f_m_quad[i]);
    P4EST_FREE(n_sj_on_f_m_quad[i]);
  }

  P4EST_FREE(sj_on_f_m_quad);
  P4EST_FREE(u_m_on_f_m_min_u_at_bndry_quad);
  P4EST_FREE(u_m_on_f_m_min_u_at_bndry_lobatto);  
  D4EST_FREE_DBYD_MAT(drst_dxyz_quad);
  
  P4EST_FREE(ones_quad);
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

static void
curved_primal_sipg_kronbichler_flux_interface
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
 void* params
)
{
  ip_flux_params_t* ip_flux_params = (ip_flux_params_t*) params;
  double sipg_kronbichler_flux_penalty_prefactor = ip_flux_params->ip_flux_penalty_prefactor;
  penalty_calc_t sipg_kronbichler_flux_penalty_calculate_fcn = ip_flux_params->ip_flux_penalty_calculate_fcn;
  
  int stride;
  int deg_p_lobatto [(P4EST_HALF)];
  int face_nodes_p_lobatto [(P4EST_HALF)];
  /* int deg_p_quad [(P4EST_HALF)]; */
  int face_nodes_p_quad [(P4EST_HALF)];

  int deg_m_lobatto [(P4EST_HALF)];
  int face_nodes_m_lobatto [(P4EST_HALF)];
  /* int deg_m_quad [(P4EST_HALF)]; */
  int face_nodes_m_quad [(P4EST_HALF)];
  
  int nodes_mortar_quad [(P4EST_HALF)];
  int nodes_mortar_lobatto [(P4EST_HALF)];
  int deg_mortar_quad [(P4EST_HALF)];
  int deg_mortar_lobatto [(P4EST_HALF)];
  int deg_p_lobatto_porder [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
  double penalty_mortar [(P4EST_HALF)];

  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);

  /* calculate degs and nodes of each face of (-) side */
  int max_volume_nodes_m_lobatto = 0;
  int total_side_nodes_m_lobatto = 0;
  int total_side_nodes_m_quad = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_lobatto[i] = e_m[i]->deg;

    int volume_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM), e_m[i]->deg);
    max_volume_nodes_m_lobatto = (volume_nodes_m_lobatto > max_volume_nodes_m_lobatto) ? volume_nodes_m_lobatto : max_volume_nodes_m_lobatto;
    
    face_nodes_m_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg);
    face_nodes_m_quad[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg_quad);
    
    total_side_nodes_m_lobatto += face_nodes_m_lobatto[i];
    total_side_nodes_m_quad += face_nodes_m_quad[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_lobatto = 0;
  int total_side_nodes_p_quad = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p_lobatto[i] = e_p_oriented[i]->deg;
    deg_p_lobatto_porder[i] = e_p[i]->deg;
    /* deg_p_quad[i] = e_p_oriented[i]->deg_quad; */

    face_nodes_p_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    face_nodes_p_quad[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg_quad);
    
    total_side_nodes_p_lobatto += face_nodes_p_lobatto[i];
    total_side_nodes_p_quad += face_nodes_p_quad[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar_quad = 0;
  int total_nodes_mortar_lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_quad[i+j] = util_max_int( e_m[i]->deg_quad,
                                           e_p_oriented[j]->deg_quad);
      deg_mortar_lobatto[i+j] = util_max_int( e_m[i]->deg,
                                              e_p_oriented[j]->deg );      
      nodes_mortar_quad[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_quad[i+j] );     
      nodes_mortar_lobatto[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_lobatto[i+j] );     
      total_nodes_mortar_quad += nodes_mortar_quad[i+j];
      total_nodes_mortar_lobatto += nodes_mortar_lobatto[i+j];
      /* penalty_mortar[i+j] = sipg_kronbichler_flux_penalty_calculate_fcn */
                            /* ( */
                             /* e_m[i]->deg, */
                             /* (e_m[i]->volume/e_m[i]->surface_area[f_m]), */
                             /* e_p_oriented[j]->deg, */
                             /* (e_p_oriented[j]->volume/e_p_oriented[j]->surface_area[f_p]), */
                             /* sipg_kronbichler_flux_penalty_prefactor */
                            /* ); */
      
    }

  int deg_mortar_quad_porder [(P4EST_HALF)];
  int nodes_mortar_quad_porder [(P4EST_HALF)];
  
  for(int i = 0; i < faces_mortar; i++){
    int inew = i;
    if (faces_mortar == (P4EST_HALF)){
      inew = d4est_reference_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, i);
    }
    deg_mortar_quad_porder[inew] = deg_mortar_quad[i];
    nodes_mortar_quad_porder[inew] = nodes_mortar_quad[i];
  }

  /* slices of scalar/vector fields of (-) onto f_m and (+) onto f_p */
  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
  double* u_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p_lobatto);
  
  /* projections of f_m/f_p on to mortar space */
  double* u_m_on_f_m_mortar = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* u_m_on_f_m_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* u_p_on_f_p_mortar = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* u_p_on_f_p_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* dudr_p_on_f_p_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_quad_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar_quad [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_quad_porder [(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_quad [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_quad [(P4EST_DIM)];
  double* drst_dxyz_m_on_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_p_on_mortar_quad_porder [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_m_on_mortar_quad, total_nodes_mortar_quad);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_p_on_mortar_quad_porder, total_nodes_mortar_quad);
  
  for (int i = 0; i < (P4EST_DIM); i++){
    dudr_p_on_f_p_porder[i] = P4EST_ALLOC(double, total_side_nodes_p_lobatto);
    dudr_p_on_f_p_mortar_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    dudr_p_on_f_p_mortar_quad_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    dudx_p_on_f_p_mortar_quad_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    dudx_p_on_f_p_mortar_quad[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
  
    dudr_m_on_f_m[i] = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
    dudr_m_on_f_m_mortar[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    dudr_m_on_f_m_mortar_quad[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    dudx_m_on_f_m_mortar_quad[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
  }

  double* ones_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  d4est_linalg_fill_vec(ones_mortar_quad, 1., total_nodes_mortar_quad);
  
  double* lifted_proj_VT_w_term1_mortar_lobatto = P4EST_ALLOC(double, max_volume_nodes_m_lobatto);
  double* proj_VT_w_term1_mortar_lobatto = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
  double* VT_w_term1_mortar_lobatto = P4EST_ALLOC(double, total_nodes_mortar_lobatto);
  double* term1_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad); 

  double* DT_lifted_proj_VT_w_term2_mortar_lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(DT_lifted_proj_VT_w_term2_mortar_lobatto, max_volume_nodes_m_lobatto);
  double* lifted_proj_VT_w_term2_mortar_lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(lifted_proj_VT_w_term2_mortar_lobatto, max_volume_nodes_m_lobatto);
  double* proj_VT_w_term2_mortar_lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(proj_VT_w_term2_mortar_lobatto, total_side_nodes_m_lobatto);
  double* VT_w_term2_mortar_lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(VT_w_term2_mortar_lobatto, total_nodes_mortar_lobatto);
  double* term2_mortar_quad [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(term2_mortar_quad, total_nodes_mortar_quad);

  double* lifted_proj_VT_w_term3_mortar_lobatto = P4EST_ALLOC(double, max_volume_nodes_m_lobatto);
  double* proj_VT_w_term3_mortar_lobatto = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
  double* VT_w_term3_mortar_lobatto = P4EST_ALLOC(double, total_nodes_mortar_lobatto);
  double* term3_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);

  
  double* sj_on_f_m_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* n_on_f_m_mortar_quad [(P4EST_DIM)];
  double* n_sj_on_f_m_mortar_quad [(P4EST_DIM)];

  double* tmp = P4EST_ALLOC(double, total_side_nodes_p_quad);
 
  for (int i = 0; i < (P4EST_DIM); i++) {
    n_on_f_m_mortar_quad[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    n_sj_on_f_m_mortar_quad[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
  }  
  
  stride = 0;
  for (int i = 0; i < faces_m; i++){
    d4est_operators_apply_slicer
      (
       d4est_ops,
       &(e_m[i]->u_elem[0]),
       (P4EST_DIM),
       f_m,
       e_m[i]->deg,
       &u_m_on_f_m[stride]
      );
    
    stride += face_nodes_m_lobatto[i];
  }
 
  stride = 0;
  for (int i = 0; i < faces_p; i++){
    d4est_operators_apply_slicer
      (
       d4est_ops,
       &(e_p_oriented[i]->u_elem[0]),
       (P4EST_DIM),
       f_p,
       e_p_oriented[i]->deg,
       tmp
      );
    
    d4est_operators_reorient_face_data
      (
       d4est_ops,
       tmp,
       ((P4EST_DIM) - 1),
       e_p_oriented[i]->deg,
       orientation,
       f_m,
       f_p,
       &u_p_on_f_p[stride]
      );
    
    stride += face_nodes_p_lobatto[i];
  }


  /* project (-)-side u trace vector onto mortar space */
  d4est_mortars_project_side_onto_mortar_space
    (
     d4est_ops,
     u_m_on_f_m,
     faces_m,
     deg_m_lobatto,
     u_m_on_f_m_mortar,
     faces_mortar,
     deg_mortar_quad
    );

  /* project (+)-side u trace vector onto mortar space */
  d4est_mortars_project_side_onto_mortar_space
    (
     d4est_ops,
     u_p_on_f_p,
     faces_p,
     deg_p_lobatto,
     u_p_on_f_p_mortar,
     faces_mortar,
     deg_mortar_quad
    );
  
  stride = 0;

  p4est_qcoord_t mortar_q0_forder [(P4EST_HALF)][(P4EST_DIM)];
  p4est_qcoord_t mortar_dq_forder;
  p4est_qcoord_t mortar_q0_porder [(P4EST_HALF)][(P4EST_DIM)];
  p4est_qcoord_t mortar_dq_porder;
  d4est_quadrature_mortar_t mortar_face_object_forder [(P4EST_HALF)];
  d4est_quadrature_mortar_t mortar_face_object_porder [(P4EST_HALF)];
  
  d4est_mortars_compute_qcoords_on_mortar
    (
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     faces_m,
     faces_mortar,
     f_m,
     mortar_q0_forder,
     &mortar_dq_forder
    );

  d4est_mortars_compute_qcoords_on_mortar
    (
     e_p[0]->tree,
     e_p[0]->q,
     e_p[0]->dq,
     faces_p,
     faces_mortar,
     f_p,
     mortar_q0_porder,
     &mortar_dq_porder
    );


  for (int f = 0; f < faces_mortar; f++){
    mortar_face_object_forder[f].dq = mortar_dq_forder;
    mortar_face_object_forder[f].tree = e_m[0]->tree;
    mortar_face_object_forder[f].mortar_side_id = mortar_side_id_m;
    mortar_face_object_forder[f].mortar_subface_id = f;   
    mortar_face_object_forder[f].face = f_m;
    mortar_face_object_forder[f].q[0] = mortar_q0_forder[f][0];
    mortar_face_object_forder[f].q[1] = mortar_q0_forder[f][1];
#if (P4EST_DIM)==3
    mortar_face_object_forder[f].q[2] = mortar_q0_forder[f][2];
#endif

    mortar_face_object_porder[f].dq = mortar_dq_porder;
    mortar_face_object_porder[f].tree = e_p[0]->tree;
    mortar_face_object_porder[f].face = f_p;
    mortar_face_object_porder[f].mortar_side_id = mortar_side_id_p;
    mortar_face_object_porder[f].mortar_subface_id = f;
    mortar_face_object_porder[f].q[0] = mortar_q0_porder[f][0];
    mortar_face_object_porder[f].q[1] = mortar_q0_porder[f][1];
#if (P4EST_DIM)==3
    mortar_face_object_porder[f].q[2] = mortar_q0_porder[f][2];
#endif    
  }

  
  
  for (int f = 0; f < faces_mortar; f++){

    d4est_quadrature_interpolate
      (
       d4est_ops,
       d4est_quad,
       d4est_geom,
       &mortar_face_object_forder,
       QUAD_OBJECT_MORTAR,
       QUAD_INTEGRAND_UNKNOWN,
       &u_m_on_f_m_mortar[stride],
       deg_mortar_quad[f],
       &u_m_on_f_m_mortar_quad[stride],
       deg_mortar_quad[f]
      );
    
    d4est_quadrature_interpolate
      (
       d4est_ops,
       d4est_quad,
       d4est_geom,
       &mortar_face_object_forder,
       QUAD_OBJECT_MORTAR,
       QUAD_INTEGRAND_UNKNOWN,
       &u_p_on_f_p_mortar[stride],
       deg_mortar_quad[f],
       &u_p_on_f_p_mortar_quad[stride],
       deg_mortar_quad[f]
      );

    stride += nodes_mortar_quad[f];
  }

  
  /* For each component of the vector */
  for (int d = 0; d < (P4EST_DIM); d++){

    stride = 0;
    for (int i = 0; i < faces_m; i++){    
      d4est_operators_apply_slicer
        (
         d4est_ops,
         &e_m[i]->dudr_elem[d][0],
         (P4EST_DIM),
         f_m,
         e_m[i]->deg,
         &dudr_m_on_f_m[d][stride]
        );
      stride += face_nodes_m_lobatto[i];
    }
    
    stride = 0;
    for (int i = 0; i < faces_p; i++){
      d4est_operators_apply_slicer
        (
         d4est_ops,
         &e_p[i]->dudr_elem[d][0],
         (P4EST_DIM),
         f_p,
         e_p[i]->deg,
         &dudr_p_on_f_p_porder[d][stride]
        );
      stride += d4est_lgl_get_nodes((P4EST_DIM)-1, e_p[i]->deg);
    }

    d4est_mortars_project_side_onto_mortar_space
      (
       d4est_ops,
       dudr_p_on_f_p_porder[d],
       faces_p,
       deg_p_lobatto_porder,
       dudr_p_on_f_p_mortar_porder[d],
       faces_mortar,
       deg_mortar_quad_porder
      );

    d4est_mortars_project_side_onto_mortar_space
      (
       d4est_ops,
       dudr_m_on_f_m[d],
       faces_m,
       deg_m_lobatto,
       dudr_m_on_f_m_mortar[d],
       faces_mortar,
       deg_mortar_quad
      );

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){

      d4est_quadrature_interpolate
        (
         d4est_ops,
         d4est_quad,
         d4est_geom,
         &mortar_face_object_forder,
         QUAD_OBJECT_MORTAR,
         QUAD_INTEGRAND_UNKNOWN,
         &dudr_m_on_f_m_mortar[d][stride],
         deg_mortar_quad[f],
         &dudr_m_on_f_m_mortar_quad[d][stride],
         deg_mortar_quad[f]
        );

      stride += nodes_mortar_quad[f];
    }

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){

      d4est_quadrature_interpolate
        (
         d4est_ops,
         d4est_quad,
         d4est_geom,
         &mortar_face_object_porder,
         QUAD_OBJECT_MORTAR,
         QUAD_INTEGRAND_UNKNOWN,
         &dudr_p_on_f_p_mortar_porder[d][stride],
         deg_mortar_quad_porder[f],
         &dudr_p_on_f_p_mortar_quad_porder[d][stride],
         deg_mortar_quad_porder[f]
        );

      stride += nodes_mortar_quad_porder[f];
    }
  }

  double* j_div_sj_on_f_m_mortar_quad = NULL;
  double* j_div_sj_on_f_p_mortar_quad_porder = NULL;
  double* j_div_sj_on_f_p_mortar_quad_porder_oriented = NULL;

  if(ip_flux_params->ip_flux_h_calc == H_EQ_J_DIV_SJ || ip_flux_params->ip_flux_h_calc == H_EQ_J_DIV_SJ_MIN){
    j_div_sj_on_f_m_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
    j_div_sj_on_f_p_mortar_quad_porder =  P4EST_ALLOC(double, total_nodes_mortar_quad);
    j_div_sj_on_f_p_mortar_quad_porder_oriented =  P4EST_ALLOC(double, total_nodes_mortar_quad);
  }
  
  d4est_mortars_compute_geometric_data_on_mortar
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     QUAD_INTEGRAND_UNKNOWN,
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     mortar_side_id_m,
     faces_m,
     faces_mortar,
     &deg_mortar_quad[0],
     f_m,
     drst_dxyz_m_on_mortar_quad,
     sj_on_f_m_mortar_quad,
     n_on_f_m_mortar_quad,
     n_sj_on_f_m_mortar_quad,
     j_div_sj_on_f_m_mortar_quad,
     COMPUTE_NORMAL_USING_JACOBIAN
    );
  

  d4est_mortars_compute_geometric_data_on_mortar
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     QUAD_INTEGRAND_UNKNOWN,
     e_p[0]->tree,
     e_p[0]->q,
     e_p[0]->dq,
     mortar_side_id_p,
     faces_p,
     faces_mortar,
     &deg_mortar_quad_porder[0],
     f_p,
     drst_dxyz_p_on_mortar_quad_porder,
     NULL,
     NULL,
     NULL,
     j_div_sj_on_f_p_mortar_quad_porder,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_linalg_fill_vec
      (
       dudx_m_on_f_m_mortar_quad[d],
       0.0,
       total_nodes_mortar_quad
      );


    d4est_linalg_fill_vec
      (
       dudx_p_on_f_p_mortar_quad_porder[d],
       0.0,
       total_nodes_mortar_quad
      );
  }
  
  for (int j = 0; j < (P4EST_DIM); j++){
    for (int i = 0; i < (P4EST_DIM); i++){   
      for (int k = 0; k < total_nodes_mortar_quad; k++){
        dudx_m_on_f_m_mortar_quad[j][k] +=
          drst_dxyz_m_on_mortar_quad[i][j][k]*dudr_m_on_f_m_mortar_quad[i][k];
        dudx_p_on_f_p_mortar_quad_porder[j][k] += drst_dxyz_p_on_mortar_quad_porder[i][j][k]*dudr_p_on_f_p_mortar_quad_porder[i][k];
      }
    }    
  }
  
  int face_mortar_stride = 0;
  for (int face = 0; face < faces_mortar; face++){
    int face_p = face;
    if (faces_mortar == (P4EST_HALF))
      face_p = d4est_reference_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, face);

    int oriented_face_mortar_stride = 0;
    for (int b = 0; b < face_p; b++){
      oriented_face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad_porder[b]);
    }


    for (int d = 0; d < (P4EST_DIM); d++){
      d4est_operators_reorient_face_data
        (
         d4est_ops,
         &dudx_p_on_f_p_mortar_quad_porder[d][oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_quad[face],
         orientation,
         f_m,
         f_p,
         &dudx_p_on_f_p_mortar_quad[d][face_mortar_stride]
        );
    }

    if (j_div_sj_on_f_p_mortar_quad_porder != NULL){
      d4est_operators_reorient_face_data
        (
         d4est_ops,
         &j_div_sj_on_f_p_mortar_quad_porder[oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_quad[face],
         orientation,
         f_m,
         f_p,
         &j_div_sj_on_f_p_mortar_quad_porder_oriented[face_mortar_stride]
        );
    }    
    
    
    face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad[face]);
  }

  double hm_min, hp_min;
  if (ip_flux_params->ip_flux_h_calc == H_EQ_J_DIV_SJ_MIN){
    hp_min = util_min_dbl_array(j_div_sj_on_f_p_mortar_quad_porder_oriented, total_nodes_mortar_quad);
    hm_min = util_min_dbl_array(j_div_sj_on_f_m_mortar_quad,total_nodes_mortar_quad);
  }
  
  stride = 0;
  double* sigma = P4EST_ALLOC(double, total_nodes_mortar_quad);
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_quad[f]; k++){
      int ks = k + stride;
      if (ip_flux_params->ip_flux_h_calc == H_EQ_VOLUME_DIV_AREA){
        mpi_abort("H_EQ_VOLUME_DIV_AREA IS NO LONGER SUPPORTED\n");
        /* sigma[ks] = penalty_mortar[f]; */
        /* printf("sigma[ks] = %f\n", sigma[ks]); */
      }
      else if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ){
        double hp = j_div_sj_on_f_p_mortar_quad_porder_oriented[ks];
        double hm = j_div_sj_on_f_m_mortar_quad[ks];
        sigma[ks] = sipg_kronbichler_flux_penalty_calculate_fcn
                    (
                     e_m[f]->deg,
                     hm,
                     e_p_oriented[f]->deg,
                     hp,
                     sipg_kronbichler_flux_penalty_prefactor
                    );

      }
      else if (ip_flux_params->ip_flux_h_calc == H_EQ_J_DIV_SJ_MIN){
        sigma[ks] = sipg_kronbichler_flux_penalty_calculate_fcn
                    (
                     e_m[f]->deg,
                     hm_min,
                     e_p_oriented[f]->deg,
                     hp_min,
                     sipg_kronbichler_flux_penalty_prefactor
                    );
      }
      else {
        mpi_abort("Select j_DIV_SJ or VOLUME_DIV_AREA for ip_flux_h_calc");
      }

    }
    stride += nodes_mortar_quad[f];
  }

  
  stride = 0;
  int stride_lobatto = 0;
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_quad[f]; k++){
      int ks = k + stride;

      term1_mortar_quad[ks] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        term1_mortar_quad[ks] += -1.*n_sj_on_f_m_mortar_quad[d][ks]
                                 *.5*(dudx_p_on_f_p_mortar_quad[d][ks] + dudx_m_on_f_m_mortar_quad[d][ks]);

#ifdef NASTY_DEBUG
        printf("dudx_p_on_f_p_mortar_quad[%d][%d] = %f\n", d, ks, dudx_p_on_f_p_mortar_quad[d][ks]);
        printf("dudx_m_on_f_m_mortar_quad[%d][%d] = %f\n", d, ks, dudx_m_on_f_m_mortar_quad[d][ks]);
        printf("n_sj_on_f_m_mortar_quad[%d][%d] = %f\n", d, ks, n_sj_on_f_m_mortar_quad[d][ks]);
#endif       
      }
        
      for (int l = 0; l < (P4EST_DIM); l++){
        term2_mortar_quad[l][ks] = 0.;
        for (int d = 0; d < (P4EST_DIM); d++){
          term2_mortar_quad[l][ks] += -.5*drst_dxyz_m_on_mortar_quad[l][d][ks]
                                      *n_sj_on_f_m_mortar_quad[d][ks]
                                      *(u_m_on_f_m_mortar_quad[ks] - u_p_on_f_p_mortar_quad[ks]);
          
        }
      }
        
      term3_mortar_quad[ks] = sj_on_f_m_mortar_quad[ks]*sigma[ks]*(u_m_on_f_m_mortar_quad[ks] - u_p_on_f_p_mortar_quad[ks]);
    }

   
      
    d4est_quadrature_apply_galerkin_integral
      (
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &mortar_face_object_forder,
       QUAD_OBJECT_MORTAR,
       QUAD_INTEGRAND_UNKNOWN,
       &term1_mortar_quad[stride],
       deg_mortar_lobatto[f],
       &ones_mortar_quad[stride],
       deg_mortar_quad[f],
       &VT_w_term1_mortar_lobatto[stride_lobatto]
      );

  
    for (int d = 0; d < (P4EST_DIM); d++){
      
      d4est_quadrature_apply_galerkin_integral
        (
         d4est_ops,
         d4est_geom,
         d4est_quad,
         &mortar_face_object_forder,
         QUAD_OBJECT_MORTAR,
         QUAD_INTEGRAND_UNKNOWN,
         &term2_mortar_quad[d][stride],
         deg_mortar_lobatto[f],
         &ones_mortar_quad[stride],
         deg_mortar_quad[f],
         &VT_w_term2_mortar_lobatto[d][stride_lobatto]
        );
      
    }

    d4est_quadrature_apply_galerkin_integral
      (
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &mortar_face_object_forder,
       QUAD_OBJECT_MORTAR,
       QUAD_INTEGRAND_UNKNOWN,
       &term3_mortar_quad[stride],
       deg_mortar_lobatto[f],
       &ones_mortar_quad[stride],
       deg_mortar_quad[f],
       &VT_w_term3_mortar_lobatto[stride_lobatto]
      );
      
    stride += nodes_mortar_quad[f];
    stride_lobatto += nodes_mortar_lobatto[f];
  }
  

  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_mortars_project_mass_mortar_onto_side
      (
       d4est_ops,
       VT_w_term2_mortar_lobatto[d],
       faces_mortar,
       deg_mortar_lobatto,
       proj_VT_w_term2_mortar_lobatto[d],
       faces_m,
       deg_m_lobatto
      );
  }

  d4est_mortars_project_mass_mortar_onto_side
    (
     d4est_ops,
     VT_w_term1_mortar_lobatto,
     faces_mortar,
     deg_mortar_lobatto,
     proj_VT_w_term1_mortar_lobatto,
     faces_m,
     deg_m_lobatto
    );

  d4est_mortars_project_mass_mortar_onto_side
    (
     d4est_ops,
     VT_w_term3_mortar_lobatto,
     faces_mortar,
     deg_mortar_lobatto,
     proj_VT_w_term3_mortar_lobatto,
     faces_m,
     deg_m_lobatto
    );

  /* copy result back to element */
  stride = 0;
  for (int f = 0; f < faces_m; f++){
    if (e_m_is_ghost[f] == 0){

      d4est_operators_apply_lift(
                                 d4est_ops,
                                 &proj_VT_w_term1_mortar_lobatto[stride],
                                 (P4EST_DIM),
                                 e_m[f]->deg,
                                 f_m,
                                 lifted_proj_VT_w_term1_mortar_lobatto);


      for (int d = 0; d < (P4EST_DIM); d++){
        d4est_operators_apply_lift(
                                   d4est_ops,
                                   &proj_VT_w_term2_mortar_lobatto[d][stride],
                                   (P4EST_DIM),
                                   e_m[f]->deg,
                                   f_m,
                                   lifted_proj_VT_w_term2_mortar_lobatto[d]);

        d4est_operators_apply_dij_transpose(d4est_ops,
                                            lifted_proj_VT_w_term2_mortar_lobatto[d],
                                            (P4EST_DIM),
                                            e_m[f]->deg,
                                            d,
                                            DT_lifted_proj_VT_w_term2_mortar_lobatto[d]
                                           );

      }
        

      d4est_operators_apply_lift(
                                 d4est_ops,
                                 &proj_VT_w_term3_mortar_lobatto[stride],
                                 (P4EST_DIM),
                                 e_m[f]->deg,
                                 f_m,
                                 lifted_proj_VT_w_term3_mortar_lobatto);



        
      int volume_nodes_m = d4est_lgl_get_nodes((P4EST_DIM), e_m[f]->deg);
      for (int i = 0; i < volume_nodes_m; i++){
        for (int d = 0; d < (P4EST_DIM); d++){
          e_m[f]->Au_elem[i] += DT_lifted_proj_VT_w_term2_mortar_lobatto[d][i];
        }
        e_m[f]->Au_elem[i] += lifted_proj_VT_w_term3_mortar_lobatto[i];
        e_m[f]->Au_elem[i] += lifted_proj_VT_w_term1_mortar_lobatto[i];
      }
    }
    stride += face_nodes_m_lobatto[f];
  }   

#ifdef NASTY_DEBUG
  printf("Interface Flux Integral\n");
  printf("Element m id = %d, f_m = %d, Element_p id = %d, f_p = %d\n", e_m[0]->id, f_m, e_p[0]->id,f_p);
  DEBUG_PRINT_ARR_DBL_SUM(proj_VT_w_term1_mortar_lobatto, total_nodes_mortar_lobatto);
  for (int d = 0; d < (P4EST_DIM); d++){
    DEBUG_PRINT_ARR_DBL_SUM(proj_VT_w_term2_mortar_lobatto[d], total_nodes_mortar_lobatto);
  }
  DEBUG_PRINT_ARR_DBL_SUM(proj_VT_w_term3_mortar_lobatto, total_nodes_mortar_lobatto);
#endif

#ifdef NASTY_DEBUG2
  printf("Element m id = %d, f_m = %d, Element_p id = %d, f_p = %d, t1 = %f, t20 = %f, t21 = %f, t22 = %f, t3 = %f\n",
         e_m[0]->id,
         f_m,
         e_p[0]->id,
         f_p,
         util_sum_array_dbl(proj_VT_w_term1_mortar_lobatto, total_nodes_mortar_lobatto),
         util_sum_array_dbl(proj_VT_w_term2_mortar_lobatto[0], total_nodes_mortar_lobatto),
         util_sum_array_dbl(proj_VT_w_term2_mortar_lobatto[1], total_nodes_mortar_lobatto),
#if (P4EST_DIM)==3
         util_sum_array_dbl(proj_VT_w_term2_mortar_lobatto[2], total_nodes_mortar_lobatto),
#else
         -1.,
#endif
         util_sum_array_dbl(proj_VT_w_term3_mortar_lobatto, total_nodes_mortar_lobatto)
        );
#endif
  
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_p_on_f_p);

  D4EST_FREE_DBYD_MAT(drst_dxyz_m_on_mortar_quad);
  D4EST_FREE_DBYD_MAT(drst_dxyz_p_on_mortar_quad_porder);

  P4EST_FREE(sigma);

  P4EST_FREE(u_m_on_f_m_mortar);
  P4EST_FREE(u_m_on_f_m_mortar_quad);
  P4EST_FREE(u_p_on_f_p_mortar);
  P4EST_FREE(u_p_on_f_p_mortar_quad);
  for (int i = 0; i < (P4EST_DIM); i++){   
    P4EST_FREE(dudr_p_on_f_p_porder[i]);
    P4EST_FREE(dudr_p_on_f_p_mortar_porder[i]);
    P4EST_FREE(dudr_p_on_f_p_mortar_quad_porder[i]);
    P4EST_FREE(dudx_p_on_f_p_mortar_quad_porder[i]);
    P4EST_FREE(dudx_p_on_f_p_mortar_quad[i]);
    P4EST_FREE(dudr_m_on_f_m[i]);
    P4EST_FREE(dudr_m_on_f_m_mortar[i]);
    P4EST_FREE(dudr_m_on_f_m_mortar_quad[i]);
    P4EST_FREE(dudx_m_on_f_m_mortar_quad[i]);
  }
  P4EST_FREE(sj_on_f_m_mortar_quad);
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(n_on_f_m_mortar_quad[i]);
    P4EST_FREE(n_sj_on_f_m_mortar_quad[i]);
  }

  if (j_div_sj_on_f_m_mortar_quad != NULL){
    P4EST_FREE(j_div_sj_on_f_m_mortar_quad);
    P4EST_FREE(j_div_sj_on_f_p_mortar_quad_porder);
    P4EST_FREE(j_div_sj_on_f_p_mortar_quad_porder_oriented);
  }
  P4EST_FREE(ones_mortar_quad);
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
  
  P4EST_FREE(tmp);  
}


static double
d4est_poisson_flux_sipg_penalty_meanp_sqr_over_meanh
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double mean_p = .5*(deg_m + deg_p);
  double mean_p_sqr = mean_p*mean_p;
  double mean_h = .5*(h_m + h_p);
  return (penalty_prefactor*mean_p_sqr)/mean_h;
}

static double
d4est_poisson_flux_sipg_penalty_maxp_sqr_over_minh
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double max_deg = (deg_m > deg_p) ? deg_m : deg_p;
  double min_h = (h_m < h_p) ? h_m : h_p;
  return (penalty_prefactor*(max_deg)*(max_deg))/min_h;
}

static void
ip_flux_params_get_string_from_h_calc
(
 h_calc_method_t h_calc,
 char string [50]
)
{
  if (h_calc == H_EQ_J_DIV_SJ){
    string = "H_EQ_J_DIV_SJ";
  }
  else if (h_calc == H_EQ_J_DIV_SJ_MIN){
    string = "H_EQ_J_DIV_SJ_MIN";
  }
  else if (h_calc == H_EQ_VOLUME_DIV_AREA){
    string = "H_EQ_VOLUME_DIV_AREA";
  }
  else {
    string = "NOT_SET";
  }
}

static void
ip_flux_params_get_string_from_bc_eval
(
 bc_eval_t bc_eval,
 char string [50]
)
{
  if (bc_eval == BC_EVAL_ON_QUADRATURE_POINTS){
    string = "BC_EVAL_ON_QUADRATURE_POINTS";
  }
  else if (bc_eval == BC_EVAL_ON_LOBATTO_POINTS){
    string = "BC_EVAL_ON_LOBATTO_POINTS";
  }
  else {
    string = "NOT_SET";
  }
}

static void
ip_flux_params_get_string_from_penalty_fcn
(
 penalty_calc_t fcn,
 char string [50]
)
{
  if (fcn == ip_flux_params_penalty_maxp_sqr_over_minh){
    string = "maxp_sqr_over_minh";
  }
  else if (fcn == ip_flux_params_penalty_meanp_sqr_over_meanh){
    string = "meanp_sqr_over_meanh";
  }
  else {
    string = "not_set";
  }
}

static penalty_calc_t
ip_flux_params_get_penalty_fcn_from_string
(
 const char* string
)
{
  if (util_match(string,"maxp_sqr_over_minh")){
    return ip_flux_params_penalty_maxp_sqr_over_minh;
  }
  else if (util_match(string,"meanp_sqr_over_meanh")){
    return ip_flux_params_penalty_meanp_sqr_over_meanh;
  }
  else {
    mpi_abort("This ip flux penalty calculation fcn does not exist");
    return NULL;
  }
}

static
int ip_flux_params_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  ip_flux_params_t* pconfig = (ip_flux_params_t*)user;
  if (util_match_couple(section,"ip_flux_params",name,"name")) {
    mpi_assert(pconfig->name[0] == '*');
    snprintf (pconfig->name, sizeof(pconfig->name), "%s", value);
  }
  
  if (util_match_couple(section,"ip_flux_params",name,"ip_flux_penalty_prefactor")) {
    mpi_assert(pconfig->ip_flux_penalty_prefactor == -1);
    pconfig->ip_flux_penalty_prefactor = atof(value);
  }
  else if (util_match_couple(section,"ip_flux_params",name,"ip_flux_penalty_calculate_fcn")) {
    mpi_assert(pconfig->ip_flux_penalty_calculate_fcn == NULL);
    pconfig->ip_flux_penalty_calculate_fcn = ip_flux_params_get_penalty_fcn_from_string(value);
  }
  else if (util_match_couple(section,"ip_flux_params",name,"ip_flux_bc_eval")) {
    mpi_assert(pconfig->ip_flux_bc_eval == BC_EVAL_NOTSET);
    if(util_match(value, "BC_EVAL_ON_QUADRATURE_POINTS")){
      pconfig->ip_flux_bc_eval = BC_EVAL_ON_QUADRATURE_POINTS;
    }
    else if (util_match(value, "BC_EVAL_ON_LOBATTO_POINTS")){
      pconfig->ip_flux_bc_eval = BC_EVAL_ON_LOBATTO_POINTS;
    }
    else {
      printf("ip_flux_params_bc_eval = %s\n", value);
      mpi_abort("ip_flux_params_bc_eval is not set to BC_EVAL_ON_LOBATTO_POINTS or BC_EVAL_ON_QUADRATURE_POINTS\n");
    }
  }
  else if (util_match_couple(section,"ip_flux_params",name,"ip_flux_h_calc")) {
    mpi_assert(pconfig->ip_flux_h_calc == H_EQ_NOTSET);
    if(util_match(value, "H_EQ_J_DIV_SJ")){
      pconfig->ip_flux_h_calc = H_EQ_J_DIV_SJ;
    }
    else if(util_match(value, "H_EQ_J_DIV_SJ_MIN")){
      pconfig->ip_flux_h_calc = H_EQ_J_DIV_SJ_MIN;
    }
    else if (util_match(value, "H_EQ_VOLUME_DIV_AREA")){
      pconfig->ip_flux_h_calc = H_EQ_VOLUME_DIV_AREA;
    }
    else {
      printf("ip_flux_params_h_calc = %s\n", value);
      mpi_abort("ip_flux_params_h_calc is not set to H_EQ_J_DIV_SJ or H_EQ_VOLUME_DIV_AREA\n");
    }
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

static void
ip_flux_params_input
(
 p4est_t* p4est,
 const char* printf_prefix,
 const char* input_file,
 ip_flux_params_t* input
)
{
  /* set defaults */
  input->name[0] = '*';
  input->ip_flux_bc_eval = BC_EVAL_NOTSET;
  input->ip_flux_h_calc = H_EQ_NOTSET;
  input->ip_flux_penalty_calculate_fcn = NULL;
  input->ip_flux_penalty_prefactor = -1.;

  if (ini_parse(input_file, ip_flux_params_input_handler, input) < 0) {
    mpi_abort("Can't load input file");
  }

  D4EST_CHECK_INPUT("ip_flux_params", input->ip_flux_penalty_prefactor, -1);
  D4EST_CHECK_INPUT("ip_flux_params", input->ip_flux_penalty_calculate_fcn, NULL);
  D4EST_CHECK_INPUT("ip_flux_params", input->ip_flux_h_calc, H_EQ_NOTSET);
  D4EST_CHECK_INPUT("ip_flux_params", input->ip_flux_bc_eval, BC_EVAL_NOTSET);
  D4EST_CHECK_INPUT("ip_flux_params", input->name[0], '*');
  
  char penalty_calculate_fcn [50];
  char h_eq [50];
  char bc_eval [50];  

  ip_flux_params_get_string_from_penalty_fcn (input->ip_flux_penalty_calculate_fcn,penalty_calculate_fcn);
  ip_flux_params_get_string_from_bc_eval (input->ip_flux_bc_eval,bc_eval);
  ip_flux_params_get_string_from_h_calc (input->ip_flux_h_calc,h_eq);
  
  if(p4est->mpirank == 0){
    printf("%s: ip_flux_params_penalty_prefactor = %f\n", printf_prefix, input->ip_flux_penalty_prefactor);
    printf("%s: ip_flux_params_penalty_calculate_fcn = %s\n", printf_prefix, penalty_calculate_fcn);
    printf("%s: ip_flux_params_h_calc = %s\n", printf_prefix, h_eq);
    printf("%s: ip_flux_params_bc_eval = %s\n", printf_prefix, bc_eval);
  }
}

ip_flux_params_t*
ip_flux_params_new
(
 p4est_t* p4est,
 const char* print_prefix,
 const char* input_file
)
{
  ip_flux_params_t* ip_flux_params = P4EST_ALLOC(ip_flux_params_t, 1); 
  ip_flux_params_input(p4est, print_prefix, input_file, ip_flux_params);
  return ip_flux_params;
}

void
ip_flux_params_destroy
(
 ip_flux_params_t* params
){
  P4EST_FREE(params);
}


d4est_mortar_fcn_ptrs_t
d4est_poisson_flux_sipg_fetch_fcns
(
 d4est_grid_fcn_t bndry_fcn,
 d4est_poisson_flux_sipg_params_t* sipg_params
)
{  
  d4est_mortar_fcn_ptrs_t sipg_flux_fcns;
  sipg_flux_fcns.flux_interface_fcn = d4est_poisson_flux_sipg_interface;
  sipg_flux_fcns.flux_boundary_fcn = d4est_poisson_flux_sipg_dirichlet;
  sipg_flux_fcns.bndry_fcn = bndry_fcn;
  sipg_flux_fcns.params = (void*)sipg_params;
  
  return sipg_flux_fcns;
}

