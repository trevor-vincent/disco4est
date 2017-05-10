#include "../Utilities/util.h"
#include "../dGMath/dgmath.h"
#include "../ElementData/curved_element_data.h"
#include "../LinearAlgebra/linalg.h"
#include "../Flux/curved_Gauss_primal_sipg_kronbichler_flux_fcns.h"



static void
curved_primal_sipg_kronbichler_flux_dirichlet
(
 curved_element_data_t* e_m,
 int f_m,
 grid_fcn_t bndry_fcn,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* geom,
 void* params
)
{
  ip_flux_params_t* ip_flux_params = (ip_flux_params_t*) params;
  double sipg_kronbichler_flux_penalty_prefactor = ip_flux_params->ip_flux_penalty_prefactor;
  penalty_calc_t sipg_kronbichler_flux_penalty_calculate_fcn = ip_flux_params->ip_flux_penalty_calculate_fcn;

  
  grid_fcn_t u_at_bndry = bndry_fcn;
  int face_nodes_m_lobatto = dgmath_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_quad = dgmath_get_nodes((P4EST_DIM) - 1, e_m->deg_integ);

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


  int volume_nodes_m_lobatto = dgmath_get_nodes((P4EST_DIM), e_m->deg);
  
  double* ones_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  linalg_fill_vec(ones_quad, 1., face_nodes_m_quad);
  
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

  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m->tree,
     e_m->q,
     e_m->dq,
     1,
     1,
     &e_m->deg_integ,
     f_m,
     drst_dxyz_quad,
     sj_on_f_m_quad,
     n_on_f_m_quad,
     n_sj_on_f_m_quad,
     J_div_SJ_quad,
     geom->geom_quad_type,
     geom,
     dgmath_jit_dbase,
     COMPUTE_NORMAL_USING_JACOBIAN
    );
  
  double* sigma = P4EST_ALLOC(double, face_nodes_m_quad);
  double h, h_min;
  if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ_MIN){
    h_min = util_min_dbl_array(J_div_SJ_quad, face_nodes_m_quad);
  }
    
  for (int i = 0; i < face_nodes_m_quad; i++){
    if (ip_flux_params->ip_flux_h_calc == H_EQ_VOLUME_DIV_AREA){
      h = (e_m->volume/e_m->surface_area[f_m]);
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
  }

  if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ || ip_flux_params->ip_flux_h_calc == H_EQ_J_DIV_SJ_MIN){
    P4EST_FREE(J_div_SJ_quad);
  }

  
  
  
  for (int d = 0; d < (P4EST_DIM); d++){
    
    dgmath_apply_slicer
      (
       dgmath_jit_dbase,
       e_m->dudr_elem[d],
       (P4EST_DIM),
       f_m,
       e_m->deg,
       dudr_m_on_f_m[d]
      );

    dgmath_interp(dgmath_jit_dbase,
                  dudr_m_on_f_m[d],
                  QUAD_LOBATTO,
                  e_m->deg,
                  dudr_m_on_f_m_quad[d],
                  geom->geom_quad_type,
                  e_m->deg_integ,
                  (P4EST_DIM)-1);
   
  }


  
  dgmath_apply_slicer(dgmath_jit_dbase, e_m->u_storage, (P4EST_DIM), f_m, e_m->deg, u_m_on_f_m);
  for (int d = 0; d < (P4EST_DIM); d++){
    linalg_fill_vec
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


  if (ip_flux_params->ip_flux_bc_eval == BC_EVAL_ON_LOBATTO_POINTS){

    double* xyz_on_f_m [(P4EST_DIM)];
    D4EST_ALLOC_DIM_VEC(xyz_on_f_m, face_nodes_m_lobatto);
    
    for (int d = 0; d < (P4EST_DIM); d++){

      dgmath_apply_slicer(dgmath_jit_dbase,
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
    }
    

    dgmath_interp(dgmath_jit_dbase,
                  u_m_on_f_m_min_u_at_bndry_lobatto,
                  QUAD_LOBATTO,
                  e_m->deg,
                  u_m_on_f_m_min_u_at_bndry_quad,
                  geom->geom_quad_type,
                  e_m->deg_integ,
                  (P4EST_DIM)-1);
    

    dgmath_interp(dgmath_jit_dbase,
                  u_at_bndry_lobatto,
                  QUAD_LOBATTO,
                  e_m->deg,
                  u_at_bndry_quad,
                  geom->geom_quad_type,
                  e_m->deg_integ,
                  (P4EST_DIM)-1);
    
    

    D4EST_FREE_DIM_VEC(xyz_on_f_m);
  }
  else if (ip_flux_params->ip_flux_bc_eval == BC_EVAL_ON_QUADRATURE_POINTS){

    double* xyz_on_f_m_quad [(P4EST_DIM)];
    D4EST_ALLOC_DIM_VEC(xyz_on_f_m_quad, face_nodes_m_quad);
    

    dgmath_interp(dgmath_jit_dbase,
                  u_m_on_f_m,
                  QUAD_LOBATTO,
                  e_m->deg,
                  u_m_on_f_m_quad,
                  geom->geom_quad_type,
                  e_m->deg_integ,
                  (P4EST_DIM)-1);
    
    
    d4est_geometry_compute_xyz_face_analytic
      (
       dgmath_jit_dbase,
       e_m->q,
       e_m->dq,
       e_m->tree,
       f_m,
       geom,
       geom->geom_quad_type,
       e_m->deg_integ,
       xyz_on_f_m_quad
      );

    
    for (int i = 0; i < face_nodes_m_quad; i++){
      u_at_bndry_quad[i]
        = u_at_bndry
        (
         xyz_on_f_m_quad[0][i],
         xyz_on_f_m_quad[1][i]
#if (P4EST_DIM)==3
                             ,
         xyz_on_f_m_quad[2][i]
#endif
        );

      //      printf("boundary = %f\n", u_at_bndry_quad[i]);
    }
    
    for (int i = 0; i < face_nodes_m_quad; i++){
      u_m_on_f_m_min_u_at_bndry_quad[i] = u_m_on_f_m_quad[i]
                                           - u_at_bndry_quad[i];
    }      


    D4EST_FREE_DIM_VEC(xyz_on_f_m_quad);
  }
  else {
    mpi_abort("Select either BC_EVAL_ON_GAUSS_POINTS or BC_EVAL_ON_LOBATTO_POINTS");
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
  }

  dgmath_apply_curved_galerkin_integral
    (
     dgmath_jit_dbase,
     term1_quad,
     e_m->deg,
     ones_quad,
     e_m->deg_integ,
     geom->geom_quad_type,
     (P4EST_DIM)-1,
     VT_w_term1_lobatto
    );

  
  for (int d = 0; d < (P4EST_DIM); d++){
    dgmath_apply_curved_galerkin_integral
      (
       dgmath_jit_dbase,
       term2_quad[d],
       e_m->deg,
       ones_quad,
       e_m->deg_integ,
       geom->geom_quad_type,
       (P4EST_DIM)-1,
       VT_w_term2_lobatto[d]
      );
  }

  
    dgmath_apply_curved_galerkin_integral
    (
     dgmath_jit_dbase,
     term3_quad,
     e_m->deg,
     ones_quad,
     e_m->deg_integ,
     geom->geom_quad_type,
     (P4EST_DIM)-1,
     VT_w_term3_lobatto
    );  
  
  dgmath_apply_LIFT(
                    dgmath_jit_dbase,
                    VT_w_term1_lobatto,
                    (P4EST_DIM),
                    e_m->deg,
                    f_m,
                    lifted_VT_w_term1_lobatto);


  for (int d = 0; d < (P4EST_DIM); d++){
    dgmath_apply_LIFT(
                      dgmath_jit_dbase,
                      VT_w_term2_lobatto[d],
                      (P4EST_DIM),
                      e_m->deg,
                      f_m,
                      lifted_VT_w_term2_lobatto[d]);

    dgmath_apply_Dij_transpose(dgmath_jit_dbase,
                               lifted_VT_w_term2_lobatto[d],
                               (P4EST_DIM),
                               e_m->deg,
                               d,
                               DT_lifted_VT_w_term2_lobatto[d]
                              );

  }
        

  dgmath_apply_LIFT(
                    dgmath_jit_dbase,
                    VT_w_term3_lobatto,
                    (P4EST_DIM),
                    e_m->deg,
                    f_m,
                    lifted_VT_w_term3_lobatto);


      
      
  int volume_nodes_m = dgmath_get_nodes((P4EST_DIM), e_m->deg);
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
  D4EST_FREE_DIM_VEC(term2_quad)
  D4EST_FREE_DIM_VEC(VT_w_term2_lobatto)
  D4EST_FREE_DIM_VEC(lifted_VT_w_term2_lobatto)
  D4EST_FREE_DIM_VEC(DT_lifted_VT_w_term2_lobatto)
  P4EST_FREE(term3_quad);
  P4EST_FREE(VT_w_term3_lobatto);
  P4EST_FREE(lifted_VT_w_term3_lobatto);
}

static void
curved_primal_sipg_kronbichler_flux_interface
(
 curved_element_data_t** e_m,
 int faces_m,
 int f_m,
 curved_element_data_t** e_p,
 int faces_p,
 int f_p,
 int* e_m_is_ghost,
 int orientation,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* geom,
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

  curved_element_data_t* e_p_oriented [(P4EST_HALF)];
  curved_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);

  /* calculate degs and nodes of each face of (-) side */
  int max_volume_nodes_m_lobatto = 0;
  int total_side_nodes_m_lobatto = 0;
  int total_side_nodes_m_quad = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_lobatto[i] = e_m[i]->deg;

    int volume_nodes_m_lobatto = dgmath_get_nodes((P4EST_DIM), e_m[i]->deg);
    max_volume_nodes_m_lobatto = (volume_nodes_m_lobatto > max_volume_nodes_m_lobatto) ? volume_nodes_m_lobatto : max_volume_nodes_m_lobatto;
    
    face_nodes_m_lobatto[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg);
    face_nodes_m_quad[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg_integ);
    
    total_side_nodes_m_lobatto += face_nodes_m_lobatto[i];
    total_side_nodes_m_quad += face_nodes_m_quad[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_lobatto = 0;
  int total_side_nodes_p_quad = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p_lobatto[i] = e_p_oriented[i]->deg;
    deg_p_lobatto_porder[i] = e_p[i]->deg;
    /* deg_p_quad[i] = e_p_oriented[i]->deg_integ; */

    face_nodes_p_lobatto[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    face_nodes_p_quad[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg_integ);
    
    total_side_nodes_p_lobatto += face_nodes_p_lobatto[i];
    total_side_nodes_p_quad += face_nodes_p_quad[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar_quad = 0;
  int total_nodes_mortar_lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_quad[i+j] = util_max_int( e_m[i]->deg_integ,
                                            e_p_oriented[j]->deg_integ);
      deg_mortar_lobatto[i+j] = util_max_int( e_m[i]->deg,
                                              e_p_oriented[j]->deg );      
      nodes_mortar_quad[i+j] = dgmath_get_nodes( (P4EST_DIM) - 1, deg_mortar_quad[i+j] );     
      nodes_mortar_lobatto[i+j] = dgmath_get_nodes( (P4EST_DIM) - 1, deg_mortar_lobatto[i+j] );     
      total_nodes_mortar_quad += nodes_mortar_quad[i+j];
      total_nodes_mortar_lobatto += nodes_mortar_lobatto[i+j];
      penalty_mortar[i+j] = sipg_kronbichler_flux_penalty_calculate_fcn
                            (
                             e_m[i]->deg,
                             (e_m[i]->volume/e_m[i]->surface_area[f_m]),
                             e_p_oriented[j]->deg,
                             (e_p_oriented[j]->volume/e_p_oriented[j]->surface_area[f_p]),
                             sipg_kronbichler_flux_penalty_prefactor
                            );
      
    }

  int deg_mortar_quad_porder [(P4EST_HALF)];
  int nodes_mortar_quad_porder [(P4EST_HALF)];
  
  for(int i = 0; i < faces_mortar; i++){
    int inew = i;
    if (faces_mortar == (P4EST_HALF)){
      inew = dgmath_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, i);
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
  linalg_fill_vec(ones_mortar_quad, 1., total_nodes_mortar_quad);
  
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
    dgmath_apply_slicer
      (
       dgmath_jit_dbase,
       &(e_m[i]->u_storage[0]),
       (P4EST_DIM),
       f_m,
       e_m[i]->deg,
       &u_m_on_f_m[stride]
      );
    
    stride += face_nodes_m_lobatto[i];
  }
 
  stride = 0;
  for (int i = 0; i < faces_p; i++){
    dgmath_apply_slicer
      (
       dgmath_jit_dbase,
       &(e_p_oriented[i]->u_storage[0]),
       (P4EST_DIM),
       f_p,
       e_p_oriented[i]->deg,
       tmp
      );
    
    dgmath_reorient_face_data
      (
       dgmath_jit_dbase,
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
  dgmath_project_side_onto_mortar_space
    (
     dgmath_jit_dbase,
     u_m_on_f_m,
     faces_m,
     deg_m_lobatto,
     u_m_on_f_m_mortar,
     faces_mortar,
     deg_mortar_quad
    );

  /* project (+)-side u trace vector onto mortar space */
  dgmath_project_side_onto_mortar_space
    (
     dgmath_jit_dbase,
     u_p_on_f_p,
     faces_p,
     deg_p_lobatto,
     u_p_on_f_p_mortar,
     faces_mortar,
     deg_mortar_quad
    );
  
  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    dgmath_interp(dgmath_jit_dbase,
                  &u_m_on_f_m_mortar[stride],
                  QUAD_LOBATTO,
                  deg_mortar_quad[f],
                  &u_m_on_f_m_mortar_quad[stride],
                  geom->geom_quad_type,
                  deg_mortar_quad[f],
                  (P4EST_DIM)-1);

    dgmath_interp(dgmath_jit_dbase,
                  &u_p_on_f_p_mortar[stride],
                  QUAD_LOBATTO,
                  deg_mortar_quad[f],
                  &u_p_on_f_p_mortar_quad[stride],
                  geom->geom_quad_type,
                  deg_mortar_quad[f],
                  (P4EST_DIM)-1);
      
    stride += nodes_mortar_quad[f];
  }

  
  /* For each component of the vector */
  for (int d = 0; d < (P4EST_DIM); d++){

    stride = 0;
    for (int i = 0; i < faces_m; i++){    

      dgmath_apply_slicer
        (
         dgmath_jit_dbase,
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
      dgmath_apply_slicer
        (
         dgmath_jit_dbase,
         &e_p[i]->dudr_elem[d][0],
         (P4EST_DIM),
         f_p,
         e_p[i]->deg,
         &dudr_p_on_f_p_porder[d][stride]
        );
      stride += dgmath_get_nodes((P4EST_DIM)-1, e_p[i]->deg);
    }


    dgmath_project_side_onto_mortar_space
      (
       dgmath_jit_dbase,
       dudr_p_on_f_p_porder[d],
       faces_p,
       deg_p_lobatto_porder,
       dudr_p_on_f_p_mortar_porder[d],
       faces_mortar,
       deg_mortar_quad_porder
      );
  

    dgmath_project_side_onto_mortar_space
      (
       dgmath_jit_dbase,
       dudr_m_on_f_m[d],
       faces_m,
       deg_m_lobatto,
       dudr_m_on_f_m_mortar[d],
       faces_mortar,
       deg_mortar_quad
      );


    stride = 0;
    for (int f = 0; f < faces_mortar; f++){

    dgmath_interp(dgmath_jit_dbase,
                  &dudr_m_on_f_m_mortar[d][stride],
                  QUAD_LOBATTO,
                  deg_mortar_quad[f],
                  &dudr_m_on_f_m_mortar_quad[d][stride],
                  geom->geom_quad_type,
                  deg_mortar_quad[f],
                  (P4EST_DIM)-1);

      
      stride += nodes_mortar_quad[f];
    }

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){

      dgmath_interp(dgmath_jit_dbase,
                    &dudr_p_on_f_p_mortar_porder[d][stride],
                    QUAD_LOBATTO,
                    deg_mortar_quad_porder[f],
                    &dudr_p_on_f_p_mortar_quad_porder[d][stride],
                    geom->geom_quad_type,
                    deg_mortar_quad_porder[f],
                    (P4EST_DIM)-1);
        
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
  
  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     faces_m,
     faces_mortar,
     &deg_mortar_quad[0],
     f_m,
     drst_dxyz_m_on_mortar_quad,
     sj_on_f_m_mortar_quad,
     n_on_f_m_mortar_quad,
     n_sj_on_f_m_mortar_quad,
     j_div_sj_on_f_m_mortar_quad,
     geom->geom_quad_type,
     geom,
     dgmath_jit_dbase,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_p[0]->tree,
     e_p[0]->q,
     e_p[0]->dq,
     faces_p,
     faces_mortar,
     &deg_mortar_quad_porder[0],
     f_p,
     drst_dxyz_p_on_mortar_quad_porder,
     NULL,
     NULL,
     NULL,
     j_div_sj_on_f_p_mortar_quad_porder,
     geom->geom_quad_type,
     geom,
     dgmath_jit_dbase,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  for (int d = 0; d < (P4EST_DIM); d++){
    linalg_fill_vec
      (
       dudx_m_on_f_m_mortar_quad[d],
       0.0,
       total_nodes_mortar_quad
      );


    linalg_fill_vec
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
      face_p = dgmath_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, face);

    int oriented_face_mortar_stride = 0;
    for (int b = 0; b < face_p; b++){
      oriented_face_mortar_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_quad_porder[b]);
    }


    for (int d = 0; d < (P4EST_DIM); d++){
      dgmath_reorient_face_data
        (
         dgmath_jit_dbase,
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
    dgmath_reorient_face_data
      (
       dgmath_jit_dbase,
       &j_div_sj_on_f_p_mortar_quad_porder[oriented_face_mortar_stride],
       (P4EST_DIM)-1,
       deg_mortar_quad[face],
       orientation,
       f_m,
       f_p,
       &j_div_sj_on_f_p_mortar_quad_porder_oriented[face_mortar_stride]
      );
    }    
    
    
    face_mortar_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_quad[face]);
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
        sigma[ks] = penalty_mortar[f];
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
        /* term1_mortar_quad[ks] += -1.*sj_on_f_m_mortar_quad[ks] */
                                  /* *n_on_f_m_mortar_quad[d][ks] */
                                  /* *.5*(dudx_p_on_f_p_mortar_quad[d][ks] + dudx_m_on_f_m_mortar_quad[d][ks]); */
        term1_mortar_quad[ks] += -1.*n_sj_on_f_m_mortar_quad[d][ks]
                                  *.5*(dudx_p_on_f_p_mortar_quad[d][ks] + dudx_m_on_f_m_mortar_quad[d][ks]);
        
      }
        
      for (int l = 0; l < (P4EST_DIM); l++){
        term2_mortar_quad[l][ks] = 0.;
        for (int d = 0; d < (P4EST_DIM); d++){
          /* term2_mortar_quad[l][ks] += -.5*sj_on_f_m_mortar_quad[ks] */
                                       /* *drst_dxyz_m_on_mortar_quad[l][d][ks] */
                                       /* *n_on_f_m_mortar_quad[d][ks] */
                                       /* *(u_m_on_f_m_mortar_quad[ks] - u_p_on_f_p_mortar_quad[ks]); */
          term2_mortar_quad[l][ks] += -.5*drst_dxyz_m_on_mortar_quad[l][d][ks]
                                       *n_sj_on_f_m_mortar_quad[d][ks]
                                       *(u_m_on_f_m_mortar_quad[ks] - u_p_on_f_p_mortar_quad[ks]);
          
        }
      }
        
      term3_mortar_quad[ks] = sj_on_f_m_mortar_quad[ks]*sigma[ks]*(u_m_on_f_m_mortar_quad[ks] - u_p_on_f_p_mortar_quad[ks]);
    }



    
   dgmath_apply_curved_galerkin_integral
      (
       dgmath_jit_dbase,
       &term1_mortar_quad[stride],
       deg_mortar_lobatto[f],
       &ones_mortar_quad[stride],
       deg_mortar_quad[f],
       geom->geom_quad_type,
       (P4EST_DIM)-1,
       &VT_w_term1_mortar_lobatto[stride_lobatto]
      );
    
    for (int d = 0; d < (P4EST_DIM); d++){
      dgmath_apply_curved_galerkin_integral
        (
         dgmath_jit_dbase,
         &term2_mortar_quad[d][stride],
         deg_mortar_lobatto[f],
         &ones_mortar_quad[stride],
         deg_mortar_quad[f],
         geom->geom_quad_type,
         (P4EST_DIM)-1,
         &VT_w_term2_mortar_lobatto[d][stride_lobatto]
        );
    }
      
      dgmath_apply_curved_galerkin_integral
      (
       dgmath_jit_dbase,
       &term3_mortar_quad[stride],
       deg_mortar_lobatto[f],
       &ones_mortar_quad[stride],
       deg_mortar_quad[f],
       geom->geom_quad_type,
       (P4EST_DIM)-1,
       &VT_w_term3_mortar_lobatto[stride_lobatto]
      );

    stride += nodes_mortar_quad[f];
    stride_lobatto += nodes_mortar_lobatto[f];
  }
  

  for (int d = 0; d < (P4EST_DIM); d++){
    dgmath_project_mass_mortar_onto_side
      (
       dgmath_jit_dbase,
       VT_w_term2_mortar_lobatto[d],
       faces_mortar,
       deg_mortar_lobatto,
       proj_VT_w_term2_mortar_lobatto[d],
       faces_m,
       deg_m_lobatto
      );
  }

  dgmath_project_mass_mortar_onto_side
    (
     dgmath_jit_dbase,
     VT_w_term1_mortar_lobatto,
     faces_mortar,
     deg_mortar_lobatto,
     proj_VT_w_term1_mortar_lobatto,
     faces_m,
     deg_m_lobatto
    );

  dgmath_project_mass_mortar_onto_side
    (
     dgmath_jit_dbase,
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

      dgmath_apply_LIFT(
                        dgmath_jit_dbase,
                        &proj_VT_w_term1_mortar_lobatto[stride],
                        (P4EST_DIM),
                        e_m[f]->deg,
                        f_m,
                        lifted_proj_VT_w_term1_mortar_lobatto);


      for (int d = 0; d < (P4EST_DIM); d++){
        dgmath_apply_LIFT(
                          dgmath_jit_dbase,
                          &proj_VT_w_term2_mortar_lobatto[d][stride],
                          (P4EST_DIM),
                          e_m[f]->deg,
                          f_m,
                          lifted_proj_VT_w_term2_mortar_lobatto[d]);

        dgmath_apply_Dij_transpose(dgmath_jit_dbase,
                                   lifted_proj_VT_w_term2_mortar_lobatto[d],
                                   (P4EST_DIM),
                                   e_m[f]->deg,
                                   d,
                                   DT_lifted_proj_VT_w_term2_mortar_lobatto[d]
                                  );

      }
        

      dgmath_apply_LIFT(
                        dgmath_jit_dbase,
                        &proj_VT_w_term3_mortar_lobatto[stride],
                        (P4EST_DIM),
                        e_m[f]->deg,
                        f_m,
                        lifted_proj_VT_w_term3_mortar_lobatto);



        
      int volume_nodes_m = dgmath_get_nodes((P4EST_DIM), e_m[f]->deg);
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

curved_flux_fcn_ptrs_t
curved_primal_sipg_kronbichler_flux_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 ip_flux_params_t* curved_quad_sipg_params
)
{  
  curved_flux_fcn_ptrs_t curved_primal_sipg_kronbichler_flux_fcns;
  curved_primal_sipg_kronbichler_flux_fcns.flux_interface_fcn = curved_primal_sipg_kronbichler_flux_interface;
  curved_primal_sipg_kronbichler_flux_fcns.flux_boundary_fcn = curved_primal_sipg_kronbichler_flux_dirichlet;
  curved_primal_sipg_kronbichler_flux_fcns.bndry_fcn = bndry_fcn;
  curved_primal_sipg_kronbichler_flux_fcns.params = (void*)curved_quad_sipg_params;
  
  return curved_primal_sipg_kronbichler_flux_fcns;
}
