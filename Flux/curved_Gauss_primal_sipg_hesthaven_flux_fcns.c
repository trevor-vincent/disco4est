#include "../Utilities/util.h"
#include "../dGMath/d4est_operators.h"
#include "../ElementData/d4est_element_data.h"
#include "../LinearAlgebra/d4est_linalg.h"
#include "../Flux/curved_gauss_primal_sipg_hesthaven_flux_fcns.h"

static void
curved_gauss_primal_sipg_hesthaven_flux_dirichlet_withbcterms
(
 d4est_element_data_t* e_m,
 int f_m,
 grid_fcn_t bndry_fcn,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* geom,
 void* params
)
{
  ip_flux_params_t* ip_flux_params = (ip_flux_params_t*) params;
  double sipg_flux_penalty_prefactor = ip_flux_params->ip_flux_penalty_prefactor;
  penalty_calc_t sipg_flux_penalty_calculate_fcn = ip_flux_params->ip_flux_penalty_calculate_fcn;
 
  grid_fcn_t u_at_bndry = bndry_fcn;
  int face_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_gauss = d4est_lgl_get_nodes((P4EST_DIM) - 1, e_m->deg_quad);

  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* u_at_bndry_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  double* u_at_bndry_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* u_m_on_f_m_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);

  double* dudr_m_on_f_m [(P4EST_DIM)];
  double* dudr_m_on_f_m_gauss [(P4EST_DIM)];
  double* dudx_m_on_f_m_gauss [(P4EST_DIM)];
  
  double* xyz_on_f_m [(P4EST_DIM)];
  double* u_m_on_f_m_min_u_at_bndry_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  double* u_m_on_f_m_min_u_at_bndry_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* sj_on_f_m_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  double* n_on_f_m_gauss [(P4EST_DIM)];
  double* n_sj_on_f_m_gauss [(P4EST_DIM)];

  double* drst_dxyz_gauss [(P4EST_DIM)][(P4EST_DIM)]; D4EST_ALLOC_DBYD_MAT(drst_dxyz_gauss, face_nodes_m_gauss);
  
  for (int d = 0; d < (P4EST_DIM); d++) {
    xyz_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m_lobatto);
    dudr_m_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m_lobatto);
    dudr_m_on_f_m_gauss[d] = P4EST_ALLOC(double, face_nodes_m_gauss);
    dudx_m_on_f_m_gauss[d] = P4EST_ALLOC(double, face_nodes_m_gauss);
    n_on_f_m_gauss[d] = P4EST_ALLOC(double, face_nodes_m_gauss);
    n_sj_on_f_m_gauss[d] = P4EST_ALLOC(double, face_nodes_m_gauss);
  }

  int volume_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM), e_m->deg);
  
  double* ones_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  d4est_linalg_fill_vec(ones_gauss, 1., face_nodes_m_gauss);
  
  double* term2_gauss [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(term2_gauss, face_nodes_m_gauss);
  double* VT_w_term2_lobatto [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(VT_w_term2_lobatto, face_nodes_m_lobatto);
  double* lifted_VT_w_term2_lobatto [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(lifted_VT_w_term2_lobatto, volume_nodes_m_lobatto);
  double* DT_lifted_VT_w_term2_lobatto [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(DT_lifted_VT_w_term2_lobatto, volume_nodes_m_lobatto);
 
  double* term3_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  double* VT_w_term3_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* lifted_VT_w_term3_lobatto = P4EST_ALLOC(double, volume_nodes_m_lobatto);
  
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
    
    d4est_operators_interp_lobatto_to_GL
      (
       d4est_ops,
       dudr_m_on_f_m[d],
       e_m->deg,
       e_m->deg_quad,
       dudr_m_on_f_m_gauss[d],
       (P4EST_DIM)-1
      );
   
  }

  double* J_div_SJ_gauss = NULL;
  if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ){
    J_div_SJ_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  }

  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m->tree,
     e_m->q,
     e_m->dq,
     1,
     1,
     &e_m->deg_quad,
     f_m,
     drst_dxyz_gauss,
     sj_on_f_m_gauss,
     n_on_f_m_gauss,
     n_sj_on_f_m_gauss,
     J_div_SJ_gauss,
     GAUSS,
     geom,
     d4est_ops,
     COMPUTE_NORMAL_USING_JACOBIAN
    );
  
  
  double* sigma = P4EST_ALLOC(double, face_nodes_m_gauss);
  double h;
  for (int i = 0; i < face_nodes_m_gauss; i++){
    if (ip_flux_params->ip_flux_h_calc == H_EQ_VOLUME_DIV_AREA){
      h = (e_m->volume/e_m->surface_area[f_m]);
    }
    else if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ){
      h = J_div_SJ_gauss[i];
    }

    sigma[i] = sipg_flux_penalty_calculate_fcn
               (
                e_m->deg,
                h,
                e_m->deg,
                h,
                sipg_flux_penalty_prefactor
               );
  }

  if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ){
    P4EST_FREE(J_div_SJ_gauss);
  }

  



  
  d4est_operators_apply_slicer(d4est_ops, e_m->u_elem, (P4EST_DIM), f_m, e_m->deg, u_m_on_f_m);
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_linalg_fill_vec
      (
       dudx_m_on_f_m_gauss[d],
       0.0,
       face_nodes_m_gauss
      );
  }
  
  for (int j = 0; j < (P4EST_DIM); j++){
    for (int i = 0; i < (P4EST_DIM); i++){   
      for (int k = 0; k < face_nodes_m_gauss; k++){
        dudx_m_on_f_m_gauss[j][k] += drst_dxyz_gauss[i][j][k]*dudr_m_on_f_m_gauss[i][k];
      }
    }
  }

  if (ip_flux_params->ip_flux_bc_eval == BC_EVAL_ON_LOBATTO_POINTS){

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
    }
    
    d4est_operators_interp_lobatto_to_GL
      (
       d4est_ops,
       u_m_on_f_m_min_u_at_bndry_lobatto,
       e_m->deg,
       e_m->deg_quad,
       u_m_on_f_m_min_u_at_bndry_gauss,
       (P4EST_DIM)-1
      );

    d4est_operators_interp_lobatto_to_GL
      (
       d4est_ops,
       u_at_bndry_lobatto,
       e_m->deg,
       e_m->deg_quad,
       u_at_bndry_gauss,
       (P4EST_DIM)-1
      );
  
  }
  else if (ip_flux_params->ip_flux_bc_eval == BC_EVAL_ON_GAUSS_POINTS){

    double* xyz_on_f_m_gauss [(P4EST_DIM)];
    D4EST_ALLOC_DIM_VEC(xyz_on_f_m_gauss, face_nodes_m_gauss);
    
    d4est_operators_interp_lobatto_to_GL
      (
       d4est_ops,
       u_m_on_f_m,
       e_m->deg,
       e_m->deg_quad,
       u_m_on_f_m_gauss,
       (P4EST_DIM)-1
      );
    
    d4est_geometry_compute_xyz_face_analytic
      (
       d4est_ops,
       e_m->q,
       e_m->dq,
       e_m->tree,
       f_m,
       geom,
       GAUSS,
       e_m->deg_quad,
       xyz_on_f_m_gauss
      );

    
    for (int i = 0; i < face_nodes_m_gauss; i++){
      u_at_bndry_gauss[i]
        = u_at_bndry
        (
         xyz_on_f_m_gauss[0][i],
         xyz_on_f_m_gauss[1][i]
#if (P4EST_DIM)==3
                             ,
         xyz_on_f_m_gauss[2][i]
#endif
        );
    }
    
    for (int i = 0; i < face_nodes_m_gauss; i++){
      u_m_on_f_m_min_u_at_bndry_gauss[i] = u_m_on_f_m_gauss[i]
                                           - u_at_bndry_gauss[i];
    }      


    D4EST_FREE_DIM_VEC(xyz_on_f_m_gauss);
  }
  else {
    mpi_abort("Select either BC_EVAL_ON_GAUSS_POINTS or BC_EVAL_ON_LOBATTO_POINTS");
  }
  
  for(int i = 0; i < face_nodes_m_gauss; i++){
    for (int l = 0; l < (P4EST_DIM); l++){
      term2_gauss[l][i] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        /* term2_gauss[l][i] += sj_on_f_m_gauss[i] */
                             /* *drst_dxyz_gauss[l][d][i] */
                             /* *n_on_f_m_gauss[d][i] */
                             /* *u_at_bndry_gauss[i]; */
        term2_gauss[l][i] += drst_dxyz_gauss[l][d][i]
                             *n_sj_on_f_m_gauss[d][i]
                             *u_at_bndry_gauss[i];
        
      }
    }
    term3_gauss[i] = sj_on_f_m_gauss[i]
                     *sigma[i]
                     *u_m_on_f_m_min_u_at_bndry_gauss[i];
  }

  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_curvedgaussMass_ongaussNodeVec
      (
       d4est_ops,
       term2_gauss[d],
       e_m->deg,
       ones_gauss,
       e_m->deg_quad,
       (P4EST_DIM)-1,
       VT_w_term2_lobatto[d]
      );
  }

  
  d4est_operators_apply_curvedgaussMass_ongaussNodeVec
    (
     d4est_ops,
     term3_gauss,
     e_m->deg,
     ones_gauss,
     e_m->deg_quad,
     (P4EST_DIM)-1,
     VT_w_term3_lobatto
    );  


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

  int volume_nodes_m = d4est_lgl_get_nodes((P4EST_DIM), e_m->deg);
  for (int i = 0; i < volume_nodes_m; i++){
    for (int d = 0; d < (P4EST_DIM); d++){
      e_m->Au_elem[i] += DT_lifted_VT_w_term2_lobatto[d][i];
    }
    e_m->Au_elem[i] += lifted_VT_w_term3_lobatto[i];
  }

  P4EST_FREE(sigma);
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_m_on_f_m_gauss);
  P4EST_FREE(u_at_bndry_gauss);
  P4EST_FREE(u_at_bndry_lobatto);
  for (int i = 0; i < (P4EST_DIM); i++) {
    P4EST_FREE(dudr_m_on_f_m[i]);
    P4EST_FREE(dudr_m_on_f_m_gauss[i]);
    P4EST_FREE(dudx_m_on_f_m_gauss[i]);
    P4EST_FREE(xyz_on_f_m[i]);
  }

  P4EST_FREE(sj_on_f_m_gauss);
  P4EST_FREE(u_m_on_f_m_min_u_at_bndry_gauss);
  P4EST_FREE(u_m_on_f_m_min_u_at_bndry_lobatto);
  
  for (int d = 0; d < (P4EST_DIM); d++){
    P4EST_FREE(n_on_f_m_gauss[d]);
    P4EST_FREE(n_sj_on_f_m_gauss[d]);
  }
  D4EST_FREE_DBYD_MAT(drst_dxyz_gauss);

  P4EST_FREE(ones_gauss);
  D4EST_FREE_DIM_VEC(term2_gauss)
  D4EST_FREE_DIM_VEC(VT_w_term2_lobatto)
  D4EST_FREE_DIM_VEC(lifted_VT_w_term2_lobatto)
  D4EST_FREE_DIM_VEC(DT_lifted_VT_w_term2_lobatto)
  P4EST_FREE(term3_gauss);
  P4EST_FREE(VT_w_term3_lobatto);
  P4EST_FREE(lifted_VT_w_term3_lobatto);
}


static void
curved_gauss_primal_sipg_hesthaven_flux_interface
(
 d4est_element_data_t** e_m,
 int faces_m,
 int f_m,
 d4est_element_data_t** e_p,
 int faces_p,
 int f_p,
 int* e_m_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* geom,
 void* params
)
{
  ip_flux_params_t* ip_flux_params = (ip_flux_params_t*) params;
  double sipg_flux_penalty_prefactor = ip_flux_params->ip_flux_penalty_prefactor;
  penalty_calc_t sipg_flux_penalty_calculate_fcn = ip_flux_params->ip_flux_penalty_calculate_fcn;
  
  int stride;
  int deg_p_lobatto [(P4EST_HALF)];
  int face_nodes_p_lobatto [(P4EST_HALF)];
  /* int deg_p_gauss [(P4EST_HALF)]; */
  int face_nodes_p_gauss [(P4EST_HALF)];

  int deg_m_lobatto [(P4EST_HALF)];
  int face_nodes_m_lobatto [(P4EST_HALF)];
  /* int deg_m_gauss [(P4EST_HALF)]; */
  int face_nodes_m_gauss [(P4EST_HALF)];
  
  int nodes_mortar_gauss [(P4EST_HALF)];
  int nodes_mortar_lobatto [(P4EST_HALF)];
  int deg_mortar_gauss [(P4EST_HALF)];
  int deg_mortar_lobatto [(P4EST_HALF)];
  int deg_p_lobatto_porder [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
  double penalty_mortar [(P4EST_HALF)];

  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);

  /* calculate degs and nodes of each face of (-) side */
  int max_volume_nodes_m_lobatto = 0;
  int total_side_nodes_m_lobatto = 0;
  int total_side_nodes_m_gauss = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_lobatto[i] = e_m[i]->deg;

    int volume_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM), e_m[i]->deg);
    max_volume_nodes_m_lobatto = (volume_nodes_m_lobatto > max_volume_nodes_m_lobatto) ? volume_nodes_m_lobatto : max_volume_nodes_m_lobatto;
    
    face_nodes_m_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg);
    face_nodes_m_gauss[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg_quad);
    
    total_side_nodes_m_lobatto += face_nodes_m_lobatto[i];
    total_side_nodes_m_gauss += face_nodes_m_gauss[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_lobatto = 0;
  int total_side_nodes_p_gauss = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p_lobatto[i] = e_p_oriented[i]->deg;
    deg_p_lobatto_porder[i] = e_p[i]->deg;
    /* deg_p_gauss[i] = e_p_oriented[i]->deg_quad; */

    face_nodes_p_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    face_nodes_p_gauss[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg_quad);
    
    total_side_nodes_p_lobatto += face_nodes_p_lobatto[i];
    total_side_nodes_p_gauss += face_nodes_p_gauss[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar_gauss = 0;
  int total_nodes_mortar_lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_gauss[i+j] = util_max_int( e_m[i]->deg_quad,
                                            e_p_oriented[j]->deg_quad);
      deg_mortar_lobatto[i+j] = util_max_int( e_m[i]->deg,
                                              e_p_oriented[j]->deg );      
      nodes_mortar_gauss[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_gauss[i+j] );     
      nodes_mortar_lobatto[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_lobatto[i+j] );     
      total_nodes_mortar_gauss += nodes_mortar_gauss[i+j];
      total_nodes_mortar_lobatto += nodes_mortar_lobatto[i+j];
      penalty_mortar[i+j] = sipg_flux_penalty_calculate_fcn
                            (
                             e_m[i]->deg,
                             (e_m[i]->volume/e_m[i]->surface_area[f_m]),
                             e_p_oriented[j]->deg,
                             (e_p_oriented[j]->volume/e_p_oriented[j]->surface_area[f_p]),
                             sipg_flux_penalty_prefactor
                            );
      
    }

  int deg_mortar_gauss_porder [(P4EST_HALF)];
  int nodes_mortar_gauss_porder [(P4EST_HALF)];
  
  for(int i = 0; i < faces_mortar; i++){
    int inew = i;
    if (faces_mortar == (P4EST_HALF)){
      inew = d4est_reference_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, i);
    }
    deg_mortar_gauss_porder[inew] = deg_mortar_gauss[i];
    nodes_mortar_gauss_porder[inew] = nodes_mortar_gauss[i];
  }

  /* slices of scalar/vector fields of (-) onto f_m and (+) onto f_p */
  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
  double* u_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p_lobatto);
  
  /* projections of f_m/f_p on to mortar space */
  double* u_m_on_f_m_mortar = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  double* u_m_on_f_m_mortar_gauss = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  double* u_p_on_f_p_mortar = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  double* u_p_on_f_p_mortar_gauss = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  double* dudr_p_on_f_p_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_gauss_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar_gauss [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_gauss_porder [(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_gauss [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_gauss [(P4EST_DIM)];
  double* drst_dxyz_m_on_mortar_gauss [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_p_on_mortar_gauss_porder [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_m_on_mortar_gauss, total_nodes_mortar_gauss);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_p_on_mortar_gauss_porder, total_nodes_mortar_gauss);
  
  for (int i = 0; i < (P4EST_DIM); i++){
    dudr_p_on_f_p_porder[i] = P4EST_ALLOC(double, total_side_nodes_p_lobatto);
    dudr_p_on_f_p_mortar_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
    dudr_p_on_f_p_mortar_gauss_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
    dudx_p_on_f_p_mortar_gauss_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
    dudx_p_on_f_p_mortar_gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  
    dudr_m_on_f_m[i] = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
    dudr_m_on_f_m_mortar[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
    dudr_m_on_f_m_mortar_gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
    dudx_m_on_f_m_mortar_gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  }

  double* ones_mortar_gauss = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  d4est_linalg_fill_vec(ones_mortar_gauss, 1., total_nodes_mortar_gauss);
  
  double* lifted_proj_VT_w_term1_mortar_lobatto = P4EST_ALLOC(double, max_volume_nodes_m_lobatto);
  double* proj_VT_w_term1_mortar_lobatto = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
  double* VT_w_term1_mortar_lobatto = P4EST_ALLOC(double, total_nodes_mortar_lobatto);
  double* term1_mortar_gauss = P4EST_ALLOC(double, total_nodes_mortar_gauss); 

  double* DT_lifted_proj_VT_w_term2_mortar_lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(DT_lifted_proj_VT_w_term2_mortar_lobatto, max_volume_nodes_m_lobatto);
  double* lifted_proj_VT_w_term2_mortar_lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(lifted_proj_VT_w_term2_mortar_lobatto, max_volume_nodes_m_lobatto);
  double* proj_VT_w_term2_mortar_lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(proj_VT_w_term2_mortar_lobatto, total_side_nodes_m_lobatto);
  double* VT_w_term2_mortar_lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(VT_w_term2_mortar_lobatto, total_nodes_mortar_lobatto);
  double* term2_mortar_gauss [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(term2_mortar_gauss, total_nodes_mortar_gauss);

  double* lifted_proj_VT_w_term3_mortar_lobatto = P4EST_ALLOC(double, max_volume_nodes_m_lobatto);
  double* proj_VT_w_term3_mortar_lobatto = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
  double* VT_w_term3_mortar_lobatto = P4EST_ALLOC(double, total_nodes_mortar_lobatto);
  double* term3_mortar_gauss = P4EST_ALLOC(double, total_nodes_mortar_gauss);

  
  double* sj_on_f_m_mortar_gauss = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  double* n_on_f_m_mortar_gauss [(P4EST_DIM)];
  double* n_sj_on_f_m_mortar_gauss [(P4EST_DIM)];

  double* tmp = P4EST_ALLOC(double, total_side_nodes_p_gauss);
 
  for (int i = 0; i < (P4EST_DIM); i++) {
    n_on_f_m_mortar_gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
    n_sj_on_f_m_mortar_gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
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
     deg_mortar_gauss
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
     deg_mortar_gauss
    );
  
  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    d4est_operators_interp_lobatto_to_GL(d4est_ops, &u_m_on_f_m_mortar[stride], deg_mortar_gauss[f], deg_mortar_gauss[f], &u_m_on_f_m_mortar_gauss[stride], (P4EST_DIM)-1);
    d4est_operators_interp_lobatto_to_GL(d4est_ops, &u_p_on_f_p_mortar[stride], deg_mortar_gauss[f], deg_mortar_gauss[f], &u_p_on_f_p_mortar_gauss[stride], (P4EST_DIM)-1);
    stride += nodes_mortar_gauss[f];
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
      for (int d = 0; d < (P4EST_DIM); d++){
        d4est_operators_apply_slicer
          (
           d4est_ops,
           &e_p[i]->dudr_elem[d][0],
           (P4EST_DIM),
           f_p,
           e_p[i]->deg,
           &dudr_p_on_f_p_porder[d][stride]
          );
      }
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
       deg_mortar_gauss_porder
      );
  

    d4est_mortars_project_side_onto_mortar_space
      (
       d4est_ops,
       dudr_m_on_f_m[d],
       faces_m,
       deg_m_lobatto,
       dudr_m_on_f_m_mortar[d],
       faces_mortar,
       deg_mortar_gauss
      );


    stride = 0;
    for (int f = 0; f < faces_mortar; f++){
      d4est_operators_interp_lobatto_to_GL(d4est_ops, &dudr_m_on_f_m_mortar[d][stride], deg_mortar_gauss[f], deg_mortar_gauss[f], &dudr_m_on_f_m_mortar_gauss[d][stride], (P4EST_DIM)-1);
      stride += nodes_mortar_gauss[f];
    }

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){
      for (int d = 0; d < (P4EST_DIM); d++){
        d4est_operators_interp_lobatto_to_GL(d4est_ops, &dudr_p_on_f_p_mortar_porder[d][stride], deg_mortar_gauss_porder[f], deg_mortar_gauss_porder[f], &dudr_p_on_f_p_mortar_gauss_porder[d][stride], (P4EST_DIM)-1);
      }
      stride += nodes_mortar_gauss_porder[f];
    }
  }

  double* j_div_sj_on_f_m_mortar_gauss = NULL;
  double* j_div_sj_on_f_p_mortar_gauss_porder = NULL;
  double* j_div_sj_on_f_p_mortar_gauss_porder_oriented = NULL;

  if(ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ){
    j_div_sj_on_f_m_mortar_gauss = P4EST_ALLOC(double, total_nodes_mortar_gauss);
    j_div_sj_on_f_p_mortar_gauss_porder =  P4EST_ALLOC(double, total_nodes_mortar_gauss);
    j_div_sj_on_f_p_mortar_gauss_porder_oriented =  P4EST_ALLOC(double, total_nodes_mortar_gauss);
  }
  
  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     faces_m,
     faces_mortar,
     &deg_mortar_gauss[0],
     f_m,
     drst_dxyz_m_on_mortar_gauss,
     sj_on_f_m_mortar_gauss,
     n_on_f_m_mortar_gauss,
     n_sj_on_f_m_mortar_gauss,
     j_div_sj_on_f_m_mortar_gauss,
     GAUSS,
     geom,
     d4est_ops,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_p[0]->tree,
     e_p[0]->q,
     e_p[0]->dq,
     faces_p,
     faces_mortar,
     &deg_mortar_gauss_porder[0],
     f_p,
     drst_dxyz_p_on_mortar_gauss_porder,
     NULL,
     NULL,
     NULL,
     j_div_sj_on_f_p_mortar_gauss_porder,
     GAUSS,
     geom,
     d4est_ops,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_linalg_fill_vec
      (
       dudx_m_on_f_m_mortar_gauss[d],
       0.0,
       total_nodes_mortar_gauss
      );


    d4est_linalg_fill_vec
      (
       dudx_p_on_f_p_mortar_gauss_porder[d],
       0.0,
       total_nodes_mortar_gauss
      );
  }
  
  for (int j = 0; j < (P4EST_DIM); j++){
    for (int i = 0; i < (P4EST_DIM); i++){   
      for (int k = 0; k < total_nodes_mortar_gauss; k++){
        dudx_m_on_f_m_mortar_gauss[j][k] +=
          drst_dxyz_m_on_mortar_gauss[i][j][k]*dudr_m_on_f_m_mortar_gauss[i][k];
        dudx_p_on_f_p_mortar_gauss_porder[j][k] += drst_dxyz_p_on_mortar_gauss_porder[i][j][k]*dudr_p_on_f_p_mortar_gauss_porder[i][k];
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
      oriented_face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_gauss_porder[b]);
    }


    for (int d = 0; d < (P4EST_DIM); d++){
      d4est_operators_reorient_face_data
        (
         d4est_ops,
         &dudx_p_on_f_p_mortar_gauss_porder[d][oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_gauss[face],
         orientation,
         f_m,
         f_p,
         &dudx_p_on_f_p_mortar_gauss[d][face_mortar_stride]
        );
    }


    if (j_div_sj_on_f_p_mortar_gauss_porder != NULL){
    d4est_operators_reorient_face_data
      (
       d4est_ops,
       &j_div_sj_on_f_p_mortar_gauss_porder[oriented_face_mortar_stride],
       (P4EST_DIM)-1,
       deg_mortar_gauss[face],
       orientation,
       f_m,
       f_p,
       &j_div_sj_on_f_p_mortar_gauss_porder_oriented[face_mortar_stride]
      );
    }    
    
    
    face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_gauss[face]);
  }


  stride = 0;
  double* sigma = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_gauss[f]; k++){
      int ks = k + stride;
      if (ip_flux_params->ip_flux_h_calc == H_EQ_VOLUME_DIV_AREA){
        sigma[ks] = penalty_mortar[f];
      }
      else if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ){
        double hp = j_div_sj_on_f_p_mortar_gauss_porder_oriented[ks];
        double hm = j_div_sj_on_f_m_mortar_gauss[ks];

        sigma[ks] = sipg_flux_penalty_calculate_fcn
                (
                 e_m[f]->deg,
                 hm,
                 e_p_oriented[f]->deg,
                 hp,
                 sipg_flux_penalty_prefactor
                );
      }
      else {
        mpi_abort("Select j_DIV_SJ or VOLUME_DIV_AREA for ip_flux_h_calc");
      }

     
    }
    stride += nodes_mortar_gauss[f];
  }


  
  stride = 0;
  int stride_lobatto = 0;
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_gauss[f]; k++){
      int ks = k + stride;

    term1_mortar_gauss[ks] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        /* term1_mortar_gauss[ks] += -1.*sj_on_f_m_mortar_gauss[ks] */
                                  /* *n_on_f_m_mortar_gauss[d][ks] */
                                  /* *.5*(dudx_p_on_f_p_mortar_gauss[d][ks] + dudx_m_on_f_m_mortar_gauss[d][ks]); */
        term1_mortar_gauss[ks] += -1.*n_sj_on_f_m_mortar_gauss[d][ks]
                                  *.5*(dudx_p_on_f_p_mortar_gauss[d][ks] + dudx_m_on_f_m_mortar_gauss[d][ks]);
        
      }
        
      for (int l = 0; l < (P4EST_DIM); l++){
        term2_mortar_gauss[l][ks] = 0.;
        for (int d = 0; d < (P4EST_DIM); d++){
          /* term2_mortar_gauss[l][ks] += -.5*sj_on_f_m_mortar_gauss[ks] */
                                       /* *drst_dxyz_m_on_mortar_gauss[l][d][ks] */
                                       /* *n_on_f_m_mortar_gauss[d][ks] */
                                       /* *(u_m_on_f_m_mortar_gauss[ks] - u_p_on_f_p_mortar_gauss[ks]); */
          term2_mortar_gauss[l][ks] += -.5*drst_dxyz_m_on_mortar_gauss[l][d][ks]
                                       *n_sj_on_f_m_mortar_gauss[d][ks]
                                       *(u_m_on_f_m_mortar_gauss[ks] - u_p_on_f_p_mortar_gauss[ks]);
          
        }
      }
        
      term3_mortar_gauss[ks] = sj_on_f_m_mortar_gauss[ks]*sigma[ks]*(u_m_on_f_m_mortar_gauss[ks] - u_p_on_f_p_mortar_gauss[ks]);
    }

    
   d4est_operators_apply_curvedgaussMass_ongaussNodeVec
      (
       d4est_ops,
       &term1_mortar_gauss[stride],
       deg_mortar_lobatto[f],
       &ones_mortar_gauss[stride],
       deg_mortar_gauss[f],
       (P4EST_DIM)-1,
       &VT_w_term1_mortar_lobatto[stride_lobatto]
      );
    
    for (int d = 0; d < (P4EST_DIM); d++){
      d4est_operators_apply_curvedgaussMass_ongaussNodeVec
        (
         d4est_ops,
         &term2_mortar_gauss[d][stride],
         deg_mortar_lobatto[f],
         &ones_mortar_gauss[stride],
         deg_mortar_gauss[f],
         (P4EST_DIM)-1,
         &VT_w_term2_mortar_lobatto[d][stride_lobatto]
        );
    }
      
    d4est_operators_apply_curvedgaussMass_ongaussNodeVec
      (
       d4est_ops,
       &term3_mortar_gauss[stride],
       deg_mortar_lobatto[f],
       &ones_mortar_gauss[stride],
       deg_mortar_gauss[f],
       (P4EST_DIM)-1,
       &VT_w_term3_mortar_lobatto[stride_lobatto]
      );

    stride += nodes_mortar_gauss[f];
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
 
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_p_on_f_p);

  D4EST_FREE_DBYD_MAT(drst_dxyz_m_on_mortar_gauss);
  D4EST_FREE_DBYD_MAT(drst_dxyz_p_on_mortar_gauss_porder);

  P4EST_FREE(sigma);
  P4EST_FREE(u_m_on_f_m_mortar);
  P4EST_FREE(u_m_on_f_m_mortar_gauss);
  P4EST_FREE(u_p_on_f_p_mortar);
  P4EST_FREE(u_p_on_f_p_mortar_gauss);
  for (int i = 0; i < (P4EST_DIM); i++){   
    P4EST_FREE(dudr_p_on_f_p_porder[i]);
    P4EST_FREE(dudr_p_on_f_p_mortar_porder[i]);
    P4EST_FREE(dudr_p_on_f_p_mortar_gauss_porder[i]);
    P4EST_FREE(dudx_p_on_f_p_mortar_gauss_porder[i]);
    P4EST_FREE(dudx_p_on_f_p_mortar_gauss[i]);
    P4EST_FREE(dudr_m_on_f_m[i]);
    P4EST_FREE(dudr_m_on_f_m_mortar[i]);
    P4EST_FREE(dudr_m_on_f_m_mortar_gauss[i]);
    P4EST_FREE(dudx_m_on_f_m_mortar_gauss[i]);
  }
  P4EST_FREE(sj_on_f_m_mortar_gauss);
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(n_on_f_m_mortar_gauss[i]);
    P4EST_FREE(n_sj_on_f_m_mortar_gauss[i]);
  }

  P4EST_FREE(j_div_sj_on_f_m_mortar_gauss);
  P4EST_FREE(j_div_sj_on_f_p_mortar_gauss_porder);
  P4EST_FREE(j_div_sj_on_f_p_mortar_gauss_porder_oriented);
  
  P4EST_FREE(ones_mortar_gauss);
  P4EST_FREE(lifted_proj_VT_w_term1_mortar_lobatto);
  P4EST_FREE(proj_VT_w_term1_mortar_lobatto);
  P4EST_FREE(VT_w_term1_mortar_lobatto);
  P4EST_FREE(term1_mortar_gauss);

  D4EST_FREE_DIM_VEC(DT_lifted_proj_VT_w_term2_mortar_lobatto);
  D4EST_FREE_DIM_VEC(lifted_proj_VT_w_term2_mortar_lobatto);
  D4EST_FREE_DIM_VEC(proj_VT_w_term2_mortar_lobatto);
  D4EST_FREE_DIM_VEC(VT_w_term2_mortar_lobatto);
  D4EST_FREE_DIM_VEC(term2_mortar_gauss);

  P4EST_FREE(lifted_proj_VT_w_term3_mortar_lobatto);
  P4EST_FREE(proj_VT_w_term3_mortar_lobatto);
  P4EST_FREE(VT_w_term3_mortar_lobatto);
  P4EST_FREE(term3_mortar_gauss);
  
  P4EST_FREE(tmp);  
}

curved_flux_fcn_ptrs_t
curved_gauss_primal_sipg_hesthaven_flux_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 ip_flux_params_t* curved_gauss_sipg_params
)
{  
  curved_flux_fcn_ptrs_t curved_gauss_primal_sipg_hesthaven_flux_fcns;
  curved_gauss_primal_sipg_hesthaven_flux_fcns.flux_interface_fcn = curved_gauss_primal_sipg_hesthaven_flux_interface;
  curved_gauss_primal_sipg_hesthaven_flux_fcns.flux_boundary_fcn = curved_gauss_primal_sipg_hesthaven_flux_dirichlet_withbcterms;
  curved_gauss_primal_sipg_hesthaven_flux_fcns.bndry_fcn = bndry_fcn;
  curved_gauss_primal_sipg_hesthaven_flux_fcns.params = (void*)curved_gauss_sipg_params;
  
  return curved_gauss_primal_sipg_hesthaven_flux_fcns;
}

/* curved_flux_fcn_ptrs_t */
/* curved_gauss_primal_sipg_hesthaven_flux_dirichlet_fetch_fcns_testrhsbc */
/* ( */
/*  grid_fcn_t bndry_fcn, */
/*  ip_flux_params_t* curved_gauss_sipg_params */
/* ) */
/* {   */
/*   curved_flux_fcn_ptrs_t curved_gauss_primal_sipg_hesthaven_flux_fcns; */
/*   curved_gauss_primal_sipg_hesthaven_flux_fcns.flux_interface_fcn = curved_gauss_primal_sipg_hesthaven_flux_interface_testrhsbc; */
/*   curved_gauss_primal_sipg_hesthaven_flux_fcns.flux_boundary_fcn = curved_gauss_primal_sipg_hesthaven_flux_dirichlet_testrhsbc; */
/*   curved_gauss_primal_sipg_hesthaven_flux_fcns.bndry_fcn = bndry_fcn; */
/*   curved_gauss_primal_sipg_hesthaven_flux_fcns.params = (void*)curved_gauss_sipg_params; */
  
/*   return curved_gauss_primal_sipg_hesthaven_flux_fcns; */
/* } */
