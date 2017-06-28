#include "../Utilities/util.h"
#include "../dGMath/d4est_operators.h"
#include "../ElementData/d4est_element_data.h"
#include "../LinearAlgebra/d4est_linalg.h"
#include "../Flux/curved_Gauss_primal_sipg_kronbichler_flux_fcns.h"



static void
curved_gauss_primal_sipg_kronbichler_uniform_flux_dirichlet
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
  double sipg_kronbichler_uniform_flux_penalty_prefactor = ip_flux_params->ip_flux_penalty_prefactor;
  penalty_calc_t sipg_kronbichler_uniform_flux_penalty_calculate_fcn = ip_flux_params->ip_flux_penalty_calculate_fcn;

  
  grid_fcn_t u_at_bndry = bndry_fcn;
  int face_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_gauss = d4est_lgl_get_nodes((P4EST_DIM) - 1, e_m->deg_quad);

  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* u_m_on_f_m_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);

  double* dudr_m_on_f_m [(P4EST_DIM)];
  double* dudr_m_on_f_m_gauss [(P4EST_DIM)];
  double* dudx_m_on_f_m_gauss [(P4EST_DIM)];
  
  double* u_m_on_f_m_min_u_at_bndry_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  double* u_m_on_f_m_min_u_at_bndry_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* sj_on_f_m_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  double* n_on_f_m_gauss [(P4EST_DIM)];
  double* n_sj_on_f_m_gauss [(P4EST_DIM)];

  double* drst_dxyz_gauss [(P4EST_DIM)][(P4EST_DIM)]; D4EST_ALLOC_DBYD_MAT(drst_dxyz_gauss, face_nodes_m_gauss);
  
  for (int d = 0; d < (P4EST_DIM); d++) {
    dudr_m_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m_lobatto);
    dudr_m_on_f_m_gauss[d] = P4EST_ALLOC(double, face_nodes_m_gauss);
    dudx_m_on_f_m_gauss[d] = P4EST_ALLOC(double, face_nodes_m_gauss);
    n_on_f_m_gauss[d] = P4EST_ALLOC(double, face_nodes_m_gauss);
  }

  double* u_at_bndry_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* u_at_bndry_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);


  int volume_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM), e_m->deg);
  
  double* ones_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  d4est_linalg_fill_vec(ones_gauss, 1., face_nodes_m_gauss);
  
  double* term1_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  double* VT_w_term1_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* lifted_VT_w_term1_lobatto = P4EST_ALLOC(double, volume_nodes_m_lobatto);

  double* term2_gauss [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(term2_gauss, face_nodes_m_gauss);
  double* VT_w_term2_lobatto [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(VT_w_term2_lobatto, face_nodes_m_lobatto);
  double* lifted_VT_w_term2_lobatto [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(lifted_VT_w_term2_lobatto, volume_nodes_m_lobatto);
  double* DT_lifted_VT_w_term2_lobatto [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(DT_lifted_VT_w_term2_lobatto, volume_nodes_m_lobatto);
 
  double* term3_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  double* VT_w_term3_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* lifted_VT_w_term3_lobatto = P4EST_ALLOC(double, volume_nodes_m_lobatto);

  double* J_div_SJ_gauss = NULL;
  if (ip_flux_params->ip_flux_h_calc == H_EQ_J_DIV_SJ){
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

    sigma[i] = sipg_kronbichler_uniform_flux_penalty_calculate_fcn
               (
                e_m->deg,
                h,
                e_m->deg,
                h,
                sipg_kronbichler_uniform_flux_penalty_prefactor
               );
  }

  if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ){
    P4EST_FREE(J_div_SJ_gauss);
  }

  
  
  
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

    D4EST_FREE_DIM_VEC(xyz_on_f_m);
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

      //      printf("boundary = %f\n", u_at_bndry_gauss[i]);
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
     
    term1_gauss[i] = 0.;
    for (int d = 0; d < (P4EST_DIM); d++){
      /* term1_gauss[i] += -1.*sj_on_f_m_gauss[i] */
                        /* *n_on_f_m_gauss[d][i] */
                        /* *(dudx_m_on_f_m_gauss[d][i]); */
      term1_gauss[i] += -1.*n_sj_on_f_m_gauss[d][i]
                        *(dudx_m_on_f_m_gauss[d][i]);

      
    }
    
    for (int l = 0; l < (P4EST_DIM); l++){
      term2_gauss[l][i] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        /* term2_gauss[l][i] += -.5*sj_on_f_m_gauss[i] */
                             /* *drst_dxyz_gauss[l][d][i] */
                             /* *n_on_f_m_gauss[d][i] */
                             /* *2.*u_m_on_f_m_min_u_at_bndry_gauss[i]; */
        term2_gauss[l][i] += -.5*drst_dxyz_gauss[l][d][i]
                             *n_sj_on_f_m_gauss[d][i]
                             *2.*u_m_on_f_m_min_u_at_bndry_gauss[i];
      }
    }

    term3_gauss[i] = sj_on_f_m_gauss[i]
                     *sigma[i]
                     *2.*u_m_on_f_m_min_u_at_bndry_gauss[i];
  }

  d4est_operators_apply_curvedgaussMass_ongaussNodeVec
    (
     d4est_ops,
     term1_gauss,
     e_m->deg,
     ones_gauss,
     e_m->deg_quad,
     (P4EST_DIM)-1,
     VT_w_term1_lobatto
    );

  
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


      
      
  int volume_nodes_m = d4est_lgl_get_nodes((P4EST_DIM), e_m->deg);
  for (int i = 0; i < volume_nodes_m; i++){
    for (int d = 0; d < (P4EST_DIM); d++){
      e_m->Au_elem[i] += DT_lifted_VT_w_term2_lobatto[d][i];
    }
    e_m->Au_elem[i] += lifted_VT_w_term3_lobatto[i];
    e_m->Au_elem[i] += lifted_VT_w_term1_lobatto[i];
  }
  
  P4EST_FREE(u_at_bndry_lobatto);
  P4EST_FREE(u_at_bndry_gauss);  
  P4EST_FREE(sigma);
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_m_on_f_m_gauss);
  for (int i = 0; i < (P4EST_DIM); i++) {
    P4EST_FREE(dudr_m_on_f_m[i]);
    P4EST_FREE(dudr_m_on_f_m_gauss[i]);
    P4EST_FREE(dudx_m_on_f_m_gauss[i]);
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
  P4EST_FREE(term1_gauss);
  P4EST_FREE(VT_w_term1_lobatto);
  P4EST_FREE(lifted_VT_w_term1_lobatto);
  D4EST_FREE_DIM_VEC(term2_gauss)
    D4EST_FREE_DIM_VEC(VT_w_term2_lobatto)
    D4EST_FREE_DIM_VEC(lifted_VT_w_term2_lobatto)
    D4EST_FREE_DIM_VEC(DT_lifted_VT_w_term2_lobatto)
    P4EST_FREE(term3_gauss);
  P4EST_FREE(VT_w_term3_lobatto);
  P4EST_FREE(lifted_VT_w_term3_lobatto);
}

static void
curved_gauss_primal_sipg_kronbichler_uniform_flux_interface
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
  /* assume uniform refinement */
  mpi_assert(faces_m == faces_p);
  mpi_assert(e_m[0]->deg_quad == e_p[0]->deg_quad);
  mpi_assert(e_m[0]->deg == e_p[0]->deg);
  
  ip_flux_params_t* ip_flux_params = (ip_flux_params_t*) params;
  double sipg_kronbichler_uniform_flux_penalty_prefactor = ip_flux_params->ip_flux_penalty_prefactor;
  penalty_calc_t sipg_kronbichler_uniform_flux_penalty_calculate_fcn = ip_flux_params->ip_flux_penalty_calculate_fcn;

  int face_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM)-1, e_m[0]->deg);
  int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM), e_m[0]->deg);
  int face_nodes_gauss = d4est_lgl_get_nodes((P4EST_DIM)-1, e_m[0]->deg_quad);
  int volume_nodes_gauss = d4est_lgl_get_nodes((P4EST_DIM), e_m[0]->deg_quad);
  
  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_lobatto);
  double* u_m_on_f_m_gauss = P4EST_ALLOC(double, face_nodes_gauss);
  double* u_p_on_f_p = P4EST_ALLOC(double, face_nodes_lobatto);
  double* u_p_on_f_p_gauss = P4EST_ALLOC(double, face_nodes_gauss);
  double* tmp = P4EST_ALLOC(double, face_nodes_lobatto);
  double* ones_gauss = P4EST_ALLOC(double, face_nodes_gauss);
  double* term1_gauss = P4EST_ALLOC(double, face_nodes_gauss);
  double* VT_w_term1_lobatto = P4EST_ALLOC(double, face_nodes_lobatto);
  double* lifted_VT_w_term1_lobatto = P4EST_ALLOC(double, volume_nodes_lobatto);
  double* term3_gauss = P4EST_ALLOC(double, face_nodes_gauss);
  double* VT_w_term3_lobatto = P4EST_ALLOC(double, face_nodes_lobatto);
  double* lifted_VT_w_term3_lobatto = P4EST_ALLOC(double, volume_nodes_lobatto);
  double* sj_on_f_m_gauss = P4EST_ALLOC(double, face_nodes_gauss);  

  double* dudr_m_on_f_m [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudr_m_on_f_m, face_nodes_lobatto);
  double* dudr_p_on_f_p_porder [(P4EST_DIM)]; 
  D4EST_ALLOC_DIM_VEC(dudr_p_on_f_p_porder, face_nodes_lobatto);
  double* VT_w_term2_lobatto [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(VT_w_term2_lobatto, face_nodes_lobatto);
  double* lifted_VT_w_term2_lobatto [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(lifted_VT_w_term2_lobatto, volume_nodes_lobatto);
  double* DT_lifted_VT_w_term2_lobatto [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(DT_lifted_VT_w_term2_lobatto, volume_nodes_lobatto);

  double* dudr_m_on_f_m_gauss [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudr_m_on_f_m_gauss, face_nodes_gauss);
  double* dudr_p_on_f_p_gauss_porder [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudr_p_on_f_p_gauss_porder, face_nodes_gauss);
  double* n_on_f_m_gauss [(P4EST_DIM)];
  double* n_sj_on_f_m_gauss [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(n_on_f_m_gauss, face_nodes_gauss);
  D4EST_ALLOC_DIM_VEC(n_sj_on_f_m_gauss, face_nodes_gauss);
  double* dudx_m_on_f_m_gauss [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudx_m_on_f_m_gauss, face_nodes_gauss);
  double* dudx_p_on_f_p_gauss [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudx_p_on_f_p_gauss, face_nodes_gauss);
  double* dudx_p_on_f_p_gauss_porder [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudx_p_on_f_p_gauss_porder, face_nodes_gauss);
  double* term2_gauss [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(term2_gauss, face_nodes_gauss);
  
  double* drst_dxyz_m_on_gauss [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_m_on_gauss, face_nodes_gauss);
  double* drst_dxyz_p_on_gauss_porder [(P4EST_DIM)][(P4EST_DIM)];  
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_p_on_gauss_porder, face_nodes_gauss);

  
  d4est_linalg_fill_vec(ones_gauss, 1., face_nodes_gauss);

  d4est_operators_apply_slicer
    (
     d4est_ops,
     &(e_m[0]->u_elem[0]),
     (P4EST_DIM),
     f_m,
     e_m[0]->deg,
     u_m_on_f_m
    );

  d4est_operators_apply_slicer
    (
     d4est_ops,
     &(e_p[0]->u_elem[0]),
     (P4EST_DIM),
     f_p,
     e_p[0]->deg,
     tmp
    );
    
  d4est_operators_reorient_face_data
    (
     d4est_ops,
     tmp,
     ((P4EST_DIM) - 1),
     e_p[0]->deg,
     orientation,
     f_m,
     f_p,
     u_p_on_f_p
    );

  d4est_operators_interp_lobatto_to_GL(d4est_ops,
                          u_m_on_f_m,
                          e_m[0]->deg,
                          e_m[0]->deg_quad,
                          u_m_on_f_m_gauss,
                          (P4EST_DIM)-1);
    
  d4est_operators_interp_lobatto_to_GL(d4est_ops,
                          u_p_on_f_p,
                          e_p[0]->deg,
                          e_p[0]->deg_quad,
                          u_p_on_f_p_gauss,
                          (P4EST_DIM)-1);

  
  /* For each component of the vector */
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_slicer
      (
       d4est_ops,
       e_m[0]->dudr_elem[d],
       (P4EST_DIM),
       f_m,
       e_m[0]->deg,
       dudr_m_on_f_m[d]
      );
  }
    
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_slicer
      (
       d4est_ops,
       e_p[0]->dudr_elem[d],
       (P4EST_DIM),
       f_p,
       e_p[0]->deg,
       dudr_p_on_f_p_porder[d]
      );

    d4est_operators_interp_lobatto_to_GL
      (
       d4est_ops,
       dudr_m_on_f_m[d],
       e_m[0]->deg,
       e_m[0]->deg_quad,
       dudr_m_on_f_m_gauss[d],
       (P4EST_DIM)-1
      );
    
    d4est_operators_interp_lobatto_to_GL
      (
       d4est_ops,
       dudr_p_on_f_p_porder[d],
       e_p[0]->deg,
       e_p[0]->deg_quad,
       dudr_p_on_f_p_gauss_porder[d],
       (P4EST_DIM)-1
      );
  }

  double* j_div_sj_on_f_m_gauss = NULL;
  double* j_div_sj_on_f_p_gauss_porder = NULL;
  double* j_div_sj_on_f_p_gauss_porder_oriented = NULL;

  if(ip_flux_params->ip_flux_h_calc == H_EQ_J_DIV_SJ){
    j_div_sj_on_f_m_gauss = P4EST_ALLOC(double, face_nodes_gauss);
    j_div_sj_on_f_p_gauss_porder =  P4EST_ALLOC(double, face_nodes_gauss);
    j_div_sj_on_f_p_gauss_porder_oriented =  P4EST_ALLOC(double, face_nodes_gauss);
  }

  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     1,
     1,
     &e_m[0]->deg_quad,
     f_m,
     drst_dxyz_m_on_gauss,
     sj_on_f_m_gauss,
     n_on_f_m_gauss,
     n_sj_on_f_m_gauss,
     j_div_sj_on_f_m_gauss,
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
     1,
     1,
     &e_p[0]->deg_quad,
     f_p,
     drst_dxyz_p_on_gauss_porder,
     NULL,
     NULL,
     NULL,
     j_div_sj_on_f_p_gauss_porder,
     GAUSS,
     geom,
     d4est_ops,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_linalg_fill_vec
      (
       dudx_m_on_f_m_gauss[d],
       0.0,
       face_nodes_gauss
      );


    d4est_linalg_fill_vec
      (
       dudx_p_on_f_p_gauss_porder[d],
       0.0,
       face_nodes_gauss
      );
  }
  
  for (int j = 0; j < (P4EST_DIM); j++){
    for (int i = 0; i < (P4EST_DIM); i++){   
      for (int k = 0; k < face_nodes_gauss; k++){
        dudx_m_on_f_m_gauss[j][k] +=
          drst_dxyz_m_on_gauss[i][j][k]*dudr_m_on_f_m_gauss[i][k];
        dudx_p_on_f_p_gauss_porder[j][k] += drst_dxyz_p_on_gauss_porder[i][j][k]*dudr_p_on_f_p_gauss_porder[i][k];
      }
    }    
  }
  

  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_reorient_face_data
      (
       d4est_ops,
       dudx_p_on_f_p_gauss_porder[d],
       (P4EST_DIM)-1,
       e_p[0]->deg_quad,
       orientation,
       f_m,
       f_p,
       dudx_p_on_f_p_gauss[d]
      );
  }

  if (j_div_sj_on_f_p_gauss_porder != NULL){
    d4est_operators_reorient_face_data
      (
       d4est_ops,
       j_div_sj_on_f_p_gauss_porder,
       (P4EST_DIM)-1,
       e_p[0]->deg_quad,
       orientation,
       f_m,
       f_p,
       j_div_sj_on_f_p_gauss_porder_oriented
      );
  }    
    
  double* sigma = P4EST_ALLOC(double, face_nodes_gauss);
  for (int ks= 0; ks < face_nodes_gauss; ks++){
    if (ip_flux_params->ip_flux_h_calc == H_EQ_VOLUME_DIV_AREA){
      sigma[ks] = sipg_kronbichler_uniform_flux_penalty_calculate_fcn
                  (
                   e_m[0]->deg,
                   e_m[0]->volume/e_m[0]->surface_area[f_m],
                   e_p[0]->deg,
                   e_p[0]->volume/e_p[0]->surface_area[f_p],
                   sipg_kronbichler_uniform_flux_penalty_prefactor
                  );
    }
    else if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ){
      double hp = j_div_sj_on_f_p_gauss_porder_oriented[ks];
      double hm = j_div_sj_on_f_m_gauss[ks];

      sigma[ks] = sipg_kronbichler_uniform_flux_penalty_calculate_fcn
                  (
                   e_m[0]->deg,
                   hm,
                   e_p[0]->deg,
                   hp,
                   sipg_kronbichler_uniform_flux_penalty_prefactor
                  );

    }
    else {
      mpi_abort("Select j_DIV_SJ or VOLUME_DIV_AREA for ip_flux_h_calc");
    }
  }
  
  for (int ks = 0; ks < face_nodes_gauss; ks++){

    term1_gauss[ks] = 0.;
    for (int d = 0; d < (P4EST_DIM); d++){
      /* term1_gauss[ks] += -1.*sj_on_f_m_gauss[ks] */
                         /* *n_sj_on_f_m_gauss[d][ks] */
                         /* *.5*(dudx_p_on_f_p_gauss[d][ks] + dudx_m_on_f_m_gauss[d][ks]); */
      term1_gauss[ks] += -1.*n_sj_on_f_m_gauss[d][ks]
                         *.5*(dudx_p_on_f_p_gauss[d][ks] + dudx_m_on_f_m_gauss[d][ks]);
      
    }
        
    for (int l = 0; l < (P4EST_DIM); l++){
      term2_gauss[l][ks] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        /* term2_gauss[l][ks] += -.5*sj_on_f_m_gauss[ks] */
                              /* *drst_dxyz_m_on_gauss[l][d][ks] */
                              /* *n_on_f_m_gauss[d][ks] */
                              /* *(u_m_on_f_m_gauss[ks] - u_p_on_f_p_gauss[ks]); */
        term2_gauss[l][ks] += -.5*drst_dxyz_m_on_gauss[l][d][ks]
                              *n_sj_on_f_m_gauss[d][ks]
                              *(u_m_on_f_m_gauss[ks] - u_p_on_f_p_gauss[ks]);
        
      }
    }
        
    term3_gauss[ks] = sj_on_f_m_gauss[ks]*sigma[ks]*(u_m_on_f_m_gauss[ks] - u_p_on_f_p_gauss[ks]);
  }
              
  d4est_operators_apply_curvedgaussMass_ongaussNodeVec
    (
     d4est_ops,
     term1_gauss,
     e_m[0]->deg,
     ones_gauss,
     e_m[0]->deg_quad,
     (P4EST_DIM)-1,
     VT_w_term1_lobatto
    );
    
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_curvedgaussMass_ongaussNodeVec
      (
       d4est_ops,
       term2_gauss[d],
       e_m[0]->deg,
       ones_gauss,
       e_m[0]->deg_quad,
       (P4EST_DIM)-1,
       VT_w_term2_lobatto[d]
      );
  }
      
  d4est_operators_apply_curvedgaussMass_ongaussNodeVec
    (
     d4est_ops,
     term3_gauss,
     e_m[0]->deg,
     ones_gauss,
     e_m[0]->deg_quad,
     (P4EST_DIM)-1,
     VT_w_term3_lobatto
    );

  d4est_operators_apply_lift
    (
     d4est_ops,
     VT_w_term1_lobatto,
     (P4EST_DIM),
     e_m[0]->deg,
     f_m,
     lifted_VT_w_term1_lobatto
    );


  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_lift(
                      d4est_ops,
                      VT_w_term2_lobatto[d],
                      (P4EST_DIM),
                      e_m[0]->deg,
                      f_m,
                      lifted_VT_w_term2_lobatto[d]);

    d4est_operators_apply_dij_transpose
      (
       d4est_ops,
       lifted_VT_w_term2_lobatto[d],
       (P4EST_DIM),
       e_m[0]->deg,
       d,
       DT_lifted_VT_w_term2_lobatto[d]
      );

  }
        

  d4est_operators_apply_lift(
                    d4est_ops,
                    VT_w_term3_lobatto,
                    (P4EST_DIM),
                    e_m[0]->deg,
                    f_m,
                    lifted_VT_w_term3_lobatto);



        
  for (int i = 0; i < volume_nodes_lobatto; i++){
    for (int d = 0; d < (P4EST_DIM); d++){
      e_m[0]->Au_elem[i] += DT_lifted_VT_w_term2_lobatto[d][i];
    }
    e_m[0]->Au_elem[i] += lifted_VT_w_term3_lobatto[i];
    e_m[0]->Au_elem[i] += lifted_VT_w_term1_lobatto[i];
  }


  
  P4EST_FREE(sj_on_f_m_gauss);
  P4EST_FREE(lifted_VT_w_term3_lobatto);
  P4EST_FREE(VT_w_term3_lobatto);
  P4EST_FREE(term3_gauss);
  P4EST_FREE(lifted_VT_w_term1_lobatto);
  P4EST_FREE(VT_w_term1_lobatto);
  P4EST_FREE(term1_gauss);
  P4EST_FREE(sigma);
  P4EST_FREE(ones_gauss);
  P4EST_FREE(tmp);
  P4EST_FREE(u_p_on_f_p_gauss);
  P4EST_FREE(u_p_on_f_p);
  P4EST_FREE(u_m_on_f_m_gauss);
  P4EST_FREE(u_m_on_f_m);


  D4EST_FREE_DIM_VEC(dudr_m_on_f_m);
  D4EST_FREE_DIM_VEC(dudr_p_on_f_p_porder);
  D4EST_FREE_DIM_VEC(VT_w_term2_lobatto);
  D4EST_FREE_DIM_VEC(lifted_VT_w_term2_lobatto);
  D4EST_FREE_DIM_VEC(DT_lifted_VT_w_term2_lobatto);

  D4EST_FREE_DIM_VEC(dudr_m_on_f_m_gauss);
  D4EST_FREE_DIM_VEC(dudr_p_on_f_p_gauss_porder);
  D4EST_FREE_DIM_VEC(n_on_f_m_gauss);
  D4EST_FREE_DIM_VEC(dudx_m_on_f_m_gauss);
  D4EST_FREE_DIM_VEC(dudx_p_on_f_p_gauss);
  D4EST_FREE_DIM_VEC(dudx_p_on_f_p_gauss_porder);
  D4EST_FREE_DIM_VEC(term2_gauss);

  D4EST_FREE_DBYD_MAT(drst_dxyz_m_on_gauss);
  D4EST_FREE_DBYD_MAT(drst_dxyz_p_on_gauss_porder);
  
 
}

curved_flux_fcn_ptrs_t
curved_gauss_primal_sipg_kronbichler_uniform_flux_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 ip_flux_params_t* curved_gauss_sipg_params
)
{  
  curved_flux_fcn_ptrs_t curved_gauss_primal_sipg_kronbichler_uniform_flux_fcns;
  curved_gauss_primal_sipg_kronbichler_uniform_flux_fcns.flux_interface_fcn = curved_gauss_primal_sipg_kronbichler_uniform_flux_interface;
  curved_gauss_primal_sipg_kronbichler_uniform_flux_fcns.flux_boundary_fcn = curved_gauss_primal_sipg_kronbichler_uniform_flux_dirichlet;
  curved_gauss_primal_sipg_kronbichler_uniform_flux_fcns.bndry_fcn = bndry_fcn;
  curved_gauss_primal_sipg_kronbichler_uniform_flux_fcns.params = (void*)curved_gauss_sipg_params;
  
  return curved_gauss_primal_sipg_kronbichler_uniform_flux_fcns;
}
