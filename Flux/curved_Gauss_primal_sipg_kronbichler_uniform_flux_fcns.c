#include "../Utilities/util.h"
#include "../dGMath/d4est_operators.h"
#include "../ElementData/d4est_element_data.h"
#include "../LinearAlgebra/d4est_linalg.h"
#include "../Flux/curved_Gauss_primal_sipg_kronbichler_flux_fcns.h"



static void
curved_Gauss_primal_sipg_kronbichler_uniform_flux_dirichlet
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
  int face_nodes_m_Lobatto = d4est_operators_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_Gauss = d4est_operators_get_nodes((P4EST_DIM) - 1, e_m->deg_quad);

  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* u_m_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);

  double* dudr_m_on_f_m [(P4EST_DIM)];
  double* dudr_m_on_f_m_Gauss [(P4EST_DIM)];
  double* dudx_m_on_f_m_Gauss [(P4EST_DIM)];
  
  double* u_m_on_f_m_min_u_at_bndry_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* u_m_on_f_m_min_u_at_bndry_Lobatto = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* sj_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* n_on_f_m_Gauss [(P4EST_DIM)];
  double* n_sj_on_f_m_Gauss [(P4EST_DIM)];

  double* drst_dxyz_Gauss [(P4EST_DIM)][(P4EST_DIM)]; D4EST_ALLOC_DBYD_MAT(drst_dxyz_Gauss, face_nodes_m_Gauss);
  
  for (int d = 0; d < (P4EST_DIM); d++) {
    dudr_m_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m_Lobatto);
    dudr_m_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
    dudx_m_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
    n_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
  }

  double* u_at_bndry_Lobatto = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* u_at_bndry_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);


  int volume_nodes_m_Lobatto = d4est_operators_get_nodes((P4EST_DIM), e_m->deg);
  
  double* ones_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  d4est_linalg_fill_vec(ones_Gauss, 1., face_nodes_m_Gauss);
  
  double* term1_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* VT_w_term1_Lobatto = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* lifted_VT_w_term1_Lobatto = P4EST_ALLOC(double, volume_nodes_m_Lobatto);

  double* term2_Gauss [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(term2_Gauss, face_nodes_m_Gauss);
  double* VT_w_term2_Lobatto [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(VT_w_term2_Lobatto, face_nodes_m_Lobatto);
  double* lifted_VT_w_term2_Lobatto [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(lifted_VT_w_term2_Lobatto, volume_nodes_m_Lobatto);
  double* DT_lifted_VT_w_term2_Lobatto [P4EST_DIM]; D4EST_ALLOC_DIM_VEC(DT_lifted_VT_w_term2_Lobatto, volume_nodes_m_Lobatto);
 
  double* term3_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* VT_w_term3_Lobatto = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* lifted_VT_w_term3_Lobatto = P4EST_ALLOC(double, volume_nodes_m_Lobatto);

  double* J_div_SJ_Gauss = NULL;
  if (ip_flux_params->ip_flux_h_calc == H_EQ_J_DIV_SJ){
    J_div_SJ_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
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
     drst_dxyz_Gauss,
     sj_on_f_m_Gauss,
     n_on_f_m_Gauss,
     n_sj_on_f_m_Gauss,
     J_div_SJ_Gauss,
     GAUSS,
     geom,
     d4est_ops,
     COMPUTE_NORMAL_USING_JACOBIAN
    );
  
  double* sigma = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double h;
  for (int i = 0; i < face_nodes_m_Gauss; i++){
    if (ip_flux_params->ip_flux_h_calc == H_EQ_VOLUME_DIV_AREA){
      h = (e_m->volume/e_m->surface_area[f_m]);
    }
    else if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ){
      h = J_div_SJ_Gauss[i];
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
    P4EST_FREE(J_div_SJ_Gauss);
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
    
    d4est_operators_interp_GLL_to_GL
      (
       d4est_ops,
       dudr_m_on_f_m[d],
       e_m->deg,
       e_m->deg_quad,
       dudr_m_on_f_m_Gauss[d],
       (P4EST_DIM)-1
      );
   
  }


  
  d4est_operators_apply_slicer(d4est_ops, e_m->u_elem, (P4EST_DIM), f_m, e_m->deg, u_m_on_f_m);
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_linalg_fill_vec
      (
       dudx_m_on_f_m_Gauss[d],
       0.0,
       face_nodes_m_Gauss
      );
  }
  
  for (int j = 0; j < (P4EST_DIM); j++){
    for (int i = 0; i < (P4EST_DIM); i++){   
      for (int k = 0; k < face_nodes_m_Gauss; k++){
        dudx_m_on_f_m_Gauss[j][k] += drst_dxyz_Gauss[i][j][k]*dudr_m_on_f_m_Gauss[i][k];
      }
    }
  }


  if (ip_flux_params->ip_flux_bc_eval == BC_EVAL_ON_LOBATTO_POINTS){

    double* xyz_on_f_m [(P4EST_DIM)];
    D4EST_ALLOC_DIM_VEC(xyz_on_f_m, face_nodes_m_Lobatto);
    
    for (int d = 0; d < (P4EST_DIM); d++){

      d4est_operators_apply_slicer(d4est_ops,
                          e_m->xyz[d],
                          (P4EST_DIM),
                          f_m,
                          e_m->deg,
                          xyz_on_f_m[d]);

    }

    
    for (int i = 0; i < face_nodes_m_Lobatto; i++){
      u_at_bndry_Lobatto[i] = u_at_bndry
                              (
                               xyz_on_f_m[0][i],
                               xyz_on_f_m[1][i]
#if (P4EST_DIM)==3
                               ,
                               xyz_on_f_m[2][i]
#endif
                              );
      u_m_on_f_m_min_u_at_bndry_Lobatto[i] = u_m_on_f_m[i]
                                             - u_at_bndry_Lobatto[i];
    }
    
    d4est_operators_interp_GLL_to_GL
      (
       d4est_ops,
       u_m_on_f_m_min_u_at_bndry_Lobatto,
       e_m->deg,
       e_m->deg_quad,
       u_m_on_f_m_min_u_at_bndry_Gauss,
       (P4EST_DIM)-1
      );

    d4est_operators_interp_GLL_to_GL
      (
       d4est_ops,
       u_at_bndry_Lobatto,
       e_m->deg,
       e_m->deg_quad,
       u_at_bndry_Gauss,
       (P4EST_DIM)-1
      );

    D4EST_FREE_DIM_VEC(xyz_on_f_m);
  }
  else if (ip_flux_params->ip_flux_bc_eval == BC_EVAL_ON_GAUSS_POINTS){

    double* xyz_on_f_m_Gauss [(P4EST_DIM)];
    D4EST_ALLOC_DIM_VEC(xyz_on_f_m_Gauss, face_nodes_m_Gauss);
    
    d4est_operators_interp_GLL_to_GL
      (
       d4est_ops,
       u_m_on_f_m,
       e_m->deg,
       e_m->deg_quad,
       u_m_on_f_m_Gauss,
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
       xyz_on_f_m_Gauss
      );

    
    for (int i = 0; i < face_nodes_m_Gauss; i++){
      u_at_bndry_Gauss[i]
        = u_at_bndry
        (
         xyz_on_f_m_Gauss[0][i],
         xyz_on_f_m_Gauss[1][i]
#if (P4EST_DIM)==3
         ,
         xyz_on_f_m_Gauss[2][i]
#endif
        );

      //      printf("boundary = %f\n", u_at_bndry_Gauss[i]);
    }
    
    for (int i = 0; i < face_nodes_m_Gauss; i++){
      u_m_on_f_m_min_u_at_bndry_Gauss[i] = u_m_on_f_m_Gauss[i]
                                           - u_at_bndry_Gauss[i];
    }      


    D4EST_FREE_DIM_VEC(xyz_on_f_m_Gauss);
  }
  else {
    mpi_abort("Select either BC_EVAL_ON_GAUSS_POINTS or BC_EVAL_ON_LOBATTO_POINTS");
  }
  
  for(int i = 0; i < face_nodes_m_Gauss; i++){
     
    term1_Gauss[i] = 0.;
    for (int d = 0; d < (P4EST_DIM); d++){
      /* term1_Gauss[i] += -1.*sj_on_f_m_Gauss[i] */
                        /* *n_on_f_m_Gauss[d][i] */
                        /* *(dudx_m_on_f_m_Gauss[d][i]); */
      term1_Gauss[i] += -1.*n_sj_on_f_m_Gauss[d][i]
                        *(dudx_m_on_f_m_Gauss[d][i]);

      
    }
    
    for (int l = 0; l < (P4EST_DIM); l++){
      term2_Gauss[l][i] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        /* term2_Gauss[l][i] += -.5*sj_on_f_m_Gauss[i] */
                             /* *drst_dxyz_Gauss[l][d][i] */
                             /* *n_on_f_m_Gauss[d][i] */
                             /* *2.*u_m_on_f_m_min_u_at_bndry_Gauss[i]; */
        term2_Gauss[l][i] += -.5*drst_dxyz_Gauss[l][d][i]
                             *n_sj_on_f_m_Gauss[d][i]
                             *2.*u_m_on_f_m_min_u_at_bndry_Gauss[i];
      }
    }

    term3_Gauss[i] = sj_on_f_m_Gauss[i]
                     *sigma[i]
                     *2.*u_m_on_f_m_min_u_at_bndry_Gauss[i];
  }

  d4est_operators_apply_curvedGaussMass_onGaussNodeVec
    (
     d4est_ops,
     term1_Gauss,
     e_m->deg,
     ones_Gauss,
     e_m->deg_quad,
     (P4EST_DIM)-1,
     VT_w_term1_Lobatto
    );

  
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_curvedGaussMass_onGaussNodeVec
      (
       d4est_ops,
       term2_Gauss[d],
       e_m->deg,
       ones_Gauss,
       e_m->deg_quad,
       (P4EST_DIM)-1,
       VT_w_term2_Lobatto[d]
      );
  }

  
  d4est_operators_apply_curvedGaussMass_onGaussNodeVec
    (
     d4est_ops,
     term3_Gauss,
     e_m->deg,
     ones_Gauss,
     e_m->deg_quad,
     (P4EST_DIM)-1,
     VT_w_term3_Lobatto
    );  
  
  d4est_operators_apply_LIFT(
                    d4est_ops,
                    VT_w_term1_Lobatto,
                    (P4EST_DIM),
                    e_m->deg,
                    f_m,
                    lifted_VT_w_term1_Lobatto);


  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_LIFT(
                      d4est_ops,
                      VT_w_term2_Lobatto[d],
                      (P4EST_DIM),
                      e_m->deg,
                      f_m,
                      lifted_VT_w_term2_Lobatto[d]);

    d4est_operators_apply_Dij_transpose(d4est_ops,
                               lifted_VT_w_term2_Lobatto[d],
                               (P4EST_DIM),
                               e_m->deg,
                               d,
                               DT_lifted_VT_w_term2_Lobatto[d]
                              );

  }
        

  d4est_operators_apply_LIFT(
                    d4est_ops,
                    VT_w_term3_Lobatto,
                    (P4EST_DIM),
                    e_m->deg,
                    f_m,
                    lifted_VT_w_term3_Lobatto);


      
      
  int volume_nodes_m = d4est_operators_get_nodes((P4EST_DIM), e_m->deg);
  for (int i = 0; i < volume_nodes_m; i++){
    for (int d = 0; d < (P4EST_DIM); d++){
      e_m->Au_elem[i] += DT_lifted_VT_w_term2_Lobatto[d][i];
    }
    e_m->Au_elem[i] += lifted_VT_w_term3_Lobatto[i];
    e_m->Au_elem[i] += lifted_VT_w_term1_Lobatto[i];
  }
  
  P4EST_FREE(u_at_bndry_Lobatto);
  P4EST_FREE(u_at_bndry_Gauss);  
  P4EST_FREE(sigma);
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_m_on_f_m_Gauss);
  for (int i = 0; i < (P4EST_DIM); i++) {
    P4EST_FREE(dudr_m_on_f_m[i]);
    P4EST_FREE(dudr_m_on_f_m_Gauss[i]);
    P4EST_FREE(dudx_m_on_f_m_Gauss[i]);
  }

  P4EST_FREE(sj_on_f_m_Gauss);
  P4EST_FREE(u_m_on_f_m_min_u_at_bndry_Gauss);
  P4EST_FREE(u_m_on_f_m_min_u_at_bndry_Lobatto);
  
  for (int d = 0; d < (P4EST_DIM); d++){
    P4EST_FREE(n_on_f_m_Gauss[d]);
    P4EST_FREE(n_sj_on_f_m_Gauss[d]);
  }
  D4EST_FREE_DBYD_MAT(drst_dxyz_Gauss);

  P4EST_FREE(ones_Gauss);
  P4EST_FREE(term1_Gauss);
  P4EST_FREE(VT_w_term1_Lobatto);
  P4EST_FREE(lifted_VT_w_term1_Lobatto);
  D4EST_FREE_DIM_VEC(term2_Gauss)
    D4EST_FREE_DIM_VEC(VT_w_term2_Lobatto)
    D4EST_FREE_DIM_VEC(lifted_VT_w_term2_Lobatto)
    D4EST_FREE_DIM_VEC(DT_lifted_VT_w_term2_Lobatto)
    P4EST_FREE(term3_Gauss);
  P4EST_FREE(VT_w_term3_Lobatto);
  P4EST_FREE(lifted_VT_w_term3_Lobatto);
}

static void
curved_Gauss_primal_sipg_kronbichler_uniform_flux_interface
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

  int face_nodes_Lobatto = d4est_operators_get_nodes((P4EST_DIM)-1, e_m[0]->deg);
  int volume_nodes_Lobatto = d4est_operators_get_nodes((P4EST_DIM), e_m[0]->deg);
  int face_nodes_Gauss = d4est_operators_get_nodes((P4EST_DIM)-1, e_m[0]->deg_quad);
  int volume_nodes_Gauss = d4est_operators_get_nodes((P4EST_DIM), e_m[0]->deg_quad);
  
  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_Lobatto);
  double* u_m_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_Gauss);
  double* u_p_on_f_p = P4EST_ALLOC(double, face_nodes_Lobatto);
  double* u_p_on_f_p_Gauss = P4EST_ALLOC(double, face_nodes_Gauss);
  double* tmp = P4EST_ALLOC(double, face_nodes_Lobatto);
  double* ones_Gauss = P4EST_ALLOC(double, face_nodes_Gauss);
  double* term1_Gauss = P4EST_ALLOC(double, face_nodes_Gauss);
  double* VT_w_term1_Lobatto = P4EST_ALLOC(double, face_nodes_Lobatto);
  double* lifted_VT_w_term1_Lobatto = P4EST_ALLOC(double, volume_nodes_Lobatto);
  double* term3_Gauss = P4EST_ALLOC(double, face_nodes_Gauss);
  double* VT_w_term3_Lobatto = P4EST_ALLOC(double, face_nodes_Lobatto);
  double* lifted_VT_w_term3_Lobatto = P4EST_ALLOC(double, volume_nodes_Lobatto);
  double* sj_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_Gauss);  

  double* dudr_m_on_f_m [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudr_m_on_f_m, face_nodes_Lobatto);
  double* dudr_p_on_f_p_porder [(P4EST_DIM)]; 
  D4EST_ALLOC_DIM_VEC(dudr_p_on_f_p_porder, face_nodes_Lobatto);
  double* VT_w_term2_Lobatto [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(VT_w_term2_Lobatto, face_nodes_Lobatto);
  double* lifted_VT_w_term2_Lobatto [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(lifted_VT_w_term2_Lobatto, volume_nodes_Lobatto);
  double* DT_lifted_VT_w_term2_Lobatto [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(DT_lifted_VT_w_term2_Lobatto, volume_nodes_Lobatto);

  double* dudr_m_on_f_m_Gauss [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudr_m_on_f_m_Gauss, face_nodes_Gauss);
  double* dudr_p_on_f_p_Gauss_porder [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudr_p_on_f_p_Gauss_porder, face_nodes_Gauss);
  double* n_on_f_m_Gauss [(P4EST_DIM)];
  double* n_sj_on_f_m_Gauss [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(n_on_f_m_Gauss, face_nodes_Gauss);
  D4EST_ALLOC_DIM_VEC(n_sj_on_f_m_Gauss, face_nodes_Gauss);
  double* dudx_m_on_f_m_Gauss [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudx_m_on_f_m_Gauss, face_nodes_Gauss);
  double* dudx_p_on_f_p_Gauss [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudx_p_on_f_p_Gauss, face_nodes_Gauss);
  double* dudx_p_on_f_p_Gauss_porder [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(dudx_p_on_f_p_Gauss_porder, face_nodes_Gauss);
  double* term2_Gauss [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(term2_Gauss, face_nodes_Gauss);
  
  double* drst_dxyz_m_on_Gauss [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_m_on_Gauss, face_nodes_Gauss);
  double* drst_dxyz_p_on_Gauss_porder [(P4EST_DIM)][(P4EST_DIM)];  
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_p_on_Gauss_porder, face_nodes_Gauss);

  
  d4est_linalg_fill_vec(ones_Gauss, 1., face_nodes_Gauss);

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

  d4est_operators_interp_GLL_to_GL(d4est_ops,
                          u_m_on_f_m,
                          e_m[0]->deg,
                          e_m[0]->deg_quad,
                          u_m_on_f_m_Gauss,
                          (P4EST_DIM)-1);
    
  d4est_operators_interp_GLL_to_GL(d4est_ops,
                          u_p_on_f_p,
                          e_p[0]->deg,
                          e_p[0]->deg_quad,
                          u_p_on_f_p_Gauss,
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

    d4est_operators_interp_GLL_to_GL
      (
       d4est_ops,
       dudr_m_on_f_m[d],
       e_m[0]->deg,
       e_m[0]->deg_quad,
       dudr_m_on_f_m_Gauss[d],
       (P4EST_DIM)-1
      );
    
    d4est_operators_interp_GLL_to_GL
      (
       d4est_ops,
       dudr_p_on_f_p_porder[d],
       e_p[0]->deg,
       e_p[0]->deg_quad,
       dudr_p_on_f_p_Gauss_porder[d],
       (P4EST_DIM)-1
      );
  }

  double* j_div_sj_on_f_m_Gauss = NULL;
  double* j_div_sj_on_f_p_Gauss_porder = NULL;
  double* j_div_sj_on_f_p_Gauss_porder_oriented = NULL;

  if(ip_flux_params->ip_flux_h_calc == H_EQ_J_DIV_SJ){
    j_div_sj_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_Gauss);
    j_div_sj_on_f_p_Gauss_porder =  P4EST_ALLOC(double, face_nodes_Gauss);
    j_div_sj_on_f_p_Gauss_porder_oriented =  P4EST_ALLOC(double, face_nodes_Gauss);
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
     drst_dxyz_m_on_Gauss,
     sj_on_f_m_Gauss,
     n_on_f_m_Gauss,
     n_sj_on_f_m_Gauss,
     j_div_sj_on_f_m_Gauss,
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
     drst_dxyz_p_on_Gauss_porder,
     NULL,
     NULL,
     NULL,
     j_div_sj_on_f_p_Gauss_porder,
     GAUSS,
     geom,
     d4est_ops,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_linalg_fill_vec
      (
       dudx_m_on_f_m_Gauss[d],
       0.0,
       face_nodes_Gauss
      );


    d4est_linalg_fill_vec
      (
       dudx_p_on_f_p_Gauss_porder[d],
       0.0,
       face_nodes_Gauss
      );
  }
  
  for (int j = 0; j < (P4EST_DIM); j++){
    for (int i = 0; i < (P4EST_DIM); i++){   
      for (int k = 0; k < face_nodes_Gauss; k++){
        dudx_m_on_f_m_Gauss[j][k] +=
          drst_dxyz_m_on_Gauss[i][j][k]*dudr_m_on_f_m_Gauss[i][k];
        dudx_p_on_f_p_Gauss_porder[j][k] += drst_dxyz_p_on_Gauss_porder[i][j][k]*dudr_p_on_f_p_Gauss_porder[i][k];
      }
    }    
  }
  

  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_reorient_face_data
      (
       d4est_ops,
       dudx_p_on_f_p_Gauss_porder[d],
       (P4EST_DIM)-1,
       e_p[0]->deg_quad,
       orientation,
       f_m,
       f_p,
       dudx_p_on_f_p_Gauss[d]
      );
  }

  if (j_div_sj_on_f_p_Gauss_porder != NULL){
    d4est_operators_reorient_face_data
      (
       d4est_ops,
       j_div_sj_on_f_p_Gauss_porder,
       (P4EST_DIM)-1,
       e_p[0]->deg_quad,
       orientation,
       f_m,
       f_p,
       j_div_sj_on_f_p_Gauss_porder_oriented
      );
  }    
    
  double* sigma = P4EST_ALLOC(double, face_nodes_Gauss);
  for (int ks= 0; ks < face_nodes_Gauss; ks++){
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
      double hp = j_div_sj_on_f_p_Gauss_porder_oriented[ks];
      double hm = j_div_sj_on_f_m_Gauss[ks];

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
  
  for (int ks = 0; ks < face_nodes_Gauss; ks++){

    term1_Gauss[ks] = 0.;
    for (int d = 0; d < (P4EST_DIM); d++){
      /* term1_Gauss[ks] += -1.*sj_on_f_m_Gauss[ks] */
                         /* *n_sj_on_f_m_Gauss[d][ks] */
                         /* *.5*(dudx_p_on_f_p_Gauss[d][ks] + dudx_m_on_f_m_Gauss[d][ks]); */
      term1_Gauss[ks] += -1.*n_sj_on_f_m_Gauss[d][ks]
                         *.5*(dudx_p_on_f_p_Gauss[d][ks] + dudx_m_on_f_m_Gauss[d][ks]);
      
    }
        
    for (int l = 0; l < (P4EST_DIM); l++){
      term2_Gauss[l][ks] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        /* term2_Gauss[l][ks] += -.5*sj_on_f_m_Gauss[ks] */
                              /* *drst_dxyz_m_on_Gauss[l][d][ks] */
                              /* *n_on_f_m_Gauss[d][ks] */
                              /* *(u_m_on_f_m_Gauss[ks] - u_p_on_f_p_Gauss[ks]); */
        term2_Gauss[l][ks] += -.5*drst_dxyz_m_on_Gauss[l][d][ks]
                              *n_sj_on_f_m_Gauss[d][ks]
                              *(u_m_on_f_m_Gauss[ks] - u_p_on_f_p_Gauss[ks]);
        
      }
    }
        
    term3_Gauss[ks] = sj_on_f_m_Gauss[ks]*sigma[ks]*(u_m_on_f_m_Gauss[ks] - u_p_on_f_p_Gauss[ks]);
  }
              
  d4est_operators_apply_curvedGaussMass_onGaussNodeVec
    (
     d4est_ops,
     term1_Gauss,
     e_m[0]->deg,
     ones_Gauss,
     e_m[0]->deg_quad,
     (P4EST_DIM)-1,
     VT_w_term1_Lobatto
    );
    
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_curvedGaussMass_onGaussNodeVec
      (
       d4est_ops,
       term2_Gauss[d],
       e_m[0]->deg,
       ones_Gauss,
       e_m[0]->deg_quad,
       (P4EST_DIM)-1,
       VT_w_term2_Lobatto[d]
      );
  }
      
  d4est_operators_apply_curvedGaussMass_onGaussNodeVec
    (
     d4est_ops,
     term3_Gauss,
     e_m[0]->deg,
     ones_Gauss,
     e_m[0]->deg_quad,
     (P4EST_DIM)-1,
     VT_w_term3_Lobatto
    );

  d4est_operators_apply_LIFT
    (
     d4est_ops,
     VT_w_term1_Lobatto,
     (P4EST_DIM),
     e_m[0]->deg,
     f_m,
     lifted_VT_w_term1_Lobatto
    );


  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_LIFT(
                      d4est_ops,
                      VT_w_term2_Lobatto[d],
                      (P4EST_DIM),
                      e_m[0]->deg,
                      f_m,
                      lifted_VT_w_term2_Lobatto[d]);

    d4est_operators_apply_Dij_transpose
      (
       d4est_ops,
       lifted_VT_w_term2_Lobatto[d],
       (P4EST_DIM),
       e_m[0]->deg,
       d,
       DT_lifted_VT_w_term2_Lobatto[d]
      );

  }
        

  d4est_operators_apply_LIFT(
                    d4est_ops,
                    VT_w_term3_Lobatto,
                    (P4EST_DIM),
                    e_m[0]->deg,
                    f_m,
                    lifted_VT_w_term3_Lobatto);



        
  for (int i = 0; i < volume_nodes_Lobatto; i++){
    for (int d = 0; d < (P4EST_DIM); d++){
      e_m[0]->Au_elem[i] += DT_lifted_VT_w_term2_Lobatto[d][i];
    }
    e_m[0]->Au_elem[i] += lifted_VT_w_term3_Lobatto[i];
    e_m[0]->Au_elem[i] += lifted_VT_w_term1_Lobatto[i];
  }


  
  P4EST_FREE(sj_on_f_m_Gauss);
  P4EST_FREE(lifted_VT_w_term3_Lobatto);
  P4EST_FREE(VT_w_term3_Lobatto);
  P4EST_FREE(term3_Gauss);
  P4EST_FREE(lifted_VT_w_term1_Lobatto);
  P4EST_FREE(VT_w_term1_Lobatto);
  P4EST_FREE(term1_Gauss);
  P4EST_FREE(sigma);
  P4EST_FREE(ones_Gauss);
  P4EST_FREE(tmp);
  P4EST_FREE(u_p_on_f_p_Gauss);
  P4EST_FREE(u_p_on_f_p);
  P4EST_FREE(u_m_on_f_m_Gauss);
  P4EST_FREE(u_m_on_f_m);


  D4EST_FREE_DIM_VEC(dudr_m_on_f_m);
  D4EST_FREE_DIM_VEC(dudr_p_on_f_p_porder);
  D4EST_FREE_DIM_VEC(VT_w_term2_Lobatto);
  D4EST_FREE_DIM_VEC(lifted_VT_w_term2_Lobatto);
  D4EST_FREE_DIM_VEC(DT_lifted_VT_w_term2_Lobatto);

  D4EST_FREE_DIM_VEC(dudr_m_on_f_m_Gauss);
  D4EST_FREE_DIM_VEC(dudr_p_on_f_p_Gauss_porder);
  D4EST_FREE_DIM_VEC(n_on_f_m_Gauss);
  D4EST_FREE_DIM_VEC(dudx_m_on_f_m_Gauss);
  D4EST_FREE_DIM_VEC(dudx_p_on_f_p_Gauss);
  D4EST_FREE_DIM_VEC(dudx_p_on_f_p_Gauss_porder);
  D4EST_FREE_DIM_VEC(term2_Gauss);

  D4EST_FREE_DBYD_MAT(drst_dxyz_m_on_Gauss);
  D4EST_FREE_DBYD_MAT(drst_dxyz_p_on_Gauss_porder);
  
 
}

curved_flux_fcn_ptrs_t
curved_Gauss_primal_sipg_kronbichler_uniform_flux_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 ip_flux_params_t* curved_Gauss_sipg_params
)
{  
  curved_flux_fcn_ptrs_t curved_Gauss_primal_sipg_kronbichler_uniform_flux_fcns;
  curved_Gauss_primal_sipg_kronbichler_uniform_flux_fcns.flux_interface_fcn = curved_Gauss_primal_sipg_kronbichler_uniform_flux_interface;
  curved_Gauss_primal_sipg_kronbichler_uniform_flux_fcns.flux_boundary_fcn = curved_Gauss_primal_sipg_kronbichler_uniform_flux_dirichlet;
  curved_Gauss_primal_sipg_kronbichler_uniform_flux_fcns.bndry_fcn = bndry_fcn;
  curved_Gauss_primal_sipg_kronbichler_uniform_flux_fcns.params = (void*)curved_Gauss_sipg_params;
  
  return curved_Gauss_primal_sipg_kronbichler_uniform_flux_fcns;
}
