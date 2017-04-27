#include "../Utilities/util.h"
#include "../dGMath/dgmath.h"
#include "../ElementData/curved_element_data.h"
#include "../LinearAlgebra/linalg.h"
#include "../Flux/curved_Gauss_primal_sipg_kronbichler_flux_fcns.h"

static void
curved_Gauss_primal_sipg_kronbichler_flux_dirichlet
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
  int face_nodes_m_Lobatto = dgmath_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_Gauss = dgmath_get_nodes((P4EST_DIM) - 1, e_m->deg_integ);

  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* u_m_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);

  double* dudr_m_on_f_m [(P4EST_DIM)];
  double* dudr_m_on_f_m_Gauss [(P4EST_DIM)];
  double* dudx_m_on_f_m_Gauss [(P4EST_DIM)];
  
  double* u_m_on_f_m_min_u_at_bndry_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* u_m_on_f_m_min_u_at_bndry_Lobatto = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* sj_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* n_on_f_m_Gauss [(P4EST_DIM)];

  double* drst_dxyz_Gauss [(P4EST_DIM)][(P4EST_DIM)]; D4EST_ALLOC_DBYD_MAT(drst_dxyz_Gauss, face_nodes_m_Gauss);
  
  for (int d = 0; d < (P4EST_DIM); d++) {
    dudr_m_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m_Lobatto);
    dudr_m_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
    dudx_m_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
    n_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
  }

  double* u_at_bndry_Lobatto = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* u_at_bndry_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);


  int volume_nodes_m_Lobatto = dgmath_get_nodes((P4EST_DIM), e_m->deg);
  
  double* ones_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  linalg_fill_vec(ones_Gauss, 1., face_nodes_m_Gauss);
  
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
     &e_m->deg_integ,
     f_m,
     drst_dxyz_Gauss,
     sj_on_f_m_Gauss,
     n_on_f_m_Gauss,
     J_div_SJ_Gauss,
     GAUSS,
     geom,
     dgmath_jit_dbase,
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

    sigma[i] = sipg_kronbichler_flux_penalty_calculate_fcn
               (
                e_m->deg,
                h,
                e_m->deg,
                h,
                sipg_kronbichler_flux_penalty_prefactor
               );
  }

  if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ){
    P4EST_FREE(J_div_SJ_Gauss);
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
    
    dgmath_interp_GLL_to_GL
      (
       dgmath_jit_dbase,
       dudr_m_on_f_m[d],
       e_m->deg,
       e_m->deg_integ,
       dudr_m_on_f_m_Gauss[d],
       (P4EST_DIM)-1
      );
   
  }


  
  dgmath_apply_slicer(dgmath_jit_dbase, e_m->u_storage, (P4EST_DIM), f_m, e_m->deg, u_m_on_f_m);
  for (int d = 0; d < (P4EST_DIM); d++){
    linalg_fill_vec
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

      dgmath_apply_slicer(dgmath_jit_dbase,
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
    
    dgmath_interp_GLL_to_GL
      (
       dgmath_jit_dbase,
       u_m_on_f_m_min_u_at_bndry_Lobatto,
       e_m->deg,
       e_m->deg_integ,
       u_m_on_f_m_min_u_at_bndry_Gauss,
       (P4EST_DIM)-1
      );

    dgmath_interp_GLL_to_GL
      (
       dgmath_jit_dbase,
       u_at_bndry_Lobatto,
       e_m->deg,
       e_m->deg_integ,
       u_at_bndry_Gauss,
       (P4EST_DIM)-1
      );

    D4EST_FREE_DIM_VEC(xyz_on_f_m);
  }
  else if (ip_flux_params->ip_flux_bc_eval == BC_EVAL_ON_GAUSS_POINTS){

    double* xyz_on_f_m_Gauss [(P4EST_DIM)];
    D4EST_ALLOC_DIM_VEC(xyz_on_f_m_Gauss, face_nodes_m_Gauss);
    
    dgmath_interp_GLL_to_GL
      (
       dgmath_jit_dbase,
       u_m_on_f_m,
       e_m->deg,
       e_m->deg_integ,
       u_m_on_f_m_Gauss,
       (P4EST_DIM)-1
      );
    
    d4est_geometry_data_compute_xyz_face_analytic
      (
       dgmath_jit_dbase,
       e_m->q,
       e_m->dq,
       e_m->tree,
       f_m,
       geom,
       GAUSS,
       e_m->deg_integ,
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
      term1_Gauss[i] += -1.*sj_on_f_m_Gauss[i]
                        *n_on_f_m_Gauss[d][i]
                        *(dudx_m_on_f_m_Gauss[d][i]);
    }
    
    for (int l = 0; l < (P4EST_DIM); l++){
      term2_Gauss[l][i] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        term2_Gauss[l][i] += -.5*sj_on_f_m_Gauss[i]
                             *drst_dxyz_Gauss[l][d][i]
                             *n_on_f_m_Gauss[d][i]
                             *2.*u_m_on_f_m_min_u_at_bndry_Gauss[i];
      }
    }

    term3_Gauss[i] = sj_on_f_m_Gauss[i]
                     *sigma[i]
                     *2.*u_m_on_f_m_min_u_at_bndry_Gauss[i];
  }

  dgmath_apply_curvedGaussMass_onGaussNodeVec
    (
     dgmath_jit_dbase,
     term1_Gauss,
     e_m->deg,
     ones_Gauss,
     e_m->deg_integ,
     (P4EST_DIM)-1,
     VT_w_term1_Lobatto
    );

  
  for (int d = 0; d < (P4EST_DIM); d++){
    dgmath_apply_curvedGaussMass_onGaussNodeVec
      (
       dgmath_jit_dbase,
       term2_Gauss[d],
       e_m->deg,
       ones_Gauss,
       e_m->deg_integ,
       (P4EST_DIM)-1,
       VT_w_term2_Lobatto[d]
      );
  }

  
  dgmath_apply_curvedGaussMass_onGaussNodeVec
    (
     dgmath_jit_dbase,
     term3_Gauss,
     e_m->deg,
     ones_Gauss,
     e_m->deg_integ,
     (P4EST_DIM)-1,
     VT_w_term3_Lobatto
    );  
  
  dgmath_apply_LIFT(
                    dgmath_jit_dbase,
                    VT_w_term1_Lobatto,
                    (P4EST_DIM),
                    e_m->deg,
                    f_m,
                    lifted_VT_w_term1_Lobatto);


  for (int d = 0; d < (P4EST_DIM); d++){
    dgmath_apply_LIFT(
                      dgmath_jit_dbase,
                      VT_w_term2_Lobatto[d],
                      (P4EST_DIM),
                      e_m->deg,
                      f_m,
                      lifted_VT_w_term2_Lobatto[d]);

    dgmath_apply_Dij_transpose(dgmath_jit_dbase,
                               lifted_VT_w_term2_Lobatto[d],
                               (P4EST_DIM),
                               e_m->deg,
                               d,
                               DT_lifted_VT_w_term2_Lobatto[d]
                              );

  }
        

  dgmath_apply_LIFT(
                    dgmath_jit_dbase,
                    VT_w_term3_Lobatto,
                    (P4EST_DIM),
                    e_m->deg,
                    f_m,
                    lifted_VT_w_term3_Lobatto);


      
      
  int volume_nodes_m = dgmath_get_nodes((P4EST_DIM), e_m->deg);
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
curved_Gauss_primal_sipg_kronbichler_flux_interface
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
  int deg_p_Lobatto [(P4EST_HALF)];
  int face_nodes_p_Lobatto [(P4EST_HALF)];
  /* int deg_p_Gauss [(P4EST_HALF)]; */
  int face_nodes_p_Gauss [(P4EST_HALF)];

  int deg_m_Lobatto [(P4EST_HALF)];
  int face_nodes_m_Lobatto [(P4EST_HALF)];
  /* int deg_m_Gauss [(P4EST_HALF)]; */
  int face_nodes_m_Gauss [(P4EST_HALF)];
  
  int nodes_mortar_Gauss [(P4EST_HALF)];
  int nodes_mortar_Lobatto [(P4EST_HALF)];
  int deg_mortar_Gauss [(P4EST_HALF)];
  int deg_mortar_Lobatto [(P4EST_HALF)];
  int deg_p_Lobatto_porder [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
  double penalty_mortar [(P4EST_HALF)];

  curved_element_data_t* e_p_oriented [(P4EST_HALF)];
  curved_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);

  /* calculate degs and nodes of each face of (-) side */
  int max_volume_nodes_m_Lobatto = 0;
  int total_side_nodes_m_Lobatto = 0;
  int total_side_nodes_m_Gauss = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_Lobatto[i] = e_m[i]->deg;

    int volume_nodes_m_Lobatto = dgmath_get_nodes((P4EST_DIM), e_m[i]->deg);
    max_volume_nodes_m_Lobatto = (volume_nodes_m_Lobatto > max_volume_nodes_m_Lobatto) ? volume_nodes_m_Lobatto : max_volume_nodes_m_Lobatto;
    
    face_nodes_m_Lobatto[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg);
    face_nodes_m_Gauss[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg_integ);
    
    total_side_nodes_m_Lobatto += face_nodes_m_Lobatto[i];
    total_side_nodes_m_Gauss += face_nodes_m_Gauss[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_Lobatto = 0;
  int total_side_nodes_p_Gauss = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p_Lobatto[i] = e_p_oriented[i]->deg;
    deg_p_Lobatto_porder[i] = e_p[i]->deg;
    /* deg_p_Gauss[i] = e_p_oriented[i]->deg_integ; */

    face_nodes_p_Lobatto[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    face_nodes_p_Gauss[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg_integ);
    
    total_side_nodes_p_Lobatto += face_nodes_p_Lobatto[i];
    total_side_nodes_p_Gauss += face_nodes_p_Gauss[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar_Gauss = 0;
  int total_nodes_mortar_Lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_Gauss[i+j] = util_max_int( e_m[i]->deg_integ,
                                            e_p_oriented[j]->deg_integ);
      deg_mortar_Lobatto[i+j] = util_max_int( e_m[i]->deg,
                                              e_p_oriented[j]->deg );      
      nodes_mortar_Gauss[i+j] = dgmath_get_nodes( (P4EST_DIM) - 1, deg_mortar_Gauss[i+j] );     
      nodes_mortar_Lobatto[i+j] = dgmath_get_nodes( (P4EST_DIM) - 1, deg_mortar_Lobatto[i+j] );     
      total_nodes_mortar_Gauss += nodes_mortar_Gauss[i+j];
      total_nodes_mortar_Lobatto += nodes_mortar_Lobatto[i+j];
      penalty_mortar[i+j] = sipg_kronbichler_flux_penalty_calculate_fcn
                            (
                             e_m[i]->deg,
                             (e_m[i]->volume/e_m[i]->surface_area[f_m]),
                             e_p_oriented[j]->deg,
                             (e_p_oriented[j]->volume/e_p_oriented[j]->surface_area[f_p]),
                             sipg_kronbichler_flux_penalty_prefactor
                            );
      
    }

  int deg_mortar_Gauss_porder [(P4EST_HALF)];
  int nodes_mortar_Gauss_porder [(P4EST_HALF)];
  
  for(int i = 0; i < faces_mortar; i++){
    int inew = i;
    if (faces_mortar == (P4EST_HALF)){
      inew = dgmath_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, i);
    }
    deg_mortar_Gauss_porder[inew] = deg_mortar_Gauss[i];
    nodes_mortar_Gauss_porder[inew] = nodes_mortar_Gauss[i];
  }

  /* slices of scalar/vector fields of (-) onto f_m and (+) onto f_p */
  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m_Lobatto);
  double* u_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p_Lobatto);
  
  /* projections of f_m/f_p on to mortar space */
  double* u_m_on_f_m_mortar = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* u_m_on_f_m_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* u_p_on_f_p_mortar = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* u_p_on_f_p_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* dudr_p_on_f_p_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_Gauss_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar_Gauss [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_Gauss_porder [(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_Gauss [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_Gauss [(P4EST_DIM)];
  double* drst_dxyz_m_on_mortar_Gauss [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_p_on_mortar_Gauss_porder [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_m_on_mortar_Gauss, total_nodes_mortar_Gauss);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_p_on_mortar_Gauss_porder, total_nodes_mortar_Gauss);
  
  for (int i = 0; i < (P4EST_DIM); i++){
    dudr_p_on_f_p_porder[i] = P4EST_ALLOC(double, total_side_nodes_p_Lobatto);
    dudr_p_on_f_p_mortar_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    dudr_p_on_f_p_mortar_Gauss_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    dudx_p_on_f_p_mortar_Gauss_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    dudx_p_on_f_p_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  
    dudr_m_on_f_m[i] = P4EST_ALLOC(double, total_side_nodes_m_Lobatto);
    dudr_m_on_f_m_mortar[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    dudr_m_on_f_m_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    dudx_m_on_f_m_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  }

  double* ones_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  linalg_fill_vec(ones_mortar_Gauss, 1., total_nodes_mortar_Gauss);
  
  double* lifted_proj_VT_w_term1_mortar_Lobatto = P4EST_ALLOC(double, max_volume_nodes_m_Lobatto);
  double* proj_VT_w_term1_mortar_Lobatto = P4EST_ALLOC(double, total_side_nodes_m_Lobatto);
  double* VT_w_term1_mortar_Lobatto = P4EST_ALLOC(double, total_nodes_mortar_Lobatto);
  double* term1_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss); 

  double* DT_lifted_proj_VT_w_term2_mortar_Lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(DT_lifted_proj_VT_w_term2_mortar_Lobatto, max_volume_nodes_m_Lobatto);
  double* lifted_proj_VT_w_term2_mortar_Lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(lifted_proj_VT_w_term2_mortar_Lobatto, max_volume_nodes_m_Lobatto);
  double* proj_VT_w_term2_mortar_Lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(proj_VT_w_term2_mortar_Lobatto, total_side_nodes_m_Lobatto);
  double* VT_w_term2_mortar_Lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(VT_w_term2_mortar_Lobatto, total_nodes_mortar_Lobatto);
  double* term2_mortar_Gauss [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(term2_mortar_Gauss, total_nodes_mortar_Gauss);

  double* lifted_proj_VT_w_term3_mortar_Lobatto = P4EST_ALLOC(double, max_volume_nodes_m_Lobatto);
  double* proj_VT_w_term3_mortar_Lobatto = P4EST_ALLOC(double, total_side_nodes_m_Lobatto);
  double* VT_w_term3_mortar_Lobatto = P4EST_ALLOC(double, total_nodes_mortar_Lobatto);
  double* term3_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);

  
  double* sj_on_f_m_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* n_on_f_m_mortar_Gauss [(P4EST_DIM)];

  double* tmp = P4EST_ALLOC(double, total_side_nodes_p_Gauss);
 
  for (int i = 0; i < (P4EST_DIM); i++) {
    n_on_f_m_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
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
    
    stride += face_nodes_m_Lobatto[i];
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
    
    stride += face_nodes_p_Lobatto[i];
  }


  /* project (-)-side u trace vector onto mortar space */
  dgmath_project_side_onto_mortar_space
    (
     dgmath_jit_dbase,
     u_m_on_f_m,
     faces_m,
     deg_m_Lobatto,
     u_m_on_f_m_mortar,
     faces_mortar,
     deg_mortar_Gauss
    );

  /* project (+)-side u trace vector onto mortar space */
  dgmath_project_side_onto_mortar_space
    (
     dgmath_jit_dbase,
     u_p_on_f_p,
     faces_p,
     deg_p_Lobatto,
     u_p_on_f_p_mortar,
     faces_mortar,
     deg_mortar_Gauss
    );
  
  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &u_m_on_f_m_mortar[stride], deg_mortar_Gauss[f], deg_mortar_Gauss[f], &u_m_on_f_m_mortar_Gauss[stride], (P4EST_DIM)-1);
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &u_p_on_f_p_mortar[stride], deg_mortar_Gauss[f], deg_mortar_Gauss[f], &u_p_on_f_p_mortar_Gauss[stride], (P4EST_DIM)-1);
    stride += nodes_mortar_Gauss[f];
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
      

      stride += face_nodes_m_Lobatto[i];
    }
    
    stride = 0;
    for (int i = 0; i < faces_p; i++){
      for (int d = 0; d < (P4EST_DIM); d++){
        dgmath_apply_slicer
          (
           dgmath_jit_dbase,
           &e_p[i]->dudr_elem[d][0],
           (P4EST_DIM),
           f_p,
           e_p[i]->deg,
           &dudr_p_on_f_p_porder[d][stride]
          );
      }
      stride += dgmath_get_nodes((P4EST_DIM)-1, e_p[i]->deg);
    }


    dgmath_project_side_onto_mortar_space
      (
       dgmath_jit_dbase,
       dudr_p_on_f_p_porder[d],
       faces_p,
       deg_p_Lobatto_porder,
       dudr_p_on_f_p_mortar_porder[d],
       faces_mortar,
       deg_mortar_Gauss_porder
      );
  

    dgmath_project_side_onto_mortar_space
      (
       dgmath_jit_dbase,
       dudr_m_on_f_m[d],
       faces_m,
       deg_m_Lobatto,
       dudr_m_on_f_m_mortar[d],
       faces_mortar,
       deg_mortar_Gauss
      );


    stride = 0;
    for (int f = 0; f < faces_mortar; f++){
      dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &dudr_m_on_f_m_mortar[d][stride], deg_mortar_Gauss[f], deg_mortar_Gauss[f], &dudr_m_on_f_m_mortar_Gauss[d][stride], (P4EST_DIM)-1);
      stride += nodes_mortar_Gauss[f];
    }

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){
      for (int d = 0; d < (P4EST_DIM); d++){
        dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &dudr_p_on_f_p_mortar_porder[d][stride], deg_mortar_Gauss_porder[f], deg_mortar_Gauss_porder[f], &dudr_p_on_f_p_mortar_Gauss_porder[d][stride], (P4EST_DIM)-1);
      }
      stride += nodes_mortar_Gauss_porder[f];
    }
  }

  double* j_div_sj_on_f_m_mortar_Gauss = NULL;
  double* j_div_sj_on_f_p_mortar_Gauss_porder = NULL;
  double* j_div_sj_on_f_p_mortar_Gauss_porder_oriented = NULL;

  if(ip_flux_params->ip_flux_h_calc == H_EQ_J_DIV_SJ){
    j_div_sj_on_f_m_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    j_div_sj_on_f_p_mortar_Gauss_porder =  P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    j_div_sj_on_f_p_mortar_Gauss_porder_oriented =  P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  }
  
  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     faces_m,
     faces_mortar,
     &deg_mortar_Gauss[0],
     f_m,
     drst_dxyz_m_on_mortar_Gauss,
     sj_on_f_m_mortar_Gauss,
     n_on_f_m_mortar_Gauss,
     j_div_sj_on_f_m_mortar_Gauss,
     GAUSS,
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
     &deg_mortar_Gauss_porder[0],
     f_p,
     drst_dxyz_p_on_mortar_Gauss_porder,
     NULL,
     NULL,
     j_div_sj_on_f_p_mortar_Gauss_porder,
     GAUSS,
     geom,
     dgmath_jit_dbase,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  for (int d = 0; d < (P4EST_DIM); d++){
    linalg_fill_vec
      (
       dudx_m_on_f_m_mortar_Gauss[d],
       0.0,
       total_nodes_mortar_Gauss
      );


    linalg_fill_vec
      (
       dudx_p_on_f_p_mortar_Gauss_porder[d],
       0.0,
       total_nodes_mortar_Gauss
      );
  }
  
  for (int j = 0; j < (P4EST_DIM); j++){
    for (int i = 0; i < (P4EST_DIM); i++){   
      for (int k = 0; k < total_nodes_mortar_Gauss; k++){
        dudx_m_on_f_m_mortar_Gauss[j][k] +=
          drst_dxyz_m_on_mortar_Gauss[i][j][k]*dudr_m_on_f_m_mortar_Gauss[i][k];
        dudx_p_on_f_p_mortar_Gauss_porder[j][k] += drst_dxyz_p_on_mortar_Gauss_porder[i][j][k]*dudr_p_on_f_p_mortar_Gauss_porder[i][k];
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
      oriented_face_mortar_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss_porder[b]);
    }


    for (int d = 0; d < (P4EST_DIM); d++){
      dgmath_reorient_face_data
        (
         dgmath_jit_dbase,
         &dudx_p_on_f_p_mortar_Gauss_porder[d][oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_Gauss[face],
         orientation,
         f_m,
         f_p,
         &dudx_p_on_f_p_mortar_Gauss[d][face_mortar_stride]
        );
    }

    if (j_div_sj_on_f_p_mortar_Gauss_porder != NULL){
    dgmath_reorient_face_data
      (
       dgmath_jit_dbase,
       &j_div_sj_on_f_p_mortar_Gauss_porder[oriented_face_mortar_stride],
       (P4EST_DIM)-1,
       deg_mortar_Gauss[face],
       orientation,
       f_m,
       f_p,
       &j_div_sj_on_f_p_mortar_Gauss_porder_oriented[face_mortar_stride]
      );
    }    
    
    
    face_mortar_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss[face]);
  }

  stride = 0;
  double* sigma = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_Gauss[f]; k++){
      int ks = k + stride;
      if (ip_flux_params->ip_flux_h_calc == H_EQ_VOLUME_DIV_AREA){
        sigma[ks] = penalty_mortar[f];
        /* printf("sigma[ks] = %f\n", sigma[ks]); */
      }
      else if (ip_flux_params->ip_flux_h_calc ==H_EQ_J_DIV_SJ){
        double hp = j_div_sj_on_f_p_mortar_Gauss_porder_oriented[ks];
        double hm = j_div_sj_on_f_m_mortar_Gauss[ks];

        sigma[ks] = sipg_kronbichler_flux_penalty_calculate_fcn
                (
                 e_m[f]->deg,
                 hm,
                 e_p_oriented[f]->deg,
                 hp,
                 sipg_kronbichler_flux_penalty_prefactor
                );

      }
      else {
        mpi_abort("Select j_DIV_SJ or VOLUME_DIV_AREA for ip_flux_h_calc");
      }

    }
    stride += nodes_mortar_Gauss[f];
  }



  
  stride = 0;
  int stride_Lobatto = 0;
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_Gauss[f]; k++){
      int ks = k + stride;

      term1_mortar_Gauss[ks] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        term1_mortar_Gauss[ks] += -1.*sj_on_f_m_mortar_Gauss[ks]
                                  *n_on_f_m_mortar_Gauss[d][ks]
                                  *.5*(dudx_p_on_f_p_mortar_Gauss[d][ks] + dudx_m_on_f_m_mortar_Gauss[d][ks]);
      }
        
      for (int l = 0; l < (P4EST_DIM); l++){
        term2_mortar_Gauss[l][ks] = 0.;
        for (int d = 0; d < (P4EST_DIM); d++){
          term2_mortar_Gauss[l][ks] += -.5*sj_on_f_m_mortar_Gauss[ks]
                                       *drst_dxyz_m_on_mortar_Gauss[l][d][ks]
                                       *n_on_f_m_mortar_Gauss[d][ks]
                                       *(u_m_on_f_m_mortar_Gauss[ks] - u_p_on_f_p_mortar_Gauss[ks]);
        }
      }
        
      term3_mortar_Gauss[ks] = sj_on_f_m_mortar_Gauss[ks]*sigma[ks]*(u_m_on_f_m_mortar_Gauss[ks] - u_p_on_f_p_mortar_Gauss[ks]);
    }



    
   dgmath_apply_curvedGaussMass_onGaussNodeVec
      (
       dgmath_jit_dbase,
       &term1_mortar_Gauss[stride],
       deg_mortar_Lobatto[f],
       &ones_mortar_Gauss[stride],
       deg_mortar_Gauss[f],
       (P4EST_DIM)-1,
       &VT_w_term1_mortar_Lobatto[stride_Lobatto]
      );
    
    for (int d = 0; d < (P4EST_DIM); d++){
      dgmath_apply_curvedGaussMass_onGaussNodeVec
        (
         dgmath_jit_dbase,
         &term2_mortar_Gauss[d][stride],
         deg_mortar_Lobatto[f],
         &ones_mortar_Gauss[stride],
         deg_mortar_Gauss[f],
         (P4EST_DIM)-1,
         &VT_w_term2_mortar_Lobatto[d][stride_Lobatto]
        );
    }
      
    dgmath_apply_curvedGaussMass_onGaussNodeVec
      (
       dgmath_jit_dbase,
       &term3_mortar_Gauss[stride],
       deg_mortar_Lobatto[f],
       &ones_mortar_Gauss[stride],
       deg_mortar_Gauss[f],
       (P4EST_DIM)-1,
       &VT_w_term3_mortar_Lobatto[stride_Lobatto]
      );

    stride += nodes_mortar_Gauss[f];
    stride_Lobatto += nodes_mortar_Lobatto[f];
  }
  

  for (int d = 0; d < (P4EST_DIM); d++){
    dgmath_project_mass_mortar_onto_side
      (
       dgmath_jit_dbase,
       VT_w_term2_mortar_Lobatto[d],
       faces_mortar,
       deg_mortar_Lobatto,
       proj_VT_w_term2_mortar_Lobatto[d],
       faces_m,
       deg_m_Lobatto
      );
  }

  dgmath_project_mass_mortar_onto_side
    (
     dgmath_jit_dbase,
     VT_w_term1_mortar_Lobatto,
     faces_mortar,
     deg_mortar_Lobatto,
     proj_VT_w_term1_mortar_Lobatto,
     faces_m,
     deg_m_Lobatto
    );

  dgmath_project_mass_mortar_onto_side
    (
     dgmath_jit_dbase,
     VT_w_term3_mortar_Lobatto,
     faces_mortar,
     deg_mortar_Lobatto,
     proj_VT_w_term3_mortar_Lobatto,
     faces_m,
     deg_m_Lobatto
    );

  /* copy result back to element */
  stride = 0;
  for (int f = 0; f < faces_m; f++){
    if (e_m_is_ghost[f] == 0){

      dgmath_apply_LIFT(
                        dgmath_jit_dbase,
                        &proj_VT_w_term1_mortar_Lobatto[stride],
                        (P4EST_DIM),
                        e_m[f]->deg,
                        f_m,
                        lifted_proj_VT_w_term1_mortar_Lobatto);


      for (int d = 0; d < (P4EST_DIM); d++){
        dgmath_apply_LIFT(
                          dgmath_jit_dbase,
                          &proj_VT_w_term2_mortar_Lobatto[d][stride],
                          (P4EST_DIM),
                          e_m[f]->deg,
                          f_m,
                          lifted_proj_VT_w_term2_mortar_Lobatto[d]);

        dgmath_apply_Dij_transpose(dgmath_jit_dbase,
                                   lifted_proj_VT_w_term2_mortar_Lobatto[d],
                                   (P4EST_DIM),
                                   e_m[f]->deg,
                                   d,
                                   DT_lifted_proj_VT_w_term2_mortar_Lobatto[d]
                                  );

      }
        

      dgmath_apply_LIFT(
                        dgmath_jit_dbase,
                        &proj_VT_w_term3_mortar_Lobatto[stride],
                        (P4EST_DIM),
                        e_m[f]->deg,
                        f_m,
                        lifted_proj_VT_w_term3_mortar_Lobatto);



        
      int volume_nodes_m = dgmath_get_nodes((P4EST_DIM), e_m[f]->deg);
      for (int i = 0; i < volume_nodes_m; i++){
        for (int d = 0; d < (P4EST_DIM); d++){
          e_m[f]->Au_elem[i] += DT_lifted_proj_VT_w_term2_mortar_Lobatto[d][i];
        }
        e_m[f]->Au_elem[i] += lifted_proj_VT_w_term3_mortar_Lobatto[i];
        e_m[f]->Au_elem[i] += lifted_proj_VT_w_term1_mortar_Lobatto[i];
      }
    }
    stride += face_nodes_m_Lobatto[f];
  }   
 
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_p_on_f_p);

  D4EST_FREE_DBYD_MAT(drst_dxyz_m_on_mortar_Gauss);
  D4EST_FREE_DBYD_MAT(drst_dxyz_p_on_mortar_Gauss_porder);

  P4EST_FREE(sigma);

  P4EST_FREE(u_m_on_f_m_mortar);
  P4EST_FREE(u_m_on_f_m_mortar_Gauss);
  P4EST_FREE(u_p_on_f_p_mortar);
  P4EST_FREE(u_p_on_f_p_mortar_Gauss);
  for (int i = 0; i < (P4EST_DIM); i++){   
    P4EST_FREE(dudr_p_on_f_p_porder[i]);
    P4EST_FREE(dudr_p_on_f_p_mortar_porder[i]);
    P4EST_FREE(dudr_p_on_f_p_mortar_Gauss_porder[i]);
    P4EST_FREE(dudx_p_on_f_p_mortar_Gauss_porder[i]);
    P4EST_FREE(dudx_p_on_f_p_mortar_Gauss[i]);
    P4EST_FREE(dudr_m_on_f_m[i]);
    P4EST_FREE(dudr_m_on_f_m_mortar[i]);
    P4EST_FREE(dudr_m_on_f_m_mortar_Gauss[i]);
    P4EST_FREE(dudx_m_on_f_m_mortar_Gauss[i]);
  }
  P4EST_FREE(sj_on_f_m_mortar_Gauss);
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(n_on_f_m_mortar_Gauss[i]);
  }

  if (j_div_sj_on_f_m_mortar_Gauss != NULL){
  P4EST_FREE(j_div_sj_on_f_m_mortar_Gauss);
  P4EST_FREE(j_div_sj_on_f_p_mortar_Gauss_porder);
  P4EST_FREE(j_div_sj_on_f_p_mortar_Gauss_porder_oriented);
  }
  P4EST_FREE(ones_mortar_Gauss);
  P4EST_FREE(lifted_proj_VT_w_term1_mortar_Lobatto);
  P4EST_FREE(proj_VT_w_term1_mortar_Lobatto);
  P4EST_FREE(VT_w_term1_mortar_Lobatto);
  P4EST_FREE(term1_mortar_Gauss);

  D4EST_FREE_DIM_VEC(DT_lifted_proj_VT_w_term2_mortar_Lobatto);
  D4EST_FREE_DIM_VEC(lifted_proj_VT_w_term2_mortar_Lobatto);
  D4EST_FREE_DIM_VEC(proj_VT_w_term2_mortar_Lobatto);
  D4EST_FREE_DIM_VEC(VT_w_term2_mortar_Lobatto);
  D4EST_FREE_DIM_VEC(term2_mortar_Gauss);

  P4EST_FREE(lifted_proj_VT_w_term3_mortar_Lobatto);
  P4EST_FREE(proj_VT_w_term3_mortar_Lobatto);
  P4EST_FREE(VT_w_term3_mortar_Lobatto);
  P4EST_FREE(term3_mortar_Gauss);
  
  P4EST_FREE(tmp);  
}

curved_flux_fcn_ptrs_t
curved_Gauss_primal_sipg_kronbichler_flux_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 ip_flux_params_t* curved_Gauss_sipg_params
)
{  
  curved_flux_fcn_ptrs_t curved_Gauss_primal_sipg_kronbichler_flux_fcns;
  curved_Gauss_primal_sipg_kronbichler_flux_fcns.flux_interface_fcn = curved_Gauss_primal_sipg_kronbichler_flux_interface;
  curved_Gauss_primal_sipg_kronbichler_flux_fcns.flux_boundary_fcn = curved_Gauss_primal_sipg_kronbichler_flux_dirichlet;
  curved_Gauss_primal_sipg_kronbichler_flux_fcns.bndry_fcn = bndry_fcn;
  curved_Gauss_primal_sipg_kronbichler_flux_fcns.params = (void*)curved_Gauss_sipg_params;
  
  return curved_Gauss_primal_sipg_kronbichler_flux_fcns;
}
