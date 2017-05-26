#include "../Utilities/util.h"
#include "../dGMath/d4est_operators.h"
#include "../ElementData/curved_element_data.h"
#include "../LinearAlgebra/linalg.h"
#include "../Flux/curved_Gauss_sipg_flux_vector_fcns.h"

/* #define DEALIASING */

static void
curved_Gauss_sipg_flux_vector_dirichlet
(
 curved_element_data_t* e_m,
 int f_m,
 grid_fcn_t bndry_fcn,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* geom,
 void* params
)
{

  /* int deg_face_Gauss_quad = e_m->deg_quad; */
  
  ip_flux_params_t* ip_params = (ip_flux_params_t*) params;
  double sipg_flux_penalty_prefactor = ip_params->ip_flux_penalty_prefactor;
  penalty_calc_t sipg_flux_penalty_calculate_fcn = ip_params->ip_flux_penalty_calculate_fcn;
  double h = (e_m->volume/e_m->surface_area[f_m]);  
  double sigma = sipg_flux_penalty_calculate_fcn
                   (
                    e_m->deg,
                    h,
                    e_m->deg,
                    h,
                    sipg_flux_penalty_prefactor
                   );

  
  grid_fcn_t u_at_bndry = bndry_fcn;
  int face_nodes_m_Lobatto = d4est_operators_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_Gauss = d4est_operators_get_nodes((P4EST_DIM) - 1, e_m->deg_quad);

  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* u_m_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* q_m_on_f_m [(P4EST_DIM)];
  double* q_m_on_f_m_Gauss [(P4EST_DIM)];

  double* dudr_m_on_f_m [(P4EST_DIM)];
  double* dudr_m_on_f_m_Gauss [(P4EST_DIM)];
  double* dudx_m_on_f_m_Gauss [(P4EST_DIM)];
  
  double* xyz_on_f_m [(P4EST_DIM)];
  double* M_qstar_min_q_n = P4EST_ALLOC_ZERO(double, face_nodes_m_Lobatto);
  double* u_on_f_m_min_u_at_bndry_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* u_on_f_m_min_u_at_bndry_Lobatto = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* sj_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* qstar_min_q_Gauss [(P4EST_DIM)];
  double* n_on_f_m_Gauss [(P4EST_DIM)];
  double* sj_n_on_f_m_Gauss [(P4EST_DIM)];


  double* drst_dxyz_Gauss [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_Gauss, face_nodes_m_Gauss);
  
  for (int d = 0; d < (P4EST_DIM); d++) {
    xyz_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m_Lobatto);
    q_m_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m_Lobatto);
    dudr_m_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m_Lobatto);
    dudr_m_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
    dudx_m_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
    q_m_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
    n_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
    sj_n_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
    qstar_min_q_Gauss[d] = P4EST_ALLOC_ZERO(double, face_nodes_m_Gauss);
  }

  for (int d = 0; d < (P4EST_DIM); d++){

    d4est_operators_apply_slicer(d4est_ops,
                        e_m->xyz[d],
                        (P4EST_DIM),
                        f_m,
                        e_m->deg,
                        xyz_on_f_m[d]);
    
    d4est_operators_apply_slicer
      (
       d4est_ops,
       e_m->q_elem[d],
       (P4EST_DIM),
       f_m,
       e_m->deg,
       q_m_on_f_m[d]
      );
    
   d4est_operators_interp_GLL_to_GL
     (
      d4est_ops,
      q_m_on_f_m[d],
      e_m->deg,
      e_m->deg_quad,
      q_m_on_f_m_Gauss[d],
      (P4EST_DIM)-1
     );

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


 curved_data_compute_drst_dxyz_Gauss_on_mortar_using_volume_data
   (
     &e_m,
     1,
     1,
     &e_m->deg_quad,
     f_m,
     drst_dxyz_Gauss,
     sj_on_f_m_Gauss,
     n_on_f_m_Gauss,
     geom->p4est_geom,
     d4est_ops
    );

  
  /* curved_element_data_compute_physical_derivatives_on_face_Gauss_nodes */
  /*   ( */
  /*    dudr_m_on_f_m_Gauss, /\* should be of mortar length, but not individually rotated *\/ */
  /*    &e_m, */
  /*    1, */
  /*    1, */
  /*    &e_m->deg_quad, */
  /*    f_m, */
  /*    dudx_m_on_f_m_Gauss, */
  /*    geom->p4est_geom, */
  /*    d4est_ops */
  /*   ); */

  /* DEBUG_PRINT_ARR_DBL(dudr_m_on_f_m_Gauss[0], d4est_operators_get_nodes((P4EST_DIM-1), e_m->deg_quad)); */
  /* DEBUG_PRINT_ARR_DBL(dudx_m_on_f_m_Gauss[0], d4est_operators_get_nodes((P4EST_DIM-1), e_m->deg_quad)); */

  
  d4est_operators_apply_slicer(d4est_ops, e_m->u_storage, (P4EST_DIM), f_m, e_m->deg, u_m_on_f_m);
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

  for (int i = 0; i < face_nodes_m_Lobatto; i++){
    u_on_f_m_min_u_at_bndry_Lobatto[i] = u_m_on_f_m[i]
                                         - u_at_bndry
                                         (
                                          xyz_on_f_m[0][i],
                                          xyz_on_f_m[1][i]
#if (P4EST_DIM)==3
                                          ,
                                          xyz_on_f_m[2][i]
#endif
                                         );
  }

  d4est_operators_interp_GLL_to_GL
    (
     d4est_ops,
     u_on_f_m_min_u_at_bndry_Lobatto,
     e_m->deg,
     e_m->deg_quad,
     u_on_f_m_min_u_at_bndry_Gauss,
     (P4EST_DIM)-1
    );
  
  
  for(int i = 0; i < face_nodes_m_Gauss; i++){
    for (int d = 0; d < (P4EST_DIM); d++){
      qstar_min_q_Gauss[d][i] = dudx_m_on_f_m_Gauss[d][i]
                       - sigma*n_on_f_m_Gauss[d][i]*(u_on_f_m_min_u_at_bndry_Gauss[i])
                       - q_m_on_f_m_Gauss[d][i];
      sj_n_on_f_m_Gauss[d][i] = sj_on_f_m_Gauss[i]*n_on_f_m_Gauss[d][i];
    }
  }

  for (int i = 0; i < face_nodes_m_Lobatto; i++){
    e_m->M_qstar_min_q_dot_n[f_m*face_nodes_m_Lobatto + i] = 0.;
  }
  
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_curvedGaussMass_onGaussNodeVec(d4est_ops,
                                                qstar_min_q_Gauss[d],
                                                e_m->deg,
                                                sj_n_on_f_m_Gauss[d],
                                                e_m->deg_quad,
                                                (P4EST_DIM)-1,
                                                M_qstar_min_q_n);
    
    for (int i = 0; i < face_nodes_m_Lobatto; i++){
      e_m->M_qstar_min_q_dot_n[f_m*face_nodes_m_Lobatto + i] += M_qstar_min_q_n[i];

      /* printf("e_m->M_qstar_min_q_dot_n[f_m*face_nodes_m_Lobatto + i] = %.25f\n",e_m->M_qstar_min_q_dot_n[f_m*face_nodes_m_Lobatto + i]); */
    }
  }

  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_m_on_f_m_Gauss);
  for (int i = 0; i < (P4EST_DIM); i++) {
    P4EST_FREE(q_m_on_f_m[i]);
    P4EST_FREE(q_m_on_f_m_Gauss[i]);
    P4EST_FREE(dudr_m_on_f_m[i]);
    P4EST_FREE(dudr_m_on_f_m_Gauss[i]);
    P4EST_FREE(dudx_m_on_f_m_Gauss[i]);
    P4EST_FREE(xyz_on_f_m[i]);
    P4EST_FREE(qstar_min_q_Gauss[i]);
  }

  P4EST_FREE(sj_on_f_m_Gauss);
  P4EST_FREE(u_on_f_m_min_u_at_bndry_Gauss);
  P4EST_FREE(u_on_f_m_min_u_at_bndry_Lobatto);
  P4EST_FREE(M_qstar_min_q_n);
  
  for (int d = 0; d < (P4EST_DIM); d++){
    P4EST_FREE(n_on_f_m_Gauss[d]);
    P4EST_FREE(sj_n_on_f_m_Gauss[d]);
  }
  D4EST_FREE_DBYD_MAT(drst_dxyz_Gauss);
  
}

static void
curved_Gauss_sipg_flux_vector_interface
(
 curved_element_data_t** e_m,
 int faces_m,
 int f_m,
 curved_element_data_t** e_p,
 int faces_p,
 int f_p,
 int* e_m_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* geom,
 void* params
)
{
  ip_flux_params_t* ip_params = (ip_flux_params_t*) params;
  double sipg_flux_penalty_prefactor = ip_params->ip_flux_penalty_prefactor;
  penalty_calc_t sipg_flux_penalty_calculate_fcn = ip_params->ip_flux_penalty_calculate_fcn;
  
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
  int total_side_nodes_m_Lobatto = 0;
  int total_side_nodes_m_Gauss = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_Lobatto[i] = e_m[i]->deg;
    /* deg_m_Gauss[i] = e_m[i]->deg_quad; */
    
    face_nodes_m_Lobatto[i] = d4est_operators_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg);
    face_nodes_m_Gauss[i] = d4est_operators_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg_quad);
    
    total_side_nodes_m_Lobatto += face_nodes_m_Lobatto[i];
    total_side_nodes_m_Gauss += face_nodes_m_Gauss[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_Lobatto = 0;
  int total_side_nodes_p_Gauss = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p_Lobatto[i] = e_p_oriented[i]->deg;
    deg_p_Lobatto_porder[i] = e_p[i]->deg;
    /* deg_p_Gauss[i] = e_p_oriented[i]->deg_quad; */

    face_nodes_p_Lobatto[i] = d4est_operators_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    face_nodes_p_Gauss[i] = d4est_operators_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg_quad);
    
    total_side_nodes_p_Lobatto += face_nodes_p_Lobatto[i];
    total_side_nodes_p_Gauss += face_nodes_p_Gauss[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar_Gauss = 0;
  int total_nodes_mortar_Lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_Gauss[i+j] = util_max_int( e_m[i]->deg_quad,
                                            e_p_oriented[j]->deg_quad);
      deg_mortar_Lobatto[i+j] = util_max_int( e_m[i]->deg,
                                            e_p_oriented[j]->deg );      
      nodes_mortar_Gauss[i+j] = d4est_operators_get_nodes( (P4EST_DIM) - 1, deg_mortar_Gauss[i+j] );     
      nodes_mortar_Lobatto[i+j] = d4est_operators_get_nodes( (P4EST_DIM) - 1, deg_mortar_Lobatto[i+j] );     
      total_nodes_mortar_Gauss += nodes_mortar_Gauss[i+j];
      total_nodes_mortar_Lobatto += nodes_mortar_Lobatto[i+j];
      penalty_mortar[i+j] = sipg_flux_penalty_calculate_fcn
                            (
                             e_m[i]->deg,
                             (e_m[i]->volume/e_m[i]->surface_area[f_m]),
                             e_p[j]->deg,
                             (e_p_oriented[j]->volume/e_p_oriented[j]->surface_area[f_p]),
                             sipg_flux_penalty_prefactor
                            );
      
    }

  int deg_mortar_Gauss_porder [(P4EST_HALF)];
  int nodes_mortar_Gauss_porder [(P4EST_HALF)];
  
  for(int i = 0; i < faces_mortar; i++){
    int inew = i;
    if (faces_mortar == (P4EST_HALF)){
      inew = d4est_operators_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, i);
    }
    deg_mortar_Gauss_porder[inew] = deg_mortar_Gauss[i];
    nodes_mortar_Gauss_porder[inew] = nodes_mortar_Gauss[i];
  }

  /* slices of scalar/vector fields of (-) onto f_m and (+) onto f_p */
  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m_Lobatto);
  double* u_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p_Lobatto);
  double* q_p_on_f_p [(P4EST_DIM)];
  double* q_m_on_f_m [(P4EST_DIM)];

  for (int i = 0; i < (P4EST_DIM); i++) {
    q_m_on_f_m[i] = P4EST_ALLOC(double, total_side_nodes_m_Lobatto);
    q_p_on_f_p[i] = P4EST_ALLOC(double, total_side_nodes_p_Lobatto);
  }
  
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
  double* q_m_on_f_m_mortar [(P4EST_DIM)];
  double* q_m_on_f_m_mortar_Gauss [(P4EST_DIM)];

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


    q_m_on_f_m_mortar[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    q_m_on_f_m_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  }
  
  double* sj_on_f_m_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  /* double* sj_on_f_m_mortar_Gauss_m = P4EST_ALLOC(double, total_nodes_mortar_Gauss); */
  double* n_on_f_m_mortar_Gauss [(P4EST_DIM)];
  /* double* n_on_f_m_mortar_Gauss_m [(P4EST_DIM)]; */
  double* sj_n_on_f_m_mortar_Gauss [(P4EST_DIM)];

  double* tmp = P4EST_ALLOC(double, total_side_nodes_p_Gauss);
  double* qstar_min_q_mortar_Gauss [(P4EST_DIM)];
  double* proj_M_qstar_min_q_dot_n_mortar = P4EST_ALLOC_ZERO(double, total_side_nodes_m_Lobatto);
  double* M_qstar_min_q_dot_n_mortar = P4EST_ALLOC_ZERO(double, total_nodes_mortar_Gauss);
  double* M_qstar_min_q_n_mortar = P4EST_ALLOC_ZERO(double, total_nodes_mortar_Gauss);
  
  for (int i = 0; i < (P4EST_DIM); i++) {
    n_on_f_m_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    /* n_on_f_m_mortar_Gauss_m[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss); */
    sj_n_on_f_m_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    qstar_min_q_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  }  
  
  stride = 0;
  for (int i = 0; i < faces_m; i++){
    d4est_operators_apply_slicer
      (
       d4est_ops,
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
    d4est_operators_apply_slicer
      (
       d4est_ops,
       &(e_p_oriented[i]->u_storage[0]),
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
    
    stride += face_nodes_p_Lobatto[i];
  }


  /* project (-)-side u trace vector onto mortar space */
  d4est_operators_project_side_onto_mortar_space
    (
     d4est_ops,
     u_m_on_f_m,
     faces_m,
     deg_m_Lobatto,
     u_m_on_f_m_mortar,
     faces_mortar,
     deg_mortar_Gauss
    );

  /* project (+)-side u trace vector onto mortar space */
  d4est_operators_project_side_onto_mortar_space
    (
     d4est_ops,
     u_p_on_f_p,
     faces_p,
     deg_p_Lobatto,
     u_p_on_f_p_mortar,
     faces_mortar,
     deg_mortar_Gauss
    );
  
  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    d4est_operators_interp_GLL_to_GL(d4est_ops, &u_m_on_f_m_mortar[stride], deg_mortar_Gauss[f], deg_mortar_Gauss[f], &u_m_on_f_m_mortar_Gauss[stride], (P4EST_DIM)-1);
    d4est_operators_interp_GLL_to_GL(d4est_ops, &u_p_on_f_p_mortar[stride], deg_mortar_Gauss[f], deg_mortar_Gauss[f], &u_p_on_f_p_mortar_Gauss[stride], (P4EST_DIM)-1);
    stride += nodes_mortar_Gauss[f];
   }

  
  /* For each component of the vector */
  for (int d = 0; d < (P4EST_DIM); d++){

    /* project (-)-u-derivative on the (-)-side faces and project q onto the (-)-side faces */
    stride = 0;
    for (int i = 0; i < faces_m; i++){    
      d4est_operators_apply_slicer
        (
         d4est_ops,
         e_m[i]->q_elem[d],
         (P4EST_DIM),
         f_m,
         e_m[i]->deg,
         &q_m_on_f_m[d][stride]
        );

     d4est_operators_apply_slicer
        (
         d4est_ops,
         &e_m[i]->dudr_elem[d][0],
         (P4EST_DIM),
         f_m,
         e_m[i]->deg,
         &dudr_m_on_f_m[d][stride]
        );
      

      stride += face_nodes_m_Lobatto[i];
    }
    
    /* compute the (+)-u-derivative 
   * and project on the (+)-side faces and project q onto the (+)-side faces */
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
      stride += d4est_operators_get_nodes((P4EST_DIM)-1, e_p[i]->deg);
    }


    d4est_operators_project_side_onto_mortar_space
      (
       d4est_ops,
       dudr_p_on_f_p_porder[d],
       faces_p,
       deg_p_Lobatto_porder,
       dudr_p_on_f_p_mortar_porder[d],
       faces_mortar,
       deg_mortar_Gauss_porder
      );
  

    d4est_operators_project_side_onto_mortar_space
      (
       d4est_ops,
       dudr_m_on_f_m[d],
       faces_m,
       deg_m_Lobatto,
       dudr_m_on_f_m_mortar[d],
       faces_mortar,
       deg_mortar_Gauss
      );


    /* project q from the (-) side onto the mortar space */
    d4est_operators_project_side_onto_mortar_space
      (
       d4est_ops,
       q_m_on_f_m[d],
       faces_m,
       deg_m_Lobatto,
       q_m_on_f_m_mortar[d],
       faces_mortar,
       deg_mortar_Gauss
      );

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){
      d4est_operators_interp_GLL_to_GL(d4est_ops, &q_m_on_f_m_mortar[d][stride], deg_mortar_Gauss[f], deg_mortar_Gauss[f], &q_m_on_f_m_mortar_Gauss[d][stride], (P4EST_DIM)-1);
      d4est_operators_interp_GLL_to_GL(d4est_ops, &dudr_m_on_f_m_mortar[d][stride], deg_mortar_Gauss[f], deg_mortar_Gauss[f], &dudr_m_on_f_m_mortar_Gauss[d][stride], (P4EST_DIM)-1);
      stride += nodes_mortar_Gauss[f];
    }

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){
      for (int d = 0; d < (P4EST_DIM); d++){
        d4est_operators_interp_GLL_to_GL(d4est_ops, &dudr_p_on_f_p_mortar_porder[d][stride], deg_mortar_Gauss_porder[f], deg_mortar_Gauss_porder[f], &dudr_p_on_f_p_mortar_Gauss_porder[d][stride], (P4EST_DIM)-1);
      }
      stride += nodes_mortar_Gauss_porder[f];
    }
  }


 curved_data_compute_drst_dxyz_Gauss_on_mortar_using_volume_data
    (
     e_m,
     faces_m,
     faces_mortar,
     &deg_mortar_Gauss[0],
     f_m,
     drst_dxyz_m_on_mortar_Gauss,
     sj_on_f_m_mortar_Gauss,
     n_on_f_m_mortar_Gauss,
     geom->p4est_geom,
     d4est_ops
    );

  curved_data_compute_drst_dxyz_Gauss_on_mortar_using_volume_data
    (
     e_p,
     faces_p,
     faces_mortar,
     &deg_mortar_Gauss_porder[0], 
     f_p,
     drst_dxyz_p_on_mortar_Gauss_porder,
     NULL,
     (double* [(P4EST_DIM)]){NULL, NULL
#if (P4EST_DIM)==3
         , NULL
#endif
         },
     geom->p4est_geom,
     d4est_ops
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
      face_p = d4est_operators_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, face);

    int oriented_face_mortar_stride = 0;
    for (int b = 0; b < face_p; b++){
      oriented_face_mortar_stride += d4est_operators_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss_porder[b]);
    }


    for (int d = 0; d < (P4EST_DIM); d++){
      d4est_operators_reorient_face_data
        (
         d4est_ops,
         &dudx_p_on_f_p_mortar_Gauss_porder[d][oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_Gauss[face],
         orientation,
         f_m,
         f_p,
         &dudx_p_on_f_p_mortar_Gauss[d][face_mortar_stride]
        );
    }
    
    face_mortar_stride += d4est_operators_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss[face]);
  }

  
  
  
  /* double* tmpxyz [(P4EST_DIM)]; */
  /* D4EST_ALLOC_DIM_VEC(tmpxyz,total_nodes_mortar_Gauss); */
  /* curved_element_data_compute_mortar_normal_and_sj_using_face_data */
  /*   ( */
  /*    e_m, */
  /*    faces_m, */
  /*    faces_mortar, */
  /*    &deg_mortar_Gauss[0], */
  /*    f_m, */
  /*    geom->dxdr_method, */
  /*    1, */
  /*    n_on_f_m_mortar_Gauss, */
  /*    sj_on_f_m_mortar_Gauss, */
  /*    geom, */
  /*    d4est_ops, */
  /*    tmpxyz */
  /*   ); */
  /* D4EST_FREE_DIM_VEC(tmpxyz); */
  
    /* calculate symmetric interior penalty flux */
  /* if ((e_m[0]->id == 0 && e_p[0]->id == 1) || (e_m[0]->id == 1 && e_p[0]->id == 0)){ */
  /*   printf("e_m[0]->id = %d, e_p[0]->id = %d, f_m = %d, f_p = %d\n", e_m[0]->id, e_p[0]->id, f_m, f_p); */
  /* } */

    stride = 0;
    int stride_Lobatto = 0;
    for (int f = 0; f < faces_mortar; f++){
      for (int k = 0; k < nodes_mortar_Gauss[f]; k++){
        int ks = k + stride;
        double sigma = penalty_mortar[f];
        
        for (int d = 0; d < (P4EST_DIM); d++){
          qstar_min_q_mortar_Gauss[d][ks] = .5*(dudx_p_on_f_p_mortar_Gauss[d][ks] + dudx_m_on_f_m_mortar_Gauss[d][ks]);
          qstar_min_q_mortar_Gauss[d][ks] -= -sigma*n_on_f_m_mortar_Gauss[d][ks]*u_p_on_f_p_mortar_Gauss[ks];
          qstar_min_q_mortar_Gauss[d][ks] -= sigma*n_on_f_m_mortar_Gauss[d][ks]*u_m_on_f_m_mortar_Gauss[ks];

          /* if ((e_m[0]->id == 0 && e_p[0]->id == 1) || (e_m[0]->id == 1 && e_p[0]->id == 0)){ */
          /*   printf("qstar_min_q_mortar_Gauss[%d][%d] = %.25f\n",d,ks, qstar_min_q_mortar_Gauss[d][ks]); */
          /* } */
          
          qstar_min_q_mortar_Gauss[d][ks] -= q_m_on_f_m_mortar_Gauss[d][ks];
          sj_n_on_f_m_mortar_Gauss[d][ks] = sj_on_f_m_mortar_Gauss[ks]*n_on_f_m_mortar_Gauss[d][ks];

          /* printf("dudx_m, dudx_p, n, q, sj = %.25f, %.25f, %.25f, %.25f, %.25f\n", dudx_m_on_f_m_mortar_Gauss[d][ks], dudx_p_on_f_p_mortar_Gauss_reoriented[d][ks], n_on_f_m_mortar_Gauss[d][ks], q_m_on_f_m_mortar_Gauss[d][ks], sj_on_f_m_mortar_Gauss[ks]); */
        }        
      }
  
      for (int d = 0; d < (P4EST_DIM); d++){
        d4est_operators_apply_curvedGaussMass_onGaussNodeVec(
                                     d4est_ops,
                                     &qstar_min_q_mortar_Gauss[d][stride],
                                     deg_mortar_Lobatto[f],
                                     &sj_n_on_f_m_mortar_Gauss[d][stride],
                                     deg_mortar_Gauss[f],
                                     (P4EST_DIM)-1,
                                     &M_qstar_min_q_n_mortar[stride_Lobatto]
                                    );
    
        for (int i = 0; i < nodes_mortar_Lobatto[f]; i++){
          M_qstar_min_q_dot_n_mortar[stride_Lobatto + i] += M_qstar_min_q_n_mortar[stride_Lobatto + i];
        }
      }
      stride += nodes_mortar_Gauss[f];
      stride_Lobatto += nodes_mortar_Lobatto[f];
      /* stride_Lobatto += d4est_operators_get_nodes((P4EST_DIM)-1, e_m[f]->deg); */
    }

    /* if (faces_mortar == P4EST_HALF) */
      /* mpi_abort("Only testing p-nonconforming atm"); */
    
    d4est_operators_project_mass_mortar_onto_side
      (
       d4est_ops,
       M_qstar_min_q_dot_n_mortar,
       faces_mortar,
       deg_mortar_Lobatto,
       proj_M_qstar_min_q_dot_n_mortar,
       faces_m,
       deg_m_Lobatto
      );
    
  /* copy result back to element */
  stride = 0;
  for (int i = 0; i < faces_m; i++){
    if(e_m_is_ghost[i] == 0)
      linalg_copy_1st_to_2nd
        (
         &proj_M_qstar_min_q_dot_n_mortar[stride],
         &(e_m[i]->M_qstar_min_q_dot_n[f_m*face_nodes_m_Lobatto[i]]),
         face_nodes_m_Lobatto[i]
        );
    stride += face_nodes_m_Lobatto[i];
  }

  /* stride = 0; */
  /* for (int i = 0; i < faces_m; i++){ */
  /*   if(e_m_is_ghost[i] == 0) */
  /*     linalg_copy_1st_to_2nd */
  /*       ( */
  /*        &M_qstar_min_q_dot_n_mortar[stride], */
  /*        &(e_m[i]->M_qstar_min_q_dot_n[f_m*face_nodes_m_Lobatto[i]]), */
  /*        face_nodes_m_Lobatto[i] */
  /*       ); */
  /*   stride += face_nodes_m_Lobatto[i]; */
  /* } */
    
 
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_p_on_f_p);
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(q_p_on_f_p[i]);
    P4EST_FREE(q_m_on_f_m[i]);
  }


  D4EST_FREE_DBYD_MAT(drst_dxyz_m_on_mortar_Gauss);
  D4EST_FREE_DBYD_MAT(drst_dxyz_p_on_mortar_Gauss_porder);
  
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
    P4EST_FREE(q_m_on_f_m_mortar[i]);
    P4EST_FREE(q_m_on_f_m_mortar_Gauss[i]);
  }

  P4EST_FREE(sj_on_f_m_mortar_Gauss);
  /* P4EST_FREE(sj_on_f_m_mortar_Gauss_m); */
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(n_on_f_m_mortar_Gauss[i]);
    /* P4EST_FREE(n_on_f_m_mortar_Gauss_m[i]); */
    P4EST_FREE(sj_n_on_f_m_mortar_Gauss[i]);
    P4EST_FREE(qstar_min_q_mortar_Gauss[i]);
  }    
  P4EST_FREE(M_qstar_min_q_n_mortar);
  P4EST_FREE(M_qstar_min_q_dot_n_mortar);
  P4EST_FREE(proj_M_qstar_min_q_dot_n_mortar);
  P4EST_FREE(tmp);  
  
}

curved_flux_fcn_ptrs_t
curved_Gauss_sipg_flux_vector_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 ip_flux_params_t* curved_Gauss_sipg_params
)
{  
  curved_flux_fcn_ptrs_t curved_Gauss_sipg_flux_vector_fcns;
  curved_Gauss_sipg_flux_vector_fcns.flux_interface_fcn = curved_Gauss_sipg_flux_vector_interface;
  curved_Gauss_sipg_flux_vector_fcns.flux_boundary_fcn = curved_Gauss_sipg_flux_vector_dirichlet;
  curved_Gauss_sipg_flux_vector_fcns.bndry_fcn = bndry_fcn;
  curved_Gauss_sipg_flux_vector_fcns.params = (void*)curved_Gauss_sipg_params;
  
  return curved_Gauss_sipg_flux_vector_fcns;
}
