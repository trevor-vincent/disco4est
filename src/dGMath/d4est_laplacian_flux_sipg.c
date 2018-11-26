 #include <pXest.h>
#include <d4est_util.h>
#include <ini.h>
#include <d4est_operators.h>
#include <d4est_element_data.h>
#include <d4est_linalg.h>
#include <d4est_laplacian_flux_sipg.h>
#include <d4est_elliptic_data.h>
#include <d4est_mortars.h>
#include <d4est_quadrature_compactified.h>
#include <d4est_quadrature_lobatto.h>
#include <d4est_laplacian_flux.h>

void
d4est_laplacian_flux_sipg_dirichlet_aux
(
 p4est_t* p4est,
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_laplacian_flux_boundary_data_t* boundary_data,
 void* boundary_condition_fcn_data,
 void* flux_parameter_data,
 double * restrict  term1_quad,
 double * restrict  VT_w_term1_lobatto,
 double * restrict  lifted_VT_w_term1_lobatto,
 double *  term2_quad [P4EST_DIM],
 double *   VT_w_term2_lobatto [P4EST_DIM],
 double *   lifted_VT_w_term2_lobatto [P4EST_DIM],
 double *   DT_lifted_VT_w_term2_lobatto [P4EST_DIM],
 double * restrict  term3_quad,
 double * restrict  VT_w_term3_lobatto,
 double * restrict  lifted_VT_w_term3_lobatto,
 double * restrict  sigma,
 double * restrict  u_at_bndry_lobatto,
 double * restrict  u_at_bndry_lobatto_to_quad
){
  d4est_laplacian_dirichlet_bc_t* bc_data = boundary_condition_fcn_data;
  
  d4est_quadrature_mortar_t* face_object = boundary_data->face_object;
  int deg_mortar_quad = boundary_data->deg_mortar_quad;
  int face_nodes_m_lobatto = boundary_data->face_nodes_m_lobatto;
  int face_nodes_m_quad = boundary_data->face_nodes_m_quad;
  double* u_m_on_f_m_quad = boundary_data->u_m_on_f_m_quad;
  double* sj_on_f_m_quad = boundary_data->sj_on_f_m_quad;
  double* xyz_on_f_m_lobatto [(P4EST_DIM)]; 
  double* drst_dxyz_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* dudx_m_on_f_m_quad [(P4EST_DIM)];  
  double* n_on_f_m_quad [(P4EST_DIM)];
  D4EST_COPY_DBYD_MAT(boundary_data->drst_dxyz_quad, drst_dxyz_quad);
  D4EST_COPY_DIM_VEC(boundary_data->dudx_m_on_f_m_quad, dudx_m_on_f_m_quad);
  D4EST_COPY_DIM_VEC(boundary_data->n_on_f_m_quad, n_on_f_m_quad);
  D4EST_COPY_DIM_VEC(boundary_data->xyz_on_f_m_lobatto, xyz_on_f_m_lobatto);
  
  d4est_laplacian_flux_sipg_params_t* ip_flux_params = (d4est_laplacian_flux_sipg_params_t*) flux_parameter_data;

  d4est_laplacian_dirichlet_bc_t* bc_params = boundary_condition_fcn_data;
  
  double* ones_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  d4est_util_fill_array(ones_quad, 1., face_nodes_m_quad);

  double* h_quad = boundary_data->h_quad;

  d4est_mesh_size_parameters_t size_params = (ip_flux_params->size_params == NULL) ? d4est_mesh_get_size_parameters(d4est_factors) : *ip_flux_params->size_params;
  
  for (int i = 0; i < face_nodes_m_quad; i++){
    sigma[i] = ip_flux_params->sipg_penalty_fcn
               (
                e_m->deg,
                h_quad[i],
                e_m->deg,
                h_quad[i],
                ip_flux_params->sipg_penalty_prefactor
               );

  }
  
  for (int i = 0; i < face_nodes_m_lobatto; i++){
    u_at_bndry_lobatto[i] = bc_data->dirichlet_fcn
                            (
                             xyz_on_f_m_lobatto[0][i],
                             xyz_on_f_m_lobatto[1][i],
#if (P4EST_DIM)==3
                             xyz_on_f_m_lobatto[2][i],
#endif
                             bc_params->user
                            );
  }

  d4est_quadrature_interpolate
    (
     d4est_ops,
     d4est_quad,
     d4est_geom,
     &face_object,
     QUAD_OBJECT_MORTAR,
     QUAD_INTEGRAND_UNKNOWN,
     u_at_bndry_lobatto,
     e_m->deg,
     u_at_bndry_lobatto_to_quad,
     deg_mortar_quad
    );
  

  double u_m_on_f_m_min_u_at_bndry_quad;  
  for(int i = 0; i < face_nodes_m_quad; i++){

    if (bc_data->eval_method == EVAL_BNDRY_FCN_ON_LOBATTO){
       u_m_on_f_m_min_u_at_bndry_quad = u_m_on_f_m_quad[i] - u_at_bndry_lobatto_to_quad[i];
    }
    else if (bc_data->eval_method == EVAL_BNDRY_FCN_ON_QUAD){
      u_m_on_f_m_min_u_at_bndry_quad
        = u_m_on_f_m_quad[i] - bc_data->dirichlet_fcn
                            (
                             boundary_data->xyz_on_f_m_quad[0][i],
                             boundary_data->xyz_on_f_m_quad[1][i],
#if (P4EST_DIM)==3
                             boundary_data->xyz_on_f_m_quad[2][i],
#endif
                             bc_params->user
                            );
    }
    else {
      D4EST_ABORT("eval method for dirichlet bc should be set\n");
    }

    term1_quad[i] = 0.;
    for (int d = 0; d < (P4EST_DIM); d++){
      term1_quad[i] += -1.*n_on_f_m_quad[d][i]*sj_on_f_m_quad[i]
                       *(dudx_m_on_f_m_quad[d][i]);
    }
    for (int l = 0; l < (P4EST_DIM); l++){
      term2_quad[l][i] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        term2_quad[l][i] += -.5*drst_dxyz_quad[l][d][i]
                            *n_on_f_m_quad[d][i]
                            *sj_on_f_m_quad[i]
                            *2.*u_m_on_f_m_min_u_at_bndry_quad;        
      }
    }
    term3_quad[i] = sj_on_f_m_quad[i]
                    *sigma[i]
                    *u_m_on_f_m_min_u_at_bndry_quad;
    
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
     deg_mortar_quad,
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
       deg_mortar_quad,
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
     deg_mortar_quad,
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

  P4EST_FREE(ones_quad);
}


static void
d4est_laplacian_flux_sipg_dirichlet
(
 p4est_t* p4est,
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_laplacian_flux_boundary_data_t* boundary_data,
 void* boundary_condition_fcn_data,
 void* flux_parameter_data
)
{
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
  
  d4est_laplacian_flux_sipg_dirichlet_aux
    (
     p4est,
     e_m,
     f_m,
     mortar_side_id_m,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     boundary_data,
     boundary_condition_fcn_data,
     flux_parameter_data,
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
  
  int volume_nodes_m = d4est_lgl_get_nodes((P4EST_DIM), e_m->deg);

  for (int i = 0; i < volume_nodes_m; i++){
    for (int d = 0; d < (P4EST_DIM); d++){
      boundary_data->Au_m[i] += DT_lifted_VT_w_term2_lobatto[d][i];
      /* printf("lifted_VT_w_term2_lobatto[%d][%d] = %.15f\n",d,i, lifted_VT_w_term2_lobatto[d][i]); */
      /* e_m->Au_elem[i] += DT_lifted_VT_w_term2_lobatto[d][i]; */
    }
    /* e_m->Au_elem[i] += lifted_VT_w_term3_lobatto[i]; */
    boundary_data->Au_m[i] +=  lifted_VT_w_term3_lobatto[i];
    /* printf("lifted_VT_w_term3_lobatto[%d] = %.15f\n",i, lifted_VT_w_term3_lobatto[i]); */
    /* e_m->Au_elem[i] += lifted_VT_w_term1_lobatto[i]; */
    boundary_data->Au_m[i] += lifted_VT_w_term1_lobatto[i];
    /* printf("lifted_VT_w_term1_lobatto[%d] = %.15f\n",i, lifted_VT_w_term1_lobatto[i]); */
  }

  
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
d4est_laplacian_flux_sipg_robin_aux
(
 p4est_t* p4est,
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_laplacian_flux_boundary_data_t* boundary_data,
 void* boundary_condition_fcn_data,
 void* flux_parameter_data,
 double * restrict  term1_quad,
 double * restrict  VT_w_term1_lobatto,
 double * restrict  lifted_VT_w_term1_lobatto
){
  d4est_quadrature_mortar_t* face_object = boundary_data->face_object;
  int deg_mortar_quad = boundary_data->deg_mortar_quad;
  int face_nodes_m_lobatto = boundary_data->face_nodes_m_lobatto;
  int face_nodes_m_quad = boundary_data->face_nodes_m_quad;
  double* u_m_on_f_m_quad = boundary_data->u_m_on_f_m_quad;
  double* sj_on_f_m_quad = boundary_data->sj_on_f_m_quad;
  /* double* j_div_sj_quad = boundary_data->j_div_sj_quad; */
  double* xyz_on_f_m_lobatto [(P4EST_DIM)]; 
  double* drst_dxyz_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* dudx_m_on_f_m_quad [(P4EST_DIM)];  
  double* n_on_f_m_quad [(P4EST_DIM)];
  double* xyz_on_f_m_quad [(P4EST_DIM)];
  D4EST_COPY_DBYD_MAT(boundary_data->drst_dxyz_quad, drst_dxyz_quad);
  D4EST_COPY_DIM_VEC(boundary_data->dudx_m_on_f_m_quad, dudx_m_on_f_m_quad);
  D4EST_COPY_DIM_VEC(boundary_data->n_on_f_m_quad, n_on_f_m_quad);
  D4EST_COPY_DIM_VEC(boundary_data->xyz_on_f_m_lobatto, xyz_on_f_m_lobatto);
  D4EST_COPY_DIM_VEC(boundary_data->xyz_on_f_m_quad, xyz_on_f_m_quad);

  double* ones_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  d4est_util_fill_array(ones_quad, 1., face_nodes_m_quad);


  d4est_laplacian_robin_bc_t* bc = (d4est_laplacian_robin_bc_t*)boundary_condition_fcn_data;
  
  for(int i = 0; i < face_nodes_m_quad; i++){

    double robin_coeff_quad = bc->robin_coeff
                          (
                           xyz_on_f_m_quad[0][i],
                           xyz_on_f_m_quad[1][i],
#if (P4EST_DIM)==3
                           xyz_on_f_m_quad[2][i],
#endif
                           bc->user,
                           boundary_data,
                           i
                          );

    double robin_rhs_quad = bc->robin_rhs
                          (
                           xyz_on_f_m_quad[0][i],
                           xyz_on_f_m_quad[1][i],
#if (P4EST_DIM)==3
                           xyz_on_f_m_quad[2][i],
#endif
                           bc->user,
                           boundary_data,
                           i
                          );

    term1_quad[i] = sj_on_f_m_quad[i]
                    *(robin_coeff_quad*u_m_on_f_m_quad[i] - robin_rhs_quad);
    
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
     deg_mortar_quad,
     VT_w_term1_lobatto
    );

 
  d4est_operators_apply_lift(
                             d4est_ops,
                             VT_w_term1_lobatto,
                             (P4EST_DIM),
                             e_m->deg,
                             f_m,
                             lifted_VT_w_term1_lobatto);


  P4EST_FREE(ones_quad);
}

static void
d4est_laplacian_flux_sipg_robin
(
 p4est_t* p4est,
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_laplacian_flux_boundary_data_t* boundary_data,
 void* boundary_condition_fcn_data,
 void* flux_parameter_data
)
{
  int face_nodes_m_lobatto = boundary_data->face_nodes_m_lobatto;
  int face_nodes_m_quad = boundary_data->face_nodes_m_quad;
  int volume_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM), e_m->deg);
  
  double* term1_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  double* VT_w_term1_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* lifted_VT_w_term1_lobatto = P4EST_ALLOC(double, volume_nodes_m_lobatto);

  d4est_laplacian_flux_sipg_robin_aux
    (
     p4est,
     e_m,
     f_m,
     mortar_side_id_m,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     boundary_data,
     boundary_condition_fcn_data,
     flux_parameter_data,
     term1_quad,
     VT_w_term1_lobatto,
     lifted_VT_w_term1_lobatto
    );

  
  int volume_nodes_m = d4est_lgl_get_nodes((P4EST_DIM), e_m->deg);
  for (int i = 0; i < volume_nodes_m; i++){
    boundary_data->Au_m[i] += lifted_VT_w_term1_lobatto[i];
  }

  P4EST_FREE(term1_quad);
  P4EST_FREE(VT_w_term1_lobatto);
  P4EST_FREE(lifted_VT_w_term1_lobatto);
}




void
d4est_laplacian_flux_sipg_interface_aux
(
 p4est_t* p4est,
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
 d4est_mesh_data_t* d4est_factors,
 d4est_laplacian_flux_interface_data_t* mortar_data,
 void* params,
 double * restrict  lifted_proj_VT_w_term1_mortar_lobatto,
 double * restrict  proj_VT_w_term1_mortar_lobatto,
 double * restrict  VT_w_term1_mortar_lobatto,
 double * restrict  term1_mortar_quad,
 double *   DT_lifted_proj_VT_w_term2_mortar_lobatto [(P4EST_DIM)],
 double *   lifted_proj_VT_w_term2_mortar_lobatto [(P4EST_DIM)],
 double *   proj_VT_w_term2_mortar_lobatto [(P4EST_DIM)],
 double *   VT_w_term2_mortar_lobatto [(P4EST_DIM)],
 double *   term2_mortar_quad [(P4EST_DIM)],
 double * restrict  lifted_proj_VT_w_term3_mortar_lobatto,
 double * restrict  proj_VT_w_term3_mortar_lobatto,
 double * restrict  VT_w_term3_mortar_lobatto,
 double * restrict  term3_mortar_quad   
)
{
  d4est_quadrature_mortar_t* mortar_face_object = mortar_data->mortar_face_object;
  
  double * restrict  u_m_on_f_m_mortar_quad = mortar_data->u_m_on_f_m_mortar_quad;
  double * restrict  sj_on_f_m_mortar_quad = mortar_data->sj_on_f_m_mortar_quad;
  /* double * restrict  j_div_sj_on_f_m_mortar_quad = mortar_data->j_div_sj_on_f_m_mortar_quad; */
  double * restrict  u_p_on_f_p_mortar_quad = mortar_data->u_p_on_f_p_mortar_quad;
  /* double* j_div_sj_on_f_p_mortar_quad = mortar_data->j_div_sj_on_f_p_mortar_quad; */
  
  double* drst_dxyz_m_on_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_quad [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_quad [(P4EST_DIM)];
  double* n_on_f_m_mortar_quad [(P4EST_DIM)];

  D4EST_COPY_DBYD_MAT(mortar_data->drst_dxyz_m_on_mortar_quad, drst_dxyz_m_on_mortar_quad);
  D4EST_COPY_DIM_VEC(mortar_data->dudx_m_on_f_m_mortar_quad, dudx_m_on_f_m_mortar_quad);
  D4EST_COPY_DIM_VEC(mortar_data->dudx_p_on_f_p_mortar_quad, dudx_p_on_f_p_mortar_quad);
  D4EST_COPY_DIM_VEC(mortar_data->n_on_f_m_mortar_quad, n_on_f_m_mortar_quad);
  
  int faces_mortar = mortar_data->faces_mortar;
  int total_side_nodes_m_lobatto = mortar_data->total_side_nodes_m_lobatto;
  int total_nodes_mortar_lobatto = mortar_data->total_nodes_mortar_lobatto;
  int total_nodes_mortar_quad = mortar_data->total_nodes_mortar_quad;
  
  int* deg_mortar_quad = mortar_data->deg_mortar_quad;
  int* nodes_mortar_quad = mortar_data->nodes_mortar_quad;
  int* nodes_mortar_lobatto = mortar_data->nodes_mortar_lobatto;
  int* deg_mortar_lobatto = mortar_data->deg_mortar_lobatto;
  int* face_nodes_m_lobatto = mortar_data->face_nodes_m_lobatto;
  int* deg_m_lobatto = mortar_data->deg_m_lobatto;
  int* deg_p_lobatto = mortar_data->deg_p_lobatto;
  

  d4est_laplacian_flux_sipg_params_t* ip_flux_params = (d4est_laplacian_flux_sipg_params_t*) params;
    
  double* ones_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  d4est_util_fill_array(ones_mortar_quad, 1., total_nodes_mortar_quad);

  double * restrict  hm_mortar_quad = mortar_data->hm_mortar_quad;//P4EST_ALLOC(double, total_nodes_mortar_quad);
  double * restrict  hp_mortar_quad = mortar_data->hp_mortar_quad;//P4EST_ALLOC(double, total_nodes_mortar_quad);

  double debug_sigma_sum = 0.;
  int stride = 0;
  double* sigma = P4EST_ALLOC(double, total_nodes_mortar_quad);
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_quad[f]; k++){
      int ks = k + stride;
        sigma[ks] = ip_flux_params->sipg_penalty_fcn
                    (
                     (faces_m == faces_mortar) ? deg_m_lobatto[f] : deg_m_lobatto[0],
                     hm_mortar_quad[ks],
                     (faces_p == faces_mortar) ? deg_p_lobatto[f] : deg_p_lobatto[0],
                     hp_mortar_quad[ks],
                     ip_flux_params->sipg_penalty_prefactor
                    );
        debug_sigma_sum += sigma[ks];
    }
    stride += nodes_mortar_quad[f];
  }


  /* P4EST_FREE(hm_mortar_quad); */
  /* P4EST_FREE(hp_mortar_quad); */

  double* sj_n_du [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(sj_n_du, total_nodes_mortar_quad);
  
  stride = 0;
  int stride_lobatto = 0;
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_quad[f]; k++){
      int ks = k + stride;

      term1_mortar_quad[ks] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        term1_mortar_quad[ks] += -1.*n_on_f_m_mortar_quad[d][ks]*sj_on_f_m_mortar_quad[ks]
                                 *.5*(dudx_p_on_f_p_mortar_quad[d][ks] + dudx_m_on_f_m_mortar_quad[d][ks]);
      }
      
      /* printf("term1_mortar_quad[%d] = %.15f\n", ks, term1_mortar_quad[ks]); */

    for (int d = 0; d < (P4EST_DIM); d++){
        sj_n_du[d][ks] = sj_on_f_m_mortar_quad[ks]
                         *n_on_f_m_mortar_quad[d][ks]
                         *(u_m_on_f_m_mortar_quad[ks] - u_p_on_f_p_mortar_quad[ks]);
          
      }
      
      for (int l = 0; l < (P4EST_DIM); l++){
        term2_mortar_quad[l][ks] = 0.;
        for (int d = 0; d < (P4EST_DIM); d++){
          term2_mortar_quad[l][ks] += -.5*drst_dxyz_m_on_mortar_quad[l][d][ks]
                                      *sj_on_f_m_mortar_quad[ks]
                                      *n_on_f_m_mortar_quad[d][ks]
                                      *(u_m_on_f_m_mortar_quad[ks] - u_p_on_f_p_mortar_quad[ks]);
          
        }
        /* printf("term2_mortar_quad[%d] = %.15f\n", ks, term2_mortar_quad[l][ks]); */
      }
        
      term3_mortar_quad[ks] = sj_on_f_m_mortar_quad[ks]*sigma[ks]*(u_m_on_f_m_mortar_quad[ks] - u_p_on_f_p_mortar_quad[ks]);
      /* printf("term3_mortar_quad[%d] = %.15f\n", ks, term3_mortar_quad[ks]); */
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

  
    for (int d = 0; d < (P4EST_DIM); d++){
      
      d4est_quadrature_apply_galerkin_integral
        (
         d4est_ops,
         d4est_geom,
         d4est_quad,
         mortar_face_object,
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
       mortar_face_object,
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


    /* if (e_m[0]->id == 11 && f_m == 1){ */
      /* printf("orientation = %d\n", orientation); */
      /* DEBUG_PRINT_ARR_DBL(sj_n_du[0], total_nodes_mortar_quad); */
      /* DEBUG_PRINT_ARR_DBL(sj_n_du[1], total_nodes_mortar_quad); */
      /* DEBUG_PRINT_ARR_DBL(sj_n_du[2], total_nodes_mortar_quad); */
    /* } */

    D4EST_FREE_DIM_VEC(sj_n_du);
  

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

  /* if(e_m[0]->id == 8 && f_m == 2){ */
  /* if(e_m[0]->id == 11 && f_m == 1){ */
  /*   /\* printf("with_opt\n"); *\/ */
  /*   printf("total_side_nodes_p_lobatto = %d\n", total_side_nodes_m_lobatto); */
    /* printf("f_m = %d\n", f_m); */
  /*   DEBUG_PRINT_ARR_DBL(VT_w_term2_mortar_lobatto[0], total_nodes_mortar_lobatto); */
  /*   DEBUG_PRINT_ARR_DBL(VT_w_term2_mortar_lobatto[1], total_nodes_mortar_lobatto); */
  /*   DEBUG_PRINT_ARR_DBL(VT_w_term2_mortar_lobatto[2], total_nodes_mortar_lobatto); */
    
    /* DEBUG_PRINT_ARR_DBL(proj_VT_w_term2_mortar_lobatto[0], total_side_nodes_m_lobatto); */
    /* DEBUG_PRINT_ARR_DBL(proj_VT_w_term2_mortar_lobatto[1], total_side_nodes_m_lobatto); */
    /* DEBUG_PRINT_ARR_DBL(proj_VT_w_term2_mortar_lobatto[2], total_side_nodes_m_lobatto); */
  /* } */
  
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

  int volume_stride = 0;
  stride = 0;
  for (int f = 0; f < faces_m; f++){
    int volume_nodes_m = d4est_lgl_get_nodes((P4EST_DIM), e_m[f]->deg);
    if (e_m_is_ghost[f] == 0){

      d4est_operators_apply_lift(
                                 d4est_ops,
                                 &proj_VT_w_term1_mortar_lobatto[stride],
                                 (P4EST_DIM),
                                 e_m[f]->deg,
                                 f_m,
                                 &lifted_proj_VT_w_term1_mortar_lobatto[volume_stride]);


      for (int d = 0; d < (P4EST_DIM); d++){
        d4est_operators_apply_lift(
                                   d4est_ops,
                                   &proj_VT_w_term2_mortar_lobatto[d][stride],
                                   (P4EST_DIM),
                                   e_m[f]->deg,
                                   f_m,
                                   &lifted_proj_VT_w_term2_mortar_lobatto[d][volume_stride]);

        d4est_operators_apply_dij_transpose(d4est_ops,
                                            &lifted_proj_VT_w_term2_mortar_lobatto[d][volume_stride],
                                            (P4EST_DIM),
                                            e_m[f]->deg,
                                            d,
                                            &DT_lifted_proj_VT_w_term2_mortar_lobatto[d][volume_stride]
                                           );

        
        if (faces_m != faces_mortar){
            d4est_linalg_vec_scale(.5, &DT_lifted_proj_VT_w_term2_mortar_lobatto[d][volume_stride], volume_nodes_m);
        }
        

      }
        
      d4est_operators_apply_lift(
                                 d4est_ops,
                                 &proj_VT_w_term3_mortar_lobatto[stride],
                                 (P4EST_DIM),
                                 e_m[f]->deg,
                                 f_m,
                                 &lifted_proj_VT_w_term3_mortar_lobatto[volume_stride]);


      
    }

    volume_stride += volume_nodes_m;
    stride += face_nodes_m_lobatto[f];
  }
  P4EST_FREE(ones_mortar_quad);
  P4EST_FREE(sigma);
  
}

static void
d4est_laplacian_flux_sipg_interface
(
 p4est_t* p4est,
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
 d4est_mesh_data_t* d4est_factors,
 d4est_laplacian_flux_interface_data_t* mortar_data,
 void* params
)
{
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

  d4est_laplacian_flux_sipg_interface_aux
    (
     p4est,
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
     d4est_factors,
     mortar_data,
     params,
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
    
  int stride = 0;
  for (int f = 0; f < faces_m; f++){
    int volume_nodes_m = d4est_lgl_get_nodes((P4EST_DIM), e_m[f]->deg);
    /* this if statement could perhaps be reduced to just Au_m[f] != NULL */
    if (e_m_is_ghost[f] == 0 &&  mortar_data->Au_m[f] != NULL){    
      /* printf("Element = %d\n, face = %d\n", e_m[f]->id, f_m); */
      for (int i = 0; i < volume_nodes_m; i++){
        /* printf("term 2\n"); */
        for (int d = 0; d < (P4EST_DIM); d++){
          mortar_data->Au_m[f][i] +=  DT_lifted_proj_VT_w_term2_mortar_lobatto[d][i + stride];
        }

        /* printf("term 3\n"); */
        mortar_data->Au_m[f][i] += lifted_proj_VT_w_term3_mortar_lobatto[i + stride];

        /* printf("term 1\n"); */
        mortar_data->Au_m[f][i] += lifted_proj_VT_w_term1_mortar_lobatto[i + stride];
        /* if(e_m[0]->tree != e_p[0]->tree) */
          /* printf("mortar_data->Au_m[%d][%d] = %.15f\n", f,i, mortar_data->Au_m[f][i]); */
      }
    }
    stride += volume_nodes_m;
  }

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


static double
d4est_laplacian_flux_sipg_penalty_meanp_sqr_over_meanh
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
d4est_laplacian_flux_sipg_penalty_mean_p_sqr_over_h
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double mean_penalty= .5*(deg_m*deg_m/h_m + deg_p*deg_p/h_p);
  return (penalty_prefactor*mean_penalty);
}



static double
d4est_laplacian_flux_sipg_penalty_maxp_sqr_over_minh
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
d4est_laplacian_flux_sipg_params_get_string_from_penalty_fcn
(
 penalty_calc_t fcn,
 char* string
)
{
  if (fcn == d4est_laplacian_flux_sipg_penalty_maxp_sqr_over_minh){
    strcpy(string,"maxp_sqr_over_minh");
  }
  else if (fcn == d4est_laplacian_flux_sipg_penalty_meanp_sqr_over_meanh){
    strcpy(string,"meanp_sqr_over_meanh");
  }
  else if (fcn == d4est_laplacian_flux_sipg_penalty_mean_p_sqr_over_h){
    strcpy(string,"mean_p_sqr_over_h");
  }
  else {
    strcpy(string,"not_set");
  }
}

static penalty_calc_t
d4est_laplacian_flux_sipg_get_penalty_fcn_from_string
(
 const char* string
)
{
  if (d4est_util_match(string,"maxp_sqr_over_minh")){
    return d4est_laplacian_flux_sipg_penalty_maxp_sqr_over_minh;
  }
  else if (d4est_util_match(string,"meanp_sqr_over_meanh")){
    return d4est_laplacian_flux_sipg_penalty_meanp_sqr_over_meanh;
  }
  else if (d4est_util_match(string,"mean_p_sqr_over_h")){
    return d4est_laplacian_flux_sipg_penalty_mean_p_sqr_over_h;
  }
  else {
    D4EST_ABORT("This ip flux penalty calculation fcn does not exist");
    return NULL;
  }
}

static
int d4est_laplacian_flux_sipg_params_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_laplacian_flux_sipg_params_t* pconfig = (d4est_laplacian_flux_sipg_params_t*)user;

/* typedef enum { H_EQ_J_DIV_SJ_QUAD, H_EQ_J_DIV_SJ_MIN_LOBATTO, H_EQ_J_DIV_SJ_MEAN_LOBATTO, H_EQ_TREE_H, H_EQ_VOLUME_DIV_AREA, H_EQ_FACE_DIAM, H_EQ_TOTAL_VOLUME_DIV_TOTAL_AREA, H_EQ_NOTSET } h_calc_method_t;
  */
  if (d4est_util_match_couple(section,"flux",name,"sipg_penalty_prefactor")) {
    D4EST_ASSERT(pconfig->sipg_penalty_prefactor == -1);
    pconfig->sipg_penalty_prefactor = atof(value);
  }
  else if (d4est_util_match_couple(section,"flux",name,"sipg_penalty_fcn")) {
    D4EST_ASSERT(pconfig->sipg_penalty_fcn == NULL);
    pconfig->sipg_penalty_fcn = d4est_laplacian_flux_sipg_get_penalty_fcn_from_string(value);
  }
 
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void
d4est_laplacian_flux_sipg_params_input
(
 p4est_t* p4est,
 const char* printf_prefix,
 const char* input_file,
 d4est_laplacian_flux_sipg_params_t* input
)
{
  input->sipg_penalty_fcn = NULL;
  input->sipg_penalty_prefactor = -1.;

  if (ini_parse(input_file, d4est_laplacian_flux_sipg_params_input_handler, input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("flux", input->sipg_penalty_prefactor, -1);
  D4EST_CHECK_INPUT("flux", input->sipg_penalty_fcn, NULL);
  
  char penalty_calculate_fcn [50];
  char h_eq [50];

  d4est_laplacian_flux_sipg_params_get_string_from_penalty_fcn (input->sipg_penalty_fcn,penalty_calculate_fcn);


  double check_function = input->sipg_penalty_fcn(1,.1,1,.1,0.);

  D4EST_ASSERT(check_function == 0.);
  
  if(p4est->mpirank == 0){
    printf("%s: d4est_laplacian_flux_sipg_params_penalty_prefactor = %f\n", printf_prefix, input->sipg_penalty_prefactor);
    printf("%s: d4est_laplacian_flux_sipg_params_penalty_calculate_fcn = %s\n", printf_prefix, penalty_calculate_fcn);
  }
}

void
d4est_laplacian_flux_sipg_params_new
(
 p4est_t* p4est,
 const char* print_prefix,
 const char* input_file,
 d4est_laplacian_flux_data_t* d4est_laplacian_flux_data
)
{
  d4est_laplacian_flux_sipg_params_t* d4est_laplacian_flux_sipg_params = P4EST_ALLOC(d4est_laplacian_flux_sipg_params_t, 1); 
  d4est_laplacian_flux_sipg_params_input(p4est, print_prefix, input_file, d4est_laplacian_flux_sipg_params);

  d4est_laplacian_flux_sipg_params->size_params = NULL;
  
  d4est_laplacian_flux_data->flux_data = d4est_laplacian_flux_sipg_params;
  d4est_laplacian_flux_data->interface_fcn = d4est_laplacian_flux_sipg_interface;
  /* d4est_laplacian_flux_data->interface_fcn = NULL; */
  /* d4est_laplacian_flux_data->boundary_fcn = NULL; */

  if (d4est_laplacian_flux_data->bc_type == BC_DIRICHLET){
    d4est_laplacian_flux_data->boundary_fcn = d4est_laplacian_flux_sipg_dirichlet;
  }
  else if (d4est_laplacian_flux_data->bc_type == BC_ROBIN){
    d4est_laplacian_flux_data->boundary_fcn = d4est_laplacian_flux_sipg_robin;
  }
  else {
    D4EST_ABORT("Not a supported boundary condition");
  }

  d4est_laplacian_flux_data->destroy = d4est_laplacian_flux_sipg_params_destroy;
}

void
d4est_laplacian_flux_sipg_params_destroy
(
 d4est_laplacian_flux_data_t* data
){
  P4EST_FREE(data->flux_data);
  P4EST_FREE(data);
}


