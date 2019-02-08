#define _GNU_SOURCE
#include <pXest.h>
#include <d4est_util.h>
#include <ini.h>
#include <d4est_operators.h>
#include <d4est_element_data.h>
#include <d4est_linalg.h>
#include <d4est_laplacian_with_opt_flux_sipg.h>
#include <d4est_elliptic_data.h>
#include <d4est_mortars_with_opt.h>
#include <d4est_quadrature_compactified.h>
#include <d4est_quadrature_lobatto.h>
#include <d4est_laplacian_with_opt_flux.h>

void
d4est_laplacian_with_opt_flux_sipg_dirichlet_aux
(
 p4est_t* p4est,
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_laplacian_with_opt_flux_boundary_data_t* boundary_data,
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
 double * restrict  u_m_on_f_m_min_u_at_bndry_quad
){
  d4est_laplacian_with_opt_dirichlet_bc_t* bc_data = boundary_condition_fcn_data;
  
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
  D4EST_COPY_DBYD_MAT(boundary_data->drst_dxyz_quad, drst_dxyz_quad);
  D4EST_COPY_DIM_VEC(boundary_data->dudx_m_on_f_m_quad, dudx_m_on_f_m_quad);
  D4EST_COPY_DIM_VEC(boundary_data->n_on_f_m_quad, n_on_f_m_quad);
  D4EST_COPY_DIM_VEC(boundary_data->xyz_on_f_m_lobatto, xyz_on_f_m_lobatto);
  
  d4est_laplacian_with_opt_flux_sipg_params_t* ip_flux_params = (d4est_laplacian_with_opt_flux_sipg_params_t*) flux_parameter_data;

  d4est_laplacian_with_opt_dirichlet_bc_t* bc_params = boundary_condition_fcn_data;
  
  double* ones_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  d4est_util_fill_array(ones_quad, 1., face_nodes_m_quad);


  double* h_quad = boundary_data->h_quad;//P4EST_ALLOC(double, face_nodes_m_quad);

  d4est_mesh_size_parameters_t size_params = (ip_flux_params->size_params == NULL) ? d4est_mesh_get_size_parameters(d4est_factors) : *ip_flux_params->size_params;

  double multiplier = d4est_laplacian_with_opt_flux_sipg_get_multiplier(flux_parameter_data,e_m->region,e_m->region);
  
  for (int i = 0; i < face_nodes_m_quad; i++){

    sigma[i] = multiplier*ip_flux_params->sipg_penalty_fcn
               (
                e_m->deg,
                h_quad[i],
                e_m->deg,
                h_quad[i],
                ip_flux_params->sipg_penalty_prefactor
               );

  }

  if (bc_data->dirichlet_fcn != zero_fcn){
    
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
     u_m_on_f_m_min_u_at_bndry_quad,
     deg_mortar_quad
    );


  }  

   if (bc_data->eval_method == EVAL_BNDRY_FCN_ON_LOBATTO && bc_data->dirichlet_fcn != zero_fcn){
      
     for (int i = 0; i < face_nodes_m_quad; i++){
       u_m_on_f_m_min_u_at_bndry_quad[i] =  -u_m_on_f_m_min_u_at_bndry_quad[i] + u_m_on_f_m_quad[i];
     } 
    }
    else if (bc_data->eval_method == EVAL_BNDRY_FCN_ON_QUAD && bc_data->dirichlet_fcn != zero_fcn){
     for (int i = 0; i < face_nodes_m_quad; i++){
      u_m_on_f_m_min_u_at_bndry_quad[i]
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
    }
    else if (bc_data->dirichlet_fcn == zero_fcn){

      for(int i = 0; i < face_nodes_m_quad; i++){

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
                                *2.*u_m_on_f_m_quad[i];        
          }
        }
        term3_quad[i] = sj_on_f_m_quad[i]
                        *sigma[i]
                        *u_m_on_f_m_quad[i];
    
      }
  
    }
    else {
      D4EST_ABORT("Not a supported boundary sipg option");
    }

   if (bc_data->dirichlet_fcn != zero_fcn){
     for(int i = 0; i < face_nodes_m_quad; i++){

 
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
                               *2.*u_m_on_f_m_min_u_at_bndry_quad[i];        
         }
       }
       term3_quad[i] = sj_on_f_m_quad[i]
                       *sigma[i]
                       *u_m_on_f_m_min_u_at_bndry_quad[i];
    
     }
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
d4est_laplacian_with_opt_flux_sipg_dirichlet
(
 p4est_t* p4est,
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_laplacian_with_opt_flux_boundary_data_t* boundary_data,
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
  
  d4est_laplacian_with_opt_flux_sipg_dirichlet_aux
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
      /* e_m->Au_elem[i] += DT_lifted_VT_w_term2_lobatto[d][i]; */
    }
    /* e_m->Au_elem[i] += lifted_VT_w_term3_lobatto[i]; */
    boundary_data->Au_m[i] +=  lifted_VT_w_term3_lobatto[i];
    /* e_m->Au_elem[i] += lifted_VT_w_term1_lobatto[i]; */
    boundary_data->Au_m[i] += lifted_VT_w_term1_lobatto[i];
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
d4est_laplacian_with_opt_flux_sipg_robin_aux
(
 p4est_t* p4est,
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_laplacian_with_opt_flux_boundary_data_t* boundary_data,
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


  d4est_laplacian_with_opt_robin_bc_t* bc = (d4est_laplacian_with_opt_robin_bc_t*)boundary_condition_fcn_data;
  
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
d4est_laplacian_with_opt_flux_sipg_robin
(
 p4est_t* p4est,
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_laplacian_with_opt_flux_boundary_data_t* boundary_data,
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

  d4est_laplacian_with_opt_flux_sipg_robin_aux
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
d4est_laplacian_with_opt_flux_sipg_interface_aux
(
 p4est_t* p4est,
 d4est_element_data_t** e_m,
 int faces_m,
 int f_m,
 int mortar_side_id_m,
 int* e_m_is_ghost,
 d4est_element_data_t** e_p,
 int faces_p,
 int f_p,
 int mortar_side_id_p,
 int* e_p_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_laplacian_with_opt_flux_interface_data_t* mortar_data,
 void* params,
 double * restrict  lifted_proj_VT_w_term1_mortar_lobatto_m,
 double * restrict  proj_VT_w_term1_mortar_lobatto_m,
 double * restrict lifted_proj_VT_w_term1_mortar_lobatto_p,
 double * restrict  proj_VT_w_term1_mortar_lobatto_p,
 double * restrict  VT_w_term1_mortar_lobatto,
 double * restrict  term1_mortar_quad,
 double * DT_lifted_proj_VT_w_term2_mortar_lobatto_m [(P4EST_DIM)],
 double * lifted_proj_VT_w_term2_mortar_lobatto_m [(P4EST_DIM)],
 double * proj_VT_w_term2_mortar_lobatto_m [(P4EST_DIM)],
 double * DT_lifted_proj_VT_w_term2_mortar_lobatto_p [(P4EST_DIM)],
 double * lifted_proj_VT_w_term2_mortar_lobatto_p [(P4EST_DIM)],
 double * proj_VT_w_term2_mortar_lobatto_p [(P4EST_DIM)],
 double * VT_w_term2_mortar_lobatto [(P4EST_DIM)],
 double * term2_mortar_quad [(P4EST_DIM)],
 double * restrict  lifted_proj_VT_w_term3_mortar_lobatto_m,
 double * restrict  proj_VT_w_term3_mortar_lobatto_m,
 double * restrict  lifted_proj_VT_w_term3_mortar_lobatto_p,
 double * restrict  proj_VT_w_term3_mortar_lobatto_p,
 double * restrict  VT_w_term3_mortar_lobatto,
 double * restrict  term3_mortar_quad   
)
{

  d4est_laplacian_with_opt_flux_sipg_params_t* ip_flux_params = (d4est_laplacian_with_opt_flux_sipg_params_t*) params;
  
  d4est_quadrature_mortar_t* mortar_face_object = mortar_data->mortar_face_object;
  
  double * restrict  u_m_on_f_m_mortar_quad = mortar_data->u_m_on_f_m_mortar_quad;
  double * restrict  sj_on_f_m_mortar_quad = mortar_data->sj_on_f_m_mortar_quad;
  /* double * restrict  j_div_sj_on_f_m_mortar_quad = mortar_data->j_div_sj_on_f_m_mortar_quad; */
  double * restrict  u_p_on_f_p_mortar_quad = mortar_data->u_p_on_f_p_mortar_quad;
  /* double* j_div_sj_on_f_p_mortar_quad = mortar_data->j_div_sj_on_f_p_mortar_quad; */
  
  double* drst_dxyz_m_on_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_p_on_mortar_quad_porder [(P4EST_DIM)][(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_quad [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_quad [(P4EST_DIM)];
  double* n_on_f_m_mortar_quad [(P4EST_DIM)];

  D4EST_COPY_DBYD_MAT(mortar_data->drst_dxyz_m_on_mortar_quad, drst_dxyz_m_on_mortar_quad);
  D4EST_COPY_DBYD_MAT(mortar_data->drst_dxyz_p_on_mortar_quad_porder, drst_dxyz_p_on_mortar_quad_porder);
  D4EST_COPY_DIM_VEC(mortar_data->dudx_m_on_f_m_mortar_quad, dudx_m_on_f_m_mortar_quad);
  D4EST_COPY_DIM_VEC(mortar_data->dudx_p_on_f_p_mortar_quad, dudx_p_on_f_p_mortar_quad);
  D4EST_COPY_DIM_VEC(mortar_data->n_on_f_m_mortar_quad, n_on_f_m_mortar_quad);
  
  int faces_mortar = mortar_data->faces_mortar;
  int total_side_nodes_m_lobatto = mortar_data->total_side_nodes_m_lobatto;
  int total_side_nodes_p_lobatto = mortar_data->total_side_nodes_p_lobatto;
  int total_nodes_mortar_lobatto = mortar_data->total_nodes_mortar_lobatto;
  int total_nodes_mortar_quad = mortar_data->total_nodes_mortar_quad;
  
  int* deg_mortar_quad = mortar_data->deg_mortar_quad;
  int* deg_mortar_quad_porder = mortar_data->deg_mortar_quad_porder;
  int* nodes_mortar_quad = mortar_data->nodes_mortar_quad;
  int* nodes_mortar_quad_porder = mortar_data->nodes_mortar_quad_porder;
  int* nodes_mortar_lobatto = mortar_data->nodes_mortar_lobatto;
  int* deg_mortar_lobatto = mortar_data->deg_mortar_lobatto;
  int* deg_mortar_lobatto_porder = mortar_data->deg_mortar_lobatto_porder;
  int* face_nodes_m_lobatto = mortar_data->face_nodes_m_lobatto;
  int* deg_m_lobatto = mortar_data->deg_m_lobatto;
  int* deg_p_lobatto = mortar_data->deg_p_lobatto;

  double* ones_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  d4est_util_fill_array(ones_mortar_quad, 1., total_nodes_mortar_quad);

  double * restrict  hm_mortar_quad = mortar_data->hm_mortar_quad;//P4EST_ALLOC(double, total_nodes_mortar_quad);
  double * restrict  hp_mortar_quad = mortar_data->hp_mortar_quad;//P4EST_ALLOC(double, total_nodes_mortar_quad);

  /* double debug_sigma_sum = 0.; */

  double multiplier = d4est_laplacian_with_opt_flux_sipg_get_multiplier(ip_flux_params,e_m[0]->region,e_p[0]->region);
  int stride = 0;
  double* sigma = P4EST_ALLOC(double, total_nodes_mortar_quad);
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_quad[f]; k++){
      int ks = k + stride;
        sigma[ks] = multiplier*ip_flux_params->sipg_penalty_fcn
                    (
                     (faces_m == faces_mortar) ? deg_m_lobatto[f] : deg_m_lobatto[0],
                     hm_mortar_quad[ks],
                     (faces_p == faces_mortar) ? deg_p_lobatto[f] : deg_p_lobatto[0],
                     hp_mortar_quad[ks],
                     ip_flux_params->sipg_penalty_prefactor
                    );


        /* debug_sigma_sum += sigma[ks]; */
    }
    stride += nodes_mortar_quad[f];
  }


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
          term2_mortar_quad[l][ks] += -.5*drst_dxyz_m_on_mortar_quad[l][d][ks]*sj_n_du[d][ks];
                                      /* *sj_on_f_m_mortar_quad[ks] */
                                      /* *n_on_f_m_mortar_quad[d][ks] */
                                      /* *(u_m_on_f_m_mortar_quad[ks] - u_p_on_f_p_mortar_quad[ks]); */
          
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

  /* if(e_p[0]->id == 8 && f_p == 2){ */
  /*   printf("with_opt\n"); */
  /*   printf("orientation = %d\n", orientation); */
  /*   DEBUG_PRINT_ARR_DBL(sj_n_du[0], total_nodes_mortar_quad); */
  /*   DEBUG_PRINT_ARR_DBL(sj_n_du[1], total_nodes_mortar_quad); */
  /*   DEBUG_PRINT_ARR_DBL(sj_n_du[2], total_nodes_mortar_quad); */
  /* }       */
  
  int is_it_tree_boundary = (e_m[0]->tree != e_p[0]->tree);
  double* term2_mortar_quad_p_porder [P4EST_DIM];
  double* VT_w_term2_mortar_lobatto_porder [P4EST_DIM];
  if(is_it_tree_boundary){    
    double* sj_n_du_porder[(P4EST_DIM)];
    for (int d = 0; d < (P4EST_DIM); d++){
      term2_mortar_quad_p_porder[d] = P4EST_ALLOC(double, total_nodes_mortar_quad);
      VT_w_term2_mortar_lobatto_porder[d] = P4EST_ALLOC(double, total_nodes_mortar_quad);
      if (orientation != 0){
        sj_n_du_porder[d] = d4est_mortars_with_opt_reorient_if_needed
                            (
                             d4est_ops,
                             e_p,
                             e_m,
                             orientation,
                             f_p,
                             f_m,
                             sj_n_du[d],
                             deg_mortar_quad_porder,
                             deg_mortar_quad,
                             total_nodes_mortar_quad,
                             faces_mortar
                            );
      }
      else {
        sj_n_du_porder[d] = sj_n_du[d];
      }
    }
    /* for (int i = 0; i < faces_p; i++){ */
    /*   printf("e_p[%d] id = %d, face = %d, e_m[0] id = %d\n", i, e_p[i]->id, f_p, e_m[0]->id); */
    /* } */

    /* if (e_p[0]->id == 11 && f_p == 1){ */
    /*   printf("e_m[0]->id = %d, faces_m =%d\n", e_m[0]->id, faces_m); */
    /*   printf("orientation = %d\n", orientation); */
    /*   DEBUG_PRINT_ARR_DBL(sj_n_du[0], total_nodes_mortar_quad); */
    /*   DEBUG_PRINT_ARR_DBL(sj_n_du[1], total_nodes_mortar_quad); */
    /*   DEBUG_PRINT_ARR_DBL(sj_n_du[2], total_nodes_mortar_quad); */
      
    /*   DEBUG_PRINT_ARR_DBL(sj_n_du_porder[0], total_nodes_mortar_quad); */
    /*   DEBUG_PRINT_ARR_DBL(sj_n_du_porder[1], total_nodes_mortar_quad); */
    /*   DEBUG_PRINT_ARR_DBL(sj_n_du_porder[2], total_nodes_mortar_quad); */
    /* } */

    stride = 0;
    stride_lobatto = 0;
    for (int f = 0; f < faces_mortar; f++){
      for (int k = 0; k < nodes_mortar_quad_porder[f]; k++){
        int ks = k + stride;

        for (int l = 0; l < (P4EST_DIM); l++){
          term2_mortar_quad_p_porder[l][ks] = 0.;
          for (int d = 0; d < (P4EST_DIM); d++){
            term2_mortar_quad_p_porder[l][ks] += -.5*drst_dxyz_p_on_mortar_quad_porder[l][d][ks]*sj_n_du_porder[d][ks];
          }
        }
      }
      for (int d = 0; d < (P4EST_DIM); d++){      
        d4est_quadrature_apply_galerkin_integral
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           mortar_data->mortar_face_object_porder,
           QUAD_OBJECT_MORTAR,
           QUAD_INTEGRAND_UNKNOWN,
           &term2_mortar_quad_p_porder[d][stride],
           deg_mortar_lobatto_porder[f],
           &ones_mortar_quad[stride],
           deg_mortar_quad_porder[f],
           &VT_w_term2_mortar_lobatto_porder[d][stride_lobatto]
          );
      }
      stride += nodes_mortar_quad_porder[f];
      stride_lobatto += mortar_data->nodes_mortar_lobatto_porder[f];
    }

    for (int d = 0; d < (P4EST_DIM); d++){
      P4EST_FREE(term2_mortar_quad_p_porder[d]);
      if (orientation != 0){
        P4EST_FREE(sj_n_du_porder[d]);
      }
    }
    
  }

  int m_side_matches_mortar = mortar_data->m_side_matches_mortar;
  int p_side_matches_mortar = mortar_data->p_side_matches_mortar;
  
  for (int d = 0; d < (P4EST_DIM); d++){

      if (!m_side_matches_mortar){
        d4est_mortars_with_opt_project_mass_mortar_onto_side
          (
           d4est_ops,
           VT_w_term2_mortar_lobatto[d],
           faces_mortar,
           deg_mortar_lobatto,
           proj_VT_w_term2_mortar_lobatto_m[d],
           faces_m,
           deg_m_lobatto
          );
      }
      if (!p_side_matches_mortar){
        d4est_mortars_with_opt_project_mass_mortar_onto_side
          (
           d4est_ops,
           (is_it_tree_boundary) ? VT_w_term2_mortar_lobatto_porder[d] : VT_w_term2_mortar_lobatto[d], //VT_w_term2_mortar_lobatto[d],
           faces_mortar,
           (is_it_tree_boundary) ? deg_mortar_lobatto_porder : deg_mortar_lobatto,
           proj_VT_w_term2_mortar_lobatto_p[d],
           faces_p,
           (is_it_tree_boundary) ? mortar_data->deg_p_lobatto_porder : deg_p_lobatto
          );
      }
  }


  if (!m_side_matches_mortar){
  d4est_mortars_with_opt_project_mass_mortar_onto_side
    (
     d4est_ops,
     VT_w_term1_mortar_lobatto,
     faces_mortar,
     deg_mortar_lobatto,
     proj_VT_w_term1_mortar_lobatto_m,
     faces_m,
     deg_m_lobatto
    );
  }
  if (!p_side_matches_mortar){
 d4est_mortars_with_opt_project_mass_mortar_onto_side
    (
     d4est_ops,
     VT_w_term1_mortar_lobatto,
     faces_mortar,
     deg_mortar_lobatto,
     proj_VT_w_term1_mortar_lobatto_p,
     faces_p,
     deg_p_lobatto
    );
  }

  if (!m_side_matches_mortar){
    d4est_mortars_with_opt_project_mass_mortar_onto_side
      (
       d4est_ops,
       VT_w_term3_mortar_lobatto,
       faces_mortar,
       deg_mortar_lobatto,
       proj_VT_w_term3_mortar_lobatto_m,
       faces_m,
       deg_m_lobatto
      );
  }
  if (!p_side_matches_mortar){
    d4est_mortars_with_opt_project_mass_mortar_onto_side
      (
       d4est_ops,
       VT_w_term3_mortar_lobatto,
       faces_mortar,
       deg_mortar_lobatto,
       proj_VT_w_term3_mortar_lobatto_p,
       faces_p,
       deg_p_lobatto /* in the order of (-) side, so this will work */
      );
  }

  int volume_stride = 0;
  stride = 0;
  for (int f = 0; f < faces_m; f++){
    int volume_nodes_m = d4est_lgl_get_nodes((P4EST_DIM), e_m[f]->deg);
    if (e_m_is_ghost[f] == 0){

      d4est_operators_apply_lift(
                                 d4est_ops,
                                 (m_side_matches_mortar) ? &VT_w_term1_mortar_lobatto[stride] :&proj_VT_w_term1_mortar_lobatto_m[stride],
                                 (P4EST_DIM),
                                 e_m[f]->deg,
                                 f_m,
                                 &lifted_proj_VT_w_term1_mortar_lobatto_m[volume_stride]);


      for (int d = 0; d < (P4EST_DIM); d++){
        d4est_operators_apply_lift(
                                   d4est_ops,
                                   (m_side_matches_mortar) ? &VT_w_term2_mortar_lobatto[d][stride] : &proj_VT_w_term2_mortar_lobatto_m[d][stride],
                                   (P4EST_DIM),
                                   e_m[f]->deg,
                                   f_m,
                                   &lifted_proj_VT_w_term2_mortar_lobatto_m[d][volume_stride]);

        d4est_operators_apply_dij_transpose(d4est_ops,
                                            &lifted_proj_VT_w_term2_mortar_lobatto_m[d][volume_stride],
                                            (P4EST_DIM),
                                            e_m[f]->deg,
                                            d,
                                            &DT_lifted_proj_VT_w_term2_mortar_lobatto_m[d][volume_stride]
                                           );

        
        if (faces_m != faces_mortar){
            d4est_linalg_vec_scale(.5, &DT_lifted_proj_VT_w_term2_mortar_lobatto_m[d][volume_stride], volume_nodes_m);
        }
        

      }
        
      d4est_operators_apply_lift(
                                 d4est_ops,
                                 (m_side_matches_mortar) ? &VT_w_term3_mortar_lobatto[stride] :&proj_VT_w_term3_mortar_lobatto_m[stride],
                                 (P4EST_DIM),
                                 e_m[f]->deg,
                                 f_m,
                                 &lifted_proj_VT_w_term3_mortar_lobatto_m[volume_stride]);


     
    }

    volume_stride += volume_nodes_m;
    stride += face_nodes_m_lobatto[f];
  }

  /* p-side  */


  volume_stride = 0;
  stride = 0;
  int deg_p_lobatto_porder [(P4EST_HALF)];
  for (int i = 0; i < faces_p; i++){
    deg_p_lobatto_porder[i] = e_p[i]->deg;
  }
  
  for (int f = 0; f < faces_p; f++){
    int volume_nodes_p = d4est_lgl_get_nodes((P4EST_DIM), e_p[f]->deg);
    int face_nodes_p = d4est_lgl_get_nodes((P4EST_DIM)-1, e_p[f]->deg);
    if (e_p_is_ghost[f] == 0){

      double* temp_b4_proj_VT = (p_side_matches_mortar) ? VT_w_term1_mortar_lobatto : proj_VT_w_term1_mortar_lobatto_p;
      double* temp_proj_VT = d4est_mortars_with_opt_reorient_if_needed
                             (
                              d4est_ops,
                              e_p,
                              e_m,
                              orientation,
                              f_p,
                              f_m,
                              temp_b4_proj_VT,
                              deg_p_lobatto_porder,
                              deg_p_lobatto,
                              total_side_nodes_p_lobatto,
                              faces_p
                             );
      
      d4est_operators_apply_lift(
                                 d4est_ops,
                                 (temp_proj_VT == NULL) ? &temp_b4_proj_VT[stride] : &temp_proj_VT[stride],
                                 (P4EST_DIM),
                                 e_p[f]->deg,
                                 f_p,
                                 &lifted_proj_VT_w_term1_mortar_lobatto_p[volume_stride]);

      if(temp_proj_VT != NULL)
        P4EST_FREE(temp_proj_VT);
      for (int d = 0; d < (P4EST_DIM); d++){

        /* term 2 does not need to be reoriented because it was reoriented then calculated */


        double* temp_term2 = (is_it_tree_boundary) ? VT_w_term2_mortar_lobatto_porder[d] : VT_w_term2_mortar_lobatto[d];
        
        d4est_operators_apply_lift(
                                   d4est_ops,
                                   (p_side_matches_mortar) ? &temp_term2[stride] : &proj_VT_w_term2_mortar_lobatto_p[d][stride],
                                   (P4EST_DIM),
                                   e_p[f]->deg,
                                   f_p,
                                   &lifted_proj_VT_w_term2_mortar_lobatto_p[d][volume_stride]);

        d4est_operators_apply_dij_transpose(d4est_ops,
                                            &lifted_proj_VT_w_term2_mortar_lobatto_p[d][volume_stride],
                                            (P4EST_DIM),
                                            e_p[f]->deg,
                                            d,
                                            &DT_lifted_proj_VT_w_term2_mortar_lobatto_p[d][volume_stride]
                                           );
        
        if (faces_p != faces_mortar){
            d4est_linalg_vec_scale(.5, &DT_lifted_proj_VT_w_term2_mortar_lobatto_p[d][volume_stride], volume_nodes_p);
        }
      }
      
      temp_b4_proj_VT = (p_side_matches_mortar) ? VT_w_term3_mortar_lobatto :  proj_VT_w_term3_mortar_lobatto_p;
      temp_proj_VT = d4est_mortars_with_opt_reorient_if_needed
                     (
                      d4est_ops,
                      e_p,
                      e_m,
                      orientation,
                      f_p,
                      f_m,
                      temp_b4_proj_VT,
                      deg_p_lobatto_porder,
                      deg_p_lobatto,
                      total_side_nodes_p_lobatto,
                      faces_p
                     );
      
      d4est_operators_apply_lift(
                                 d4est_ops,
                                 (temp_proj_VT == NULL) ? &temp_b4_proj_VT[stride] : &temp_proj_VT[stride],
                                 (P4EST_DIM),
                                 e_p[f]->deg,
                                 f_p,
                                 &lifted_proj_VT_w_term3_mortar_lobatto_p[volume_stride]);

      if (temp_proj_VT != NULL)
        P4EST_FREE(temp_proj_VT);
    }

    volume_stride += volume_nodes_p;
    stride += face_nodes_p;
  }



  if (is_it_tree_boundary){
    for (int d = 0; d < (P4EST_DIM); d++){
      P4EST_FREE(VT_w_term2_mortar_lobatto_porder[d]);
    }
  }
  D4EST_FREE_DIM_VEC(sj_n_du);

  
  P4EST_FREE(ones_mortar_quad);
  P4EST_FREE(sigma);
  
}

static void
d4est_laplacian_with_opt_flux_sipg_interface
(
 p4est_t* p4est,
 d4est_element_data_t** e_m,
 int faces_m,
 int f_m,
 int mortar_side_id_m,
 int* e_m_is_ghost,
 d4est_element_data_t** e_p,
 int faces_p,
 int f_p,
 int mortar_side_id_p,
 int* e_p_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_laplacian_with_opt_flux_interface_data_t* mortar_data,
 void* params
)
{
  d4est_laplacian_with_opt_flux_sipg_params_t* ip_flux_params = (d4est_laplacian_with_opt_flux_sipg_params_t*) params;


  /* if flux is symmetric don't recompute p-side data */
  if(
     mortar_side_id_m != 0
     &&
     mortar_side_id_p == ip_flux_params->last_mortar_side_id_m){
    return;
  }
  
  ip_flux_params->last_mortar_side_id_m = mortar_side_id_m;

  
  int faces_mortar = mortar_data->faces_mortar;
  int total_side_nodes_m_lobatto = mortar_data->total_side_nodes_m_lobatto;
  int total_side_nodes_p_lobatto = mortar_data->total_side_nodes_p_lobatto;
  int total_nodes_mortar_lobatto = mortar_data->total_nodes_mortar_lobatto;
  int total_nodes_mortar_quad = mortar_data->total_nodes_mortar_quad;
  
  int total_volume_nodes_m_lobatto = 0;
  for (int i = 0; i < faces_m; i++){
    total_volume_nodes_m_lobatto += d4est_lgl_get_nodes((P4EST_DIM), e_m[i]->deg);
  }

  int total_volume_nodes_p_lobatto = 0;
  for (int i = 0; i < faces_p; i++){
    total_volume_nodes_p_lobatto += d4est_lgl_get_nodes((P4EST_DIM), e_p[i]->deg);
  }

  
  double* lifted_proj_VT_w_term1_mortar_lobatto_m = P4EST_ALLOC(double, total_volume_nodes_m_lobatto);
  double* proj_VT_w_term1_mortar_lobatto_m = P4EST_ALLOC(double, total_side_nodes_m_lobatto);  
  double* DT_lifted_proj_VT_w_term2_mortar_lobatto_m [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(DT_lifted_proj_VT_w_term2_mortar_lobatto_m, total_volume_nodes_m_lobatto);
  double* lifted_proj_VT_w_term2_mortar_lobatto_m [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(lifted_proj_VT_w_term2_mortar_lobatto_m, total_volume_nodes_m_lobatto);
  double* proj_VT_w_term2_mortar_lobatto_m [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(proj_VT_w_term2_mortar_lobatto_m, total_side_nodes_m_lobatto);
  double* lifted_proj_VT_w_term3_mortar_lobatto_m = P4EST_ALLOC(double, total_volume_nodes_m_lobatto);
  double* proj_VT_w_term3_mortar_lobatto_m = P4EST_ALLOC(double, total_side_nodes_m_lobatto);


  double* lifted_proj_VT_w_term1_mortar_lobatto_p = P4EST_ALLOC(double, total_volume_nodes_p_lobatto);
  double* proj_VT_w_term1_mortar_lobatto_p = P4EST_ALLOC(double, total_side_nodes_p_lobatto);  
  double* DT_lifted_proj_VT_w_term2_mortar_lobatto_p [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(DT_lifted_proj_VT_w_term2_mortar_lobatto_p, total_volume_nodes_p_lobatto);
  double* lifted_proj_VT_w_term2_mortar_lobatto_p [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(lifted_proj_VT_w_term2_mortar_lobatto_p, total_volume_nodes_p_lobatto);
  double* proj_VT_w_term2_mortar_lobatto_p [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(proj_VT_w_term2_mortar_lobatto_p, total_side_nodes_p_lobatto);
  double* lifted_proj_VT_w_term3_mortar_lobatto_p = P4EST_ALLOC(double, total_volume_nodes_p_lobatto);
  double* proj_VT_w_term3_mortar_lobatto_p = P4EST_ALLOC(double, total_side_nodes_p_lobatto);
  
  double* VT_w_term1_mortar_lobatto = P4EST_ALLOC(double, total_nodes_mortar_lobatto);
  double* VT_w_term2_mortar_lobatto [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(VT_w_term2_mortar_lobatto, total_nodes_mortar_lobatto);
  double* VT_w_term3_mortar_lobatto= P4EST_ALLOC(double, total_nodes_mortar_lobatto);
  
  double* term1_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad); 
  double* term2_mortar_quad [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(term2_mortar_quad, total_nodes_mortar_quad);
  double* term3_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  
  d4est_laplacian_with_opt_flux_sipg_interface_aux
    (
     p4est,
     e_m,
     faces_m,
     f_m,
     mortar_side_id_m,
     e_m_is_ghost,
     e_p,
     faces_p,
     f_p,
     mortar_side_id_p,
     e_p_is_ghost,
     orientation,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     mortar_data,
     params,
     lifted_proj_VT_w_term1_mortar_lobatto_m,
     proj_VT_w_term1_mortar_lobatto_m,
     lifted_proj_VT_w_term1_mortar_lobatto_p,
     proj_VT_w_term1_mortar_lobatto_p,
     VT_w_term1_mortar_lobatto,
     term1_mortar_quad,
     DT_lifted_proj_VT_w_term2_mortar_lobatto_m,
     lifted_proj_VT_w_term2_mortar_lobatto_m,
     proj_VT_w_term2_mortar_lobatto_m,
     DT_lifted_proj_VT_w_term2_mortar_lobatto_p,
     lifted_proj_VT_w_term2_mortar_lobatto_p,
     proj_VT_w_term2_mortar_lobatto_p,
     VT_w_term2_mortar_lobatto,
     term2_mortar_quad,
     lifted_proj_VT_w_term3_mortar_lobatto_m,
     proj_VT_w_term3_mortar_lobatto_m,
     lifted_proj_VT_w_term3_mortar_lobatto_p,
     proj_VT_w_term3_mortar_lobatto_p,
     VT_w_term3_mortar_lobatto,
     term3_mortar_quad    
    );
  
  
  int stride = 0;
  for (int f = 0; f < faces_m; f++){
    int volume_nodes_m = d4est_lgl_get_nodes((P4EST_DIM), e_m[f]->deg);
    /* printf("Element = %d\n, face = %d\n", e_m[f]->id, f_m); */
    if (e_m_is_ghost[f] == 0){    
      
      for (int i = 0; i < volume_nodes_m; i++){
        for (int d = 0; d < (P4EST_DIM); d++){
          mortar_data->Au_m[f][i] +=  DT_lifted_proj_VT_w_term2_mortar_lobatto_m[d][i + stride];
        }

        mortar_data->Au_m[f][i] += lifted_proj_VT_w_term3_mortar_lobatto_m[i + stride];

        mortar_data->Au_m[f][i] += lifted_proj_VT_w_term1_mortar_lobatto_m[i + stride];
        /* if(e_m[0]->tree != e_p[0]->tree) */
        /* printf("mortar_data->Au_m[%d][%d] = %.15f\n", f, i, mortar_data->Au_m[f][i]); */
      }
    }
    stride += volume_nodes_m;
  }


  stride = 0;
  for (int f = 0; f < faces_p; f++){
    int volume_nodes_p = d4est_lgl_get_nodes((P4EST_DIM), e_p[f]->deg);
    if (e_p_is_ghost[f] == 0){      
      /* printf("Element = %d\n, face = %d\n", e_p[f]->id, f_p); */
      for (int i = 0; i < volume_nodes_p; i++){
        /* printf("term 2\n"); */
        for (int d = 0; d < (P4EST_DIM); d++){
          mortar_data->Au_p[f][i] +=  DT_lifted_proj_VT_w_term2_mortar_lobatto_p[d][i + stride];
        }

        /* printf("term 3\n"); */
        /* need negative here because normal sign switches */
        mortar_data->Au_p[f][i] -= lifted_proj_VT_w_term3_mortar_lobatto_p[i + stride];

        /* printf("term 1\n"); */
        /* need negative here because normal sign switches */
        mortar_data->Au_p[f][i] -= lifted_proj_VT_w_term1_mortar_lobatto_p[i + stride];
        /* if(e_m[0]->tree != e_p[0]->tree) */
        /* printf("mortar_data->Au_p[%d][%d] = %.15f\n", f, i, mortar_data->Au_p[f][i]); */
      }
    }
    stride += volume_nodes_p;
  }
  
  
  P4EST_FREE(lifted_proj_VT_w_term1_mortar_lobatto_m);
  P4EST_FREE(proj_VT_w_term1_mortar_lobatto_m);
  D4EST_FREE_DIM_VEC(DT_lifted_proj_VT_w_term2_mortar_lobatto_m);
  D4EST_FREE_DIM_VEC(lifted_proj_VT_w_term2_mortar_lobatto_m);
  D4EST_FREE_DIM_VEC(proj_VT_w_term2_mortar_lobatto_m);
  P4EST_FREE(lifted_proj_VT_w_term3_mortar_lobatto_m);
  P4EST_FREE(proj_VT_w_term3_mortar_lobatto_m);

  P4EST_FREE(lifted_proj_VT_w_term1_mortar_lobatto_p);
  P4EST_FREE(proj_VT_w_term1_mortar_lobatto_p);
  D4EST_FREE_DIM_VEC(DT_lifted_proj_VT_w_term2_mortar_lobatto_p);
  D4EST_FREE_DIM_VEC(lifted_proj_VT_w_term2_mortar_lobatto_p);
  D4EST_FREE_DIM_VEC(proj_VT_w_term2_mortar_lobatto_p);
  P4EST_FREE(lifted_proj_VT_w_term3_mortar_lobatto_p);
  P4EST_FREE(proj_VT_w_term3_mortar_lobatto_p);
  
  P4EST_FREE(VT_w_term1_mortar_lobatto);
  D4EST_FREE_DIM_VEC(VT_w_term2_mortar_lobatto);
  P4EST_FREE(VT_w_term3_mortar_lobatto);

  P4EST_FREE(term1_mortar_quad);
  D4EST_FREE_DIM_VEC(term2_mortar_quad);
  P4EST_FREE(term3_mortar_quad);  
}


static double
d4est_laplacian_with_opt_flux_sipg_penalty_meanp_sqr_over_meanh
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
d4est_laplacian_with_opt_flux_sipg_penalty_mean_p_sqr_over_h
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
d4est_laplacian_with_opt_flux_sipg_penalty_maxp_sqr_over_minh
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
d4est_laplacian_with_opt_flux_sipg_params_get_string_from_penalty_fcn
(
 penalty_calc_t fcn,
 char* string
)
{
  if (fcn == d4est_laplacian_with_opt_flux_sipg_penalty_maxp_sqr_over_minh){
    strcpy(string,"maxp_sqr_over_minh");
  }
  else if (fcn == d4est_laplacian_with_opt_flux_sipg_penalty_meanp_sqr_over_meanh){
    strcpy(string,"meanp_sqr_over_meanh");
  }
  else if (fcn == d4est_laplacian_with_opt_flux_sipg_penalty_mean_p_sqr_over_h){
    strcpy(string,"mean_p_sqr_over_h");
  }
  else {
    strcpy(string,"not_set");
  }
}

static penalty_calc_t
d4est_laplacian_with_opt_flux_sipg_get_penalty_fcn_from_string
(
 const char* string
)
{
  if (d4est_util_match(string,"maxp_sqr_over_minh")){
    return d4est_laplacian_with_opt_flux_sipg_penalty_maxp_sqr_over_minh;
  }
  else if (d4est_util_match(string,"meanp_sqr_over_meanh")){
    return d4est_laplacian_with_opt_flux_sipg_penalty_meanp_sqr_over_meanh;
  }
  else if (d4est_util_match(string,"mean_p_sqr_over_h")){
    return d4est_laplacian_with_opt_flux_sipg_penalty_mean_p_sqr_over_h;
  }
  else {
    D4EST_ABORT("This ip flux penalty calculation fcn does not exist");
    return NULL;
  }
}


static
int d4est_laplacian_with_opt_flux_sipg_region_multiplier_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_laplacian_with_opt_flux_sipg_params_t* pconfig = (d4est_laplacian_with_opt_flux_sipg_params_t*)user;

  int hit = 0;
  for (int i = 0; i < pconfig->sipg_number_of_regions; i++){
    char* mult_name;
    asprintf(&mult_name,"sipg_region%d_multiplier", i);

    if (d4est_util_match_couple(section,"flux",name,mult_name)) {
      D4EST_ASSERT(pconfig->sipg_region_multipliers[i] == -1);
      pconfig->sipg_region_multipliers[i] = atof(value);
      hit++;
    }
    free(mult_name);
  }
  if (hit)
    return 1;
  else
    return 0;
}

static
int d4est_laplacian_with_opt_flux_sipg_region_boundary_multiplier_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_laplacian_with_opt_flux_sipg_params_t* pconfig = (d4est_laplacian_with_opt_flux_sipg_params_t*)user;

  int hit = 0;

  for (int i = 0; i < pconfig->sipg_number_of_region_boundaries; i++){
    char* mult_name;
    asprintf(&mult_name,"sipg_region_boundary%d_multiplier", i);

    if (d4est_util_match_couple(section,"flux",name,mult_name)) {
      D4EST_ASSERT(pconfig->sipg_region_boundary_multipliers[i] == -1);
      pconfig->sipg_region_boundary_multipliers[i] = atof(value);
      hit++;
    }
    free(mult_name);
  }
  
  if (hit)
    return 1;
  else
    return 0;
}


static
int d4est_laplacian_with_opt_flux_sipg_params_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_laplacian_with_opt_flux_sipg_params_t* pconfig = (d4est_laplacian_with_opt_flux_sipg_params_t*)user;

  if (d4est_util_match_couple(section,"flux",name,"sipg_penalty_prefactor")) {
    D4EST_ASSERT(pconfig->sipg_penalty_prefactor == -1);
    pconfig->sipg_penalty_prefactor = atof(value);
  }
  else if (d4est_util_match_couple(section,"flux",name,"sipg_penalty_fcn")) {
    D4EST_ASSERT(pconfig->sipg_penalty_fcn == NULL);
    pconfig->sipg_penalty_fcn = d4est_laplacian_with_opt_flux_sipg_get_penalty_fcn_from_string(value);
  }
  else if (d4est_util_match_couple(section,"flux",name,"sipg_use_region_multipliers")) {
    D4EST_ASSERT(pconfig->sipg_use_region_multipliers == 0);
    pconfig->sipg_use_region_multipliers = atoi(value);
  }
  else if (d4est_util_match_couple(section,"flux",name,"sipg_use_region_boundary_multipliers")) {
    D4EST_ASSERT(pconfig->sipg_use_region_boundary_multipliers == 0);
    pconfig->sipg_use_region_boundary_multipliers = atoi(value);
  }
  else if (d4est_util_match_couple(section,"flux",name,"sipg_number_of_regions")) {
    D4EST_ASSERT(pconfig->sipg_number_of_regions == -1);
    pconfig->sipg_number_of_regions = atoi(value);
  } 
  else if (d4est_util_match_couple(section,"flux",name,"sipg_number_of_region_boundaries")) {
    D4EST_ASSERT(pconfig->sipg_number_of_region_boundaries == -1);
    pconfig->sipg_number_of_region_boundaries = atoi(value);
  }  
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

double
d4est_laplacian_with_opt_flux_sipg_get_multiplier
(
 d4est_laplacian_with_opt_flux_sipg_params_t* input,
 int region_m,
 int region_p
){
  if (!input->sipg_use_region_boundary_multipliers &&
      !input->sipg_use_region_multipliers){
    return 1.0;
  }

  double multiplier = 1.0;
  if (region_m == region_p){
    if (input->sipg_use_region_multipliers){
      multiplier = input->sipg_region_multipliers[region_m];
    }
  }
  else {
    int smaller_region = (region_m < region_p) ? region_m : region_p;
    if (input->sipg_use_region_boundary_multipliers){
      multiplier = input->sipg_region_boundary_multipliers[smaller_region];
    }
    else if (input->sipg_use_region_multipliers){
      double mult_m = input->sipg_region_multipliers[region_m];
      double mult_p = input->sipg_region_multipliers[region_p];
      multiplier = .5*(mult_m + mult_p);
    }
    else {
      multiplier = 1.0;
    }
  }

  return multiplier;
}

void
d4est_laplacian_with_opt_flux_sipg_params_input
(
 p4est_t* p4est,
 const char* printf_prefix,
 const char* input_file,
 d4est_laplacian_with_opt_flux_sipg_params_t* input
)
{
  input->sipg_penalty_fcn = NULL;
  input->sipg_penalty_prefactor = -1.;
  input->sipg_use_region_multipliers = 0;
  input->sipg_use_region_boundary_multipliers = 0;
  input->sipg_number_of_regions = -1;
  input->sipg_number_of_region_boundaries = -1;

  if (ini_parse(input_file, d4est_laplacian_with_opt_flux_sipg_params_input_handler, input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("flux", input->sipg_penalty_prefactor, -1);
  D4EST_CHECK_INPUT("flux", input->sipg_penalty_fcn, NULL);

  if(input->sipg_use_region_multipliers){
    
    D4EST_CHECK_INPUT("flux", input->sipg_number_of_regions, -1);

    if (input->sipg_number_of_regions > 10){
      D4EST_ABORT("sipg_number_of_regions is capped at 10 in the code for now");
    }

    for (int i = 0; i < 10; i++){
      input->sipg_region_multipliers[i] = -1;
    }
    
    if (ini_parse(
                  input_file,
                  d4est_laplacian_with_opt_flux_sipg_region_multiplier_input_handler,
                  input) < 0) {
      D4EST_ABORT("Can't load input file");
    }
    
    for (int i = 0; i < input->sipg_number_of_regions; i++){
      D4EST_CHECK_INPUT("flux", input->sipg_region_multipliers[i], -1);
    }
  }

  if(input->sipg_use_region_boundary_multipliers){
    
    D4EST_CHECK_INPUT("flux", input->sipg_number_of_region_boundaries, -1);

    if (input->sipg_number_of_region_boundaries > 10){
      D4EST_ABORT("sipg number_of_region_boundaries capped at 10 now");
    }
    
    for (int i = 0; i < 10; i++){
      input->sipg_region_boundary_multipliers[i] = -1;
    }
    
    if (ini_parse(
                  input_file,
                  d4est_laplacian_with_opt_flux_sipg_region_boundary_multiplier_input_handler,
                  input) < 0) {
      D4EST_ABORT("Can't load input file");
    }
    
    for (int i = 0; i < input->sipg_number_of_region_boundaries; i++){
      D4EST_CHECK_INPUT("flux", input->sipg_region_boundary_multipliers[i], -1);
    }
  }
  
  char penalty_calculate_fcn [50];
  char h_eq [50];

  d4est_laplacian_with_opt_flux_sipg_params_get_string_from_penalty_fcn (input->sipg_penalty_fcn,penalty_calculate_fcn);

  double check_function = input->sipg_penalty_fcn(1,.1,1,.1,0.);
  D4EST_ASSERT(check_function == 0.);
  
  if(p4est->mpirank == 0){
    zlog_category_t* c_def = zlog_get_category("d4est_laplacian_with_opt_flux_sipg_params");
    
    zlog_info(c_def,"sipg_params_penalty_prefactor = %f\n", input->sipg_penalty_prefactor);
    zlog_info(c_def,"sipg_params_penalty_calculate_fcn = %s\n", penalty_calculate_fcn);
    for (int i = 0; i < input->sipg_number_of_regions; i++){
      for (int j = 0; j < input->sipg_number_of_regions; j++){
        zlog_info(c_def,"multiplier for face between regions %d, %d = %f\n", i,j, d4est_laplacian_with_opt_flux_sipg_get_multiplier(input,i,j));
      }
    }
  }
}

void
d4est_laplacian_with_opt_flux_sipg_params_new
(
 p4est_t* p4est,
 const char* print_prefix,
 const char* input_file,
 d4est_laplacian_with_opt_flux_data_t* d4est_laplacian_with_opt_flux_data
)
{
  d4est_laplacian_with_opt_flux_sipg_params_t* d4est_laplacian_with_opt_flux_sipg_params = P4EST_ALLOC(d4est_laplacian_with_opt_flux_sipg_params_t, 1); 
  d4est_laplacian_with_opt_flux_sipg_params_input(p4est, print_prefix, input_file, d4est_laplacian_with_opt_flux_sipg_params);

  d4est_laplacian_with_opt_flux_sipg_params->size_params = NULL;
  
  d4est_laplacian_with_opt_flux_data->flux_data = d4est_laplacian_with_opt_flux_sipg_params;
  d4est_laplacian_with_opt_flux_data->interface_fcn = d4est_laplacian_with_opt_flux_sipg_interface;
  d4est_laplacian_with_opt_flux_data->skip_p_side = 1;
  d4est_laplacian_with_opt_flux_data->last_mortar_side_id_m = 1;
  
  if (d4est_laplacian_with_opt_flux_data->bc_type == BC_DIRICHLET){
    d4est_laplacian_with_opt_flux_data->boundary_fcn = d4est_laplacian_with_opt_flux_sipg_dirichlet;
  }
  else if (d4est_laplacian_with_opt_flux_data->bc_type == BC_ROBIN){
    d4est_laplacian_with_opt_flux_data->boundary_fcn = d4est_laplacian_with_opt_flux_sipg_robin;
  }
  else {
    D4EST_ABORT("Not a supported boundary condition");
  }

  d4est_laplacian_with_opt_flux_data->destroy = d4est_laplacian_with_opt_flux_sipg_params_destroy;
}

void
d4est_laplacian_with_opt_flux_sipg_params_destroy
(
 d4est_laplacian_with_opt_flux_data_t* data
){
  P4EST_FREE(data->flux_data);
  P4EST_FREE(data);
}

