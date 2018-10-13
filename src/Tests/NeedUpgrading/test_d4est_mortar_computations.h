#ifndef TEST_D4EST_MORTAR_COMPUTATIONS_H
#define TEST_D4EST_MORTAR_COMPUTATIONS_H 

typedef struct {

  double sj_and_n_interface_error;
  double ;
  
} testd4est_mortars_computations_data_t;

static void
testd4est_laplacian_flux_boundary_debug
(
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_xyz_fcn_t boundary_condition,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_laplacian_flux_boundary_data_t* boundary_data,
 d4est_laplacian_flux_sipg_debug_data_t* debug_data,
 void* params
 
)
{
  d4est_quadrature_mortar_t* face_object = boundary_data->face_object;
  int face_nodes_m_lobatto = boundary_data->face_nodes_m_lobatto;
  int face_nodes_m_quad = boundary_data->face_nodes_m_quad;
  double* u_m_on_f_m_quad = boundary_data->u_m_on_f_m_quad;
  double* sj_on_f_m_quad = boundary_data->sj_on_f_m_quad;
  double* j_div_sj_quad = boundary_data->j_div_sj_quad;
  double* xyz_on_f_m_lobatto [(P4EST_DIM)]; 
  double* drst_dxyz_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* dudx_m_on_f_m_quad [(P4EST_DIM)];  
  double* n_sj_on_f_m_quad [(P4EST_DIM)];
  D4EST_COPY_DBYD_MAT(boundary_data->drst_dxyz_quad, drst_dxyz_quad);
  D4EST_COPY_DIM_VEC(boundary_data->dudx_m_on_f_m_quad, dudx_m_on_f_m_quad);
  D4EST_COPY_DIM_VEC(boundary_data->n_sj_on_f_m_quad, n_sj_on_f_m_quad);
  D4EST_COPY_DIM_VEC(boundary_data->xyz_on_f_m_lobatto, xyz_on_f_m_lobatto);

 double* sigma = P4EST_ALLOC(double, face_nodes_m_quad);
  double h, h_min;
  if (ip_flux_params->sipg_flux_h ==H_EQ_J_DIV_SJ_MIN){
    h_min = d4est_util_min_dbl_array(j_div_sj_quad, face_nodes_m_quad);
  }
  for (int i = 0; i < face_nodes_m_quad; i++){
    int is_it_min = (ip_flux_params->sipg_flux_h == H_EQ_J_DIV_SJ_MIN);
    double h = (is_it_min) ? h_min : j_div_sj_quad[i];
    sigma[i] = sipg_kronbichler_flux_penalty_calculate_fcn
               (
                e_m->deg,
                h,
                e_m->deg,
                h,
                sipg_kronbichler_flux_penalty_prefactor
               );   
  }

  double* u_at_bndry_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* u_at_bndry_lobatto_to_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  

  for (int i = 0; i < face_nodes_m_lobatto; i++){
    u_at_bndry_lobatto[i] = boundary_condition
                            (
                             xyz_on_f_m_lobatto[0][i],
                             xyz_on_f_m_lobatto[1][i],
#if (P4EST_DIM)==3
                             xyz_on_f_m_lobatto[2][i],
#endif
                             ip_flux_params->user
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
     e_m->deg_quad
    );
  
  for(int i = 0; i < face_nodes_m_quad; i++){
    double u_m_on_f_m_min_u_at_bndry_quad
      = u_m_on_f_m_quad[i] - u_at_bndry_lobatto_to_quad[i];    
    term3_quad[i] = sj_on_f_m_quad[i]
                    *sigma[i]
                    *2.*u_m_on_f_m_min_u_at_bndry_quad;
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

  P4EST_FREE(u_at_bndry_lobatto);
  P4EST_FREE(u_at_bndry_lobatto_to_quad);
  P4EST_FREE(sigma);
  P4EST_FREE(term3_quad);  
}

static void
testd4est_flux_interface_debug
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
 d4est_laplacian_flux_interface_data_t* mortar_data,
 d4est_laplacian_flux_sipg_debug_data_t* debug_data,
 void* params
)
{
  d4est_quadrature_mortar_t* mortar_face_object = mortar_data->mortar_face_object;
  
  int faces_mortar = mortar_data->faces_mortar;
  int total_side_nodes_m_lobatto = mortar_data->total_side_nodes_m_lobatto;
  int total_nodes_mortar_lobatto = mortar_data->total_nodes_mortar_lobatto;
  int total_nodes_mortar_quad = mortar_data->total_nodes_mortar_quad;
  
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
  
  int* deg_mortar_quad = mortar_data->deg_mortar_quad;
  int* nodes_mortar_quad = mortar_data->nodes_mortar_quad;
  int* nodes_mortar_lobatto = mortar_data->nodes_mortar_lobatto;
  int* deg_mortar_lobatto = mortar_data->deg_mortar_lobatto;
  int* face_nodes_m_lobatto = mortar_data->deg_mortar_lobatto;
  int* deg_m_lobatto = mortar_data->deg_m_lobatto;
  int* deg_p_lobatto = mortar_data->deg_p_lobatto;

 
  double* sjvol_p_on_f_p_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* sjvol_p_on_f_p_mortar_quad_porder = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* nvol_p_on_f_p_mortar_quad [(P4EST_DIM)];
  double* nvol_p_on_f_p_mortar_quad_porder [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(nvol_p_on_f_p_mortar_quad,total_nodes_mortar_quad);
  D4EST_ALLOC_DIM_VEC(nvol_p_on_f_p_mortar_quad_porder,total_nodes_mortar_quad);

  testd4est_mortars_computations_data_t* data = params;
  
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
     NULL,
     sjvol_p_on_f_p_mortar_quad_porder,
     nvol_p_on_f_p_mortar_quad_porder,
     NULL,
     NULL,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

    d4est_operators_reorient_face_data
      (
       d4est_ops,
       &sjvol_p_on_f_p_mortar_quad_porder[oriented_face_mortar_stride],
       (P4EST_DIM)-1,
       deg_mortar_quad[face],
       orientation,
       f_m,
       f_p,
       &sjvol_p_on_f_p_mortar_quad[face_mortar_stride]
      );

    for (int d = 0; d < (P4EST_DIM); d++){
      d4est_operators_reorient_face_data
        (
         d4est_ops,
         &n_p_on_f_p_mortar_quad_porder[d][oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_quad[face],
         orientation,
         f_m,
         f_p,
         &n_p_on_f_p_mortar_quad[d][face_mortar_stride]
        );

      d4est_linalg_vec_scale(-1., &n_p_on_f_p_mortar_quad[d][face_mortar_stride], d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad[face]));      
    }

    double max_error = d4est_util_max_error(sj_p_on_f_p_mortar_quad, sj_m_on_f_m_mortar_quad, total_nodes_mortar_quad);

    for (int d = 0; d < (P4EST_DIM); d++){
      max_error += d4est_util_max_error(n_m_on_f_m_mortar_quad[d], n_p_on_f_p_mortar_quad[d], total_nodes_mortar_quad);
    }
    
    data->sj_and_n_interface_error += max_error;

  P4EST_FREE(sjvol_p_on_f_p_mortar_quad);
  P4EST_FREE(sjvol_p_on_f_p_mortar_quad_porder);
  D4EST_FREE_DIM_VEC(nvol_p_on_f_p_mortar_quad);
  D4EST_FREE_DIM_VEC(nvol_p_on_f_p_mortar_quad_porder);    
}

void
testd4est_mortars_computations_sj_and_n_on_interface
(
 double err_tol
)
{
  d4est_laplacian_flux_data_t* d4est_laplacian_flux_data = P4EST_ALLOC(d4est_laplacian_flux_data_t,1);
  testd4est_mortars_computations_data_t* data = P4EST_ALLOC(testd4est_mortars_computations_data_t, 1);
  data->sj_and_n_interface_error = 0.;
  
  d4est_laplacian_flux_data->user = data;
  d4est_laplacian_flux_data->interface_fcn = testd4est_mortars_computations_interface;
  d4est_laplacian_flux_data->boundary_fcn = NULL;
  d4est_laplacian_flux_data->boundary_condition = NULL;
  d4est_laplacian_flux_data->destroy = NULL;

  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data; 
        int deg = ed->deg;
        int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM),deg);
        ed->Au_elem = &(prob_vecs->Au[ed->nodal_stride]);
        d4est_util_fill_array(&(ed->u_elem[0]), 0., volume_nodes_lobatto);
        for (int i = 0; i < (P4EST_DIM); i++){
          d4est_util_fill_array(&ed->dudr_elem[i][0], 0., volume_nodes_lobatto);
        }
      }
    }
  
  d4est_mortars_fcn_ptrs_t flux_fcns = d4est_laplacian_flux_fetch_fcns(d4est_laplacian_flux_data);
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
  
  int err = (data->sj_and_n_interface_error < err_tol);
  P4EST_FREE(d4est_laplacian_flux_data);
  D4EST_ASSERT(err);
}




#endif
