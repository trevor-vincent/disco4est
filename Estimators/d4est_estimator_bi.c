#include <pXest.h>
#include <util.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_element_data.h>
#include <d4est_linalg.h>
#include <d4est_mesh.h>
#include <d4est_estimator_bi.h>
#include <d4est_mortars.h>

void
d4est_estimator_bi_compute
(
 p4est_t* p4est,
 d4est_elliptic_problem_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 d4est_estimator_bi_penalty_data_t bi_penalty_data,
 d4est_grid_fcn_t u_bndry_fcn,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 double (*get_diam)(d4est_element_data_t*)
)
{
  d4est_elliptic_eqns_build_residual
    (
     p4est,
     ghost,
     fcns,
     vecs,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad
    );
  
  d4est_mesh_compute_l2_norm_sqr
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     vecs->Au,
     vecs->local_nodes,
     STORE_LOCALLY
    );

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
        double* eta2 = &(ed->local_estimator);
        /* handle ||R||^2 * h^2/p^2 term */
        double h = get_diam(ed);
        *eta2 *= h*h/(deg*deg);


        d4est_linalg_copy_1st_to_2nd
          (
           &(vecs->u[ed->nodal_stride]),
           &(ed->u_elem)[0],
           volume_nodes_lobatto
          );
    
        for (int i = 0; i < (P4EST_DIM); i++){
          d4est_operators_apply_dij(d4est_ops, &(vecs->u[ed->nodal_stride]), (P4EST_DIM), ed->deg, i, &ed->dudr_elem[i][0]);
        }

        
      }

    }

  d4est_mortar_fcn_ptrs_t d4est_estimator_bi_flux_fcn_ptrs
    = d4est_estimator_bi_dirichlet_fetch_fcns
    (
     u_bndry_fcn,
     &bi_penalty_data
    );
  
  d4est_mortar_compute_flux_on_local_elements
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &d4est_estimator_bi_flux_fcn_ptrs,
     EXCHANGE_GHOST_DATA
    );
}



static void
d4est_estimator_bi_dirichlet
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
  d4est_grid_fcn_t u_at_bndry = bndry_fcn;
  d4est_estimator_bi_penalty_data_t* penalty_data = params;
  int face_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_quad = d4est_lgl_get_nodes((P4EST_DIM) - 1, e_m->deg_quad);

  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* u_on_f_m_min_u_at_bndry_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* u_on_f_m_min_u_at_bndry_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  

  double* MJe2 = P4EST_ALLOC(double, face_nodes_m_quad);
  double* sj_on_f_m_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  double* j_div_sj_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  double* Je2 = P4EST_ALLOC(double, face_nodes_m_quad);

  double* xyz_on_f_m [(P4EST_DIM)];
  double* n_on_f_m_quad [(P4EST_DIM)];
  /* double* sj_n_on_f_m_quad [(P4EST_DIM)]; */

  
  for (int d = 0; d < (P4EST_DIM); d++) {
    xyz_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m_lobatto);
    n_on_f_m_quad[d] = P4EST_ALLOC(double, face_nodes_m_quad);


    d4est_operators_apply_slicer(d4est_ops,
                        e_m->xyz[d],
                        (P4EST_DIM),
                        f_m,
                        e_m->deg,
                        xyz_on_f_m[d]);
    
  }

  
  
  d4est_operators_apply_slicer(d4est_ops,
                      e_m->u_elem,
                      (P4EST_DIM),
                      f_m,
                      e_m->deg,
                      u_m_on_f_m);


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
     NULL,
     sj_on_f_m_quad,
     n_on_f_m_quad,
     NULL,
     j_div_sj_quad,
     COMPUTE_NORMAL_USING_JACOBIAN
    );
  
  double* Je2_prefactor = P4EST_ALLOC(double, face_nodes_m_quad);
  double h, h_min;
  if (penalty_data->sipg_flux_h == H_EQ_J_DIV_SJ_MIN){
    h_min = util_min_dbl_array(j_div_sj_quad, face_nodes_m_quad);
  }
    
  for (int i = 0; i < face_nodes_m_quad; i++){
    if (penalty_data->sipg_flux_h == H_EQ_J_DIV_SJ){
      h = j_div_sj_quad[i];
    }
    else if (penalty_data->sipg_flux_h == H_EQ_J_DIV_SJ_MIN){
      h = h_min;
    }
    else {
      mpi_abort("[D4EST_ERROR]: NOT SUPPORTED IP_FLUX_H_CALC");
    }
    Je2_prefactor[i] = penalty_data->u_dirichlet_penalty_fcn
                       (
                        e_m->deg,
                        h,
                        e_m->deg,
                        h,
                        penalty_data->penalty_prefactor
                       );
  }


  for (int i = 0; i < face_nodes_m_lobatto; i++){
    u_on_f_m_min_u_at_bndry_lobatto[i] = u_m_on_f_m[i]
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


  d4est_quadrature_interpolate
    (
     d4est_ops,
     d4est_quad,
     d4est_geom,
     &face_object,
     QUAD_OBJECT_MORTAR,
     QUAD_INTEGRAND_UNKNOWN,
     u_on_f_m_min_u_at_bndry_lobatto,
     e_m->deg,
     u_on_f_m_min_u_at_bndry_quad,
     e_m->deg_quad
    );

    
  
  for (int dim = 0; dim < (P4EST_DIM); dim++){

    for(int i = 0; i < face_nodes_m_quad; i++){
      Je2[i] = n_on_f_m_quad[dim][i]*Je2_prefactor[i]*(u_on_f_m_min_u_at_bndry_quad[i]);
    }

    double Je2MJe2 = d4est_quadrature_innerproduct
                     (
                      d4est_ops,
                      d4est_geom,
                      d4est_quad,
                      &face_object,
                      QUAD_OBJECT_MORTAR,
                      QUAD_INTEGRAND_UNKNOWN,
                      Je2,
                      Je2,
                      sj_on_f_m_quad,
                      e_m->deg_quad
                     );


    
    e_m->local_estimator += Je2MJe2;

  }

  P4EST_FREE(Je2_prefactor);
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(sj_on_f_m_quad);
  P4EST_FREE(j_div_sj_quad);
  P4EST_FREE(Je2);
  P4EST_FREE(MJe2);
  P4EST_FREE(u_on_f_m_min_u_at_bndry_lobatto);
  P4EST_FREE(u_on_f_m_min_u_at_bndry_quad);
  
  for (int d = 0; d < (P4EST_DIM); d++){
    P4EST_FREE(xyz_on_f_m[d]);
    P4EST_FREE(n_on_f_m_quad[d]);
  }
}

static void
d4est_estimator_bi_interface
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
  int stride;
  int deg_p_lobatto [(P4EST_HALF)];
  int face_nodes_p_lobatto [(P4EST_HALF)];
  int face_nodes_p_quad [(P4EST_HALF)];
  int deg_m_lobatto [(P4EST_HALF)];
  int face_nodes_m_lobatto [(P4EST_HALF)];
  int face_nodes_m_quad [(P4EST_HALF)];
  int nodes_mortar_quad [(P4EST_HALF)];
  int nodes_mortar_lobatto [(P4EST_HALF)];
  int deg_mortar_quad [(P4EST_HALF)];
  int deg_mortar_lobatto [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
  /* double Je1_prefactor_mortar [(P4EST_HALF)]; */
  /* double Je2_prefactor_mortar [(P4EST_HALF)]; */
  int deg_p_lobatto_porder [(P4EST_HALF)];
  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);
  d4est_estimator_bi_penalty_data_t* penalty_data = params;
  
  
  int total_side_nodes_m_lobatto = 0;
  int total_side_nodes_m_quad = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_lobatto[i] = e_m[i]->deg;
    /* deg_m_quad[i] = e_m[i]->deg_quad; */
    
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
  double* sj_on_f_m_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* n_on_f_m_mortar_quad [(P4EST_DIM)];
  /* double* sj_n_on_f_m_mortar_quad [(P4EST_DIM)];   */


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
    n_on_f_m_mortar_quad[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    /* sj_n_on_f_m_mortar_quad[i] = P4EST_ALLOC(double, total_nodes_mortar_quad); */
  }
 
  double* Je1 = P4EST_ALLOC_ZERO(double, total_nodes_mortar_quad);
  double* MJe1 = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* Je2 [(P4EST_DIM)]; 
  for (int d = 0; d < (P4EST_DIM); d++) {
    Je2[d] = P4EST_ALLOC(double, total_nodes_mortar_quad);
  }
  double* MJe2 = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* tmp = P4EST_ALLOC(double, total_side_nodes_p_quad);

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

  P4EST_FREE(tmp);

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

    /* project (-)-u-derivative on the (-)-side faces and project q onto the (-)-side faces */
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
    
    /* compute the (+)-u-derivative 
   * and project on the (+)-side faces and project q onto the (+)-side faces */
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

    /* project the derivatives from (-) and (+) sides onto the mortar space */
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

  if(penalty_data->sipg_flux_h == H_EQ_J_DIV_SJ
     || penalty_data->sipg_flux_h == H_EQ_J_DIV_SJ_MIN
    ){
    j_div_sj_on_f_m_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
    j_div_sj_on_f_p_mortar_quad_porder =  P4EST_ALLOC(double, total_nodes_mortar_quad);
    j_div_sj_on_f_p_mortar_quad_porder_oriented =  P4EST_ALLOC(double, total_nodes_mortar_quad);
  }
  else {
    mpi_abort("[D4EST_ERROR]: ip_flux_params->sipg_flux_h can only be H_EQ_J_DIV_SJ(_MIN)\n");
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
     NULL,
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
  if (penalty_data->sipg_flux_h == H_EQ_J_DIV_SJ_MIN){
    hp_min = util_min_dbl_array(j_div_sj_on_f_p_mortar_quad_porder_oriented, total_nodes_mortar_quad);
    hm_min = util_min_dbl_array(j_div_sj_on_f_m_mortar_quad,total_nodes_mortar_quad);
  }
  double* Je1_prefactor = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* Je2_prefactor = P4EST_ALLOC(double, total_nodes_mortar_quad);
  
  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_quad[f]; k++){
      int ks = k + stride;
      int is_it_min = (penalty_data->sipg_flux_h == H_EQ_J_DIV_SJ_MIN);
      double hp = (is_it_min) ? hp_min : j_div_sj_on_f_p_mortar_quad_porder_oriented[ks];
      double hm = (is_it_min) ? hm_min : j_div_sj_on_f_m_mortar_quad[ks];
      
      Je1_prefactor[ks] = penalty_data->gradu_penalty_fcn
                           (
                            e_m[f]->deg,
                            hm,
                            e_p_oriented[f]->deg,
                            hp,
                            penalty_data->penalty_prefactor
                           ); 

      Je2_prefactor[ks] = penalty_data->u_penalty_fcn
                          (
                           e_m[f]->deg,
                           hm,
                           e_p_oriented[f]->deg,
                           hp,
                           penalty_data->penalty_prefactor
                          );    
      

    }
    stride += nodes_mortar_quad[f];
  }

  

  /* calculate symmetric interior penalty flux */
  int k;
  int f;
  int ks;
  double n_ks;
  double sj_ks;
  stride = 0;
  for (f = 0; f < faces_mortar; f++){
    for (k = 0; k < nodes_mortar_quad[f]; k++){
      ks = k + stride;
      sj_ks = sj_on_f_m_mortar_quad[ks];
      Je1[ks] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        n_ks = n_on_f_m_mortar_quad[d][ks];        
        Je1[ks] += Je1_prefactor[ks]*n_ks*
                   (dudx_m_on_f_m_mortar_quad[d][ks] - dudx_p_on_f_p_mortar_quad[d][ks]);
        Je2[d][ks] = n_ks*u_m_on_f_m_mortar_quad[ks];
        Je2[d][ks] -= n_ks*u_p_on_f_p_mortar_quad[ks];
        Je2[d][ks] *= Je2_prefactor[ks];
      }
    }
    stride += nodes_mortar_quad[f];
  }


    
  /* the contribution in every direction must be added up due to it being a vector norm */
  stride = 0;
  for (f = 0; f < faces_mortar; f++){
    for (int d = 0; d < (P4EST_DIM); d++){

      double Je2MJe2 = d4est_quadrature_innerproduct
                       (
                        d4est_ops,
                        d4est_geom,
                        d4est_quad,
                        &mortar_face_object_forder,
                        QUAD_OBJECT_MORTAR,
                        QUAD_INTEGRAND_UNKNOWN,
                        &Je2[d][stride],
                        &Je2[d][stride],
                        &sj_on_f_m_mortar_quad[stride],
                        deg_mortar_quad[f]
                       );

      
      /* even if it's a ghost we can still add it to the ghosts estimator because we will not send it back! */

      if(faces_m == (P4EST_HALF)){
        e_m[f]->local_estimator += Je2MJe2;
      }
      else{
        e_m[0]->local_estimator += Je2MJe2;
      }
    }
    stride += nodes_mortar_quad[f];
  }


    
  stride = 0;
  for (f = 0; f < faces_mortar; f++){  

    double Je1MJe1 = d4est_quadrature_innerproduct
                     (
                      d4est_ops,
                      d4est_geom,
                      d4est_quad,
                      &mortar_face_object_forder,
                      QUAD_OBJECT_MORTAR,
                      QUAD_INTEGRAND_UNKNOWN,
                      &Je1[stride],
                      &Je1[stride],
                      &sj_on_f_m_mortar_quad[stride],
                      deg_mortar_quad[f]
                     );
    
   
    /* even if it's a ghost we can still add it to the ghosts estimator because we will not send it back! */
    if(faces_m == (P4EST_HALF)){
      e_m[f]->local_estimator += Je1MJe1;
    }
    else{
      e_m[0]->local_estimator += Je1MJe1;
    }
    stride += nodes_mortar_quad[f];
  }


  D4EST_FREE_DBYD_MAT(drst_dxyz_m_on_mortar_quad);
  D4EST_FREE_DBYD_MAT(drst_dxyz_p_on_mortar_quad_porder);

  P4EST_FREE(Je1_prefactor);
  P4EST_FREE(Je2_prefactor);
  P4EST_FREE(j_div_sj_on_f_m_mortar_quad);
  P4EST_FREE(j_div_sj_on_f_p_mortar_quad_porder);
  P4EST_FREE(j_div_sj_on_f_p_mortar_quad_porder_oriented);
  
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_p_on_f_p);
  P4EST_FREE(u_m_on_f_m_mortar);
  P4EST_FREE(u_m_on_f_m_mortar_quad);
  P4EST_FREE(u_p_on_f_p_mortar);
  P4EST_FREE(u_p_on_f_p_mortar_quad);
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(Je2[i]);
    P4EST_FREE(dudr_p_on_f_p_porder[i]);
    P4EST_FREE(dudr_p_on_f_p_mortar_porder[i]);
    P4EST_FREE(dudr_p_on_f_p_mortar_quad_porder[i]);
    P4EST_FREE(dudx_p_on_f_p_mortar_quad_porder[i]);
    P4EST_FREE(dudx_p_on_f_p_mortar_quad[i]);
    P4EST_FREE(dudr_m_on_f_m[i]);
    P4EST_FREE(dudr_m_on_f_m_mortar[i]);
    P4EST_FREE(dudr_m_on_f_m_mortar_quad[i]);
    P4EST_FREE(dudx_m_on_f_m_mortar_quad[i]);
    P4EST_FREE(n_on_f_m_mortar_quad[i]);
  }
  P4EST_FREE(sj_on_f_m_mortar_quad);
  P4EST_FREE(Je1);
  P4EST_FREE(MJe1);
  P4EST_FREE(MJe2);

}



d4est_mortar_fcn_ptrs_t
d4est_estimator_bi_dirichlet_fetch_fcns
(
 d4est_grid_fcn_t bndry_fcn,
 d4est_estimator_bi_penalty_data_t* penalty_data
)
{
  d4est_mortar_fcn_ptrs_t d4est_estimator_bi_fcns;
  d4est_estimator_bi_fcns.flux_interface_fcn
    = d4est_estimator_bi_interface;

  d4est_estimator_bi_fcns.flux_boundary_fcn
    = d4est_estimator_bi_dirichlet;

  d4est_estimator_bi_fcns.bndry_fcn = bndry_fcn;
 
  d4est_estimator_bi_fcns.params = penalty_data;
  /* d4est_estimator_bi_u_prefactor_calculate_fcn = u_penalty_fcn; */
  /* d4est_estimator_bi_u_dirichlet_prefactor_calculate_fcn = u_dirichlet_penalty_fcn; */
  /* d4est_estimator_bi_gradu_prefactor_calculate_fcn = gradu_penalty_fcn; */
  
  return d4est_estimator_bi_fcns;
}
