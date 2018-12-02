static int field_size_of_ghost_fcn
(
 d4est_ghost_t* d4est_ghost,
 int gid,
 int tn,
 void* user_ctx
)
{
  return sizeof(d4est_mortar_side_data_t)*(P4EST_FACES);
}

static int field_size_of_mirror_fcn
(
 d4est_ghost_t* d4est_ghost,
 int mirror,
 int tn,
 void* user_ctx
)
{
 return sizeof(d4est_mortar_side_data_t)*(P4EST_FACES);
}

static int field_stride_of_mirror_fcn
(
 d4est_ghost_t* d4est_ghost,
 int mirr,
 int tn,
 void* user_ctx
)
{
 return sizeof(d4est_mortar_side_data_t)*(P4EST_FACES)*mirr;
}


d4est_solver_schwarz_mortar_data_t*
d4est_solver_schwarz_mortar_data_init
(
 p4est_t* p4est,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* d4est_ghost
)
{
  d4est_solver_schwarz_mortar_data_t* schwarz_mortar_data =
    P4EST_ALLOC(d4est_solver_schwarz_mortar_data_t, 1);
  
  d4est_mortar_side_data_t* mirror_data = P4EST_ALLOC(d4est_mortar_side_data_t, (P4EST_FACES)*num_mirror_elements);

  for (int i = 0; i < num_mirror_elements; i++){
    /* get mirror data */
    d4est_element_data_t* mirror_ed = d4est_ghost->mirror_elements[i];
    for (int f = 0; f < (P4EST_FACES); f++){
      int mortar_side_id = mirror_ed->face_belongs_to_which_mortar[f];
      mirror_data[i*(P4EST_FACES) + f] = d4est_factors->mortar_side_data[mortar_side_id];
    }
  }
  
 schwarz_mortar_data->mortar_side_ghost_data =
    d4est_ghost_data_ext_init
    (
     p4est,
     d4est_ghost,
     1,
     field_size_of_ghost_fcn,
     field_size_of_mirror_fcn,
     field_stride_of_mirror_fcn,
     NULL
    );

  d4est_ghost_data_ext_exchange
    (
     p4est,
     d4est_ghost,
     schwarz_mortar_data->mortar_side_ghost_data,
     (char**)&mirror_data
    );


  /* loop through data and set strides */
  
  int mortar_quad_ghost_size = 0;
  int boundary_quad_ghost_size = 0;
  
  int boundary_scalar_stride = 0;
  int mortar_scalar_stride = 0;

  for (int gid = 0; gid < num_ghost_elements; gid++){
    
    d4est_mortar_side_data_t* ghost_mortar_side_data =
      d4est_ghost_data_ext_get_field_on_element
      (
       &d4est_ghost->ghost_elements[gid],
       0,
       *d4est_ghost_data_ext
      );
    
    for (int f = 0; f < (P4EST_FACES); f++){
      int total_mortar_nodes_quad = 0;
 
      if (mortar_side_data->faces_p == 0){
        mortar_side_data->boundary_quad_stride = boundary_quad_stride;
        mortar_side_data->mortar_quad_stride = -1;
        total_mortar_nodes_quad += d4est_lgl_get_nodes((P4EST_DIM)-1, mortar_side_data->deg_quad_m[0]);
        boundary_quad_stride += total_mortar_nodes_quad;
      }
      else {
        for (int i = 0; i < faces_m; i++){
          for (int j = 0; j < faces_p; j++){
            /* find max degree for each face pair of the two sides*/
            
            int jnew = j;
            if (faces_p == (P4EST_HALF)){
              jnew = d4est_reference_reorient_face_order((P4EST_DIM)-1,
                                                         ghost_mortar_side_data->f_m,
                                                         ghost_mortar_side_data->f_p,
                                                         ghost_mortar_side_data->orientation, j);
            }
            deg_mortar_quad = d4est_util_max_int(ghost_mortar_side_data->deg_m_quad[i],
                                                 ghost_mortar_side_data->deg_p_quad[jnew]);
            mortar_nodes_quad = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_quad);
            total_mortar_nodes_quad += mortar_nodes_quad;
          }
        }
        mortar_side_data->mortar_quad_stride = mortar_quad_stride;
        mortar_side_data->boundary_quad_stride = -1;
      }
      
      mortar_quad_stride += total_mortar_nodes_quad;
      ghost_mortar_side_data++;
    }    
  }

  int mortar_quad_ghost_size = mortar_quad_stride;
  int boundary_quad_ghost_size = boundary_quad_stride;
   /* allocate mortar shit */
  schwarz_mortar_data->sj_m_mortar_quad = P4EST_ALLOC_ZERO(double, mortar_quad_ghost_size);
  schwarz_mortar_data->n_m_mortar_quad = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*mortar_quad_ghost_size);
  schwarz_mortar_data->hp_mortar_quad = P4EST_ALLOC_ZERO(double, mortar_quad_ghost_size);
  schwarz_mortar_data->hm_mortar_quad = P4EST_ALLOC_ZERO(double, mortar_quad_ghost_size);
  schwarz_mortar_data->drst_dxyz_m_mortar_quad = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*(P4EST_DIM)*mortar_quad_ghost_size);
  schwarz_mortar_data->drst_dxyz_p_mortar_quad_porder = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*(P4EST_DIM)*mortar_quad_ghost_size);
  schwarz_mortar_data->xyz_on_f_m_quad = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*boundary_quad_ghost_size);
  schwarz_mortar_data->xyz_on_f_m_lobatto = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*boundary_quad_ghost_size);

  for (int gid = 0; gid < num_ghost_elements; gid++){

    d4est_mortar_side_data_t* ghost_mortar_side_data =
      d4est_ghost_data_ext_get_field_on_element
      (
       &d4est_ghost->ghost_elements[gid],
       0,
       *d4est_ghost_data_ext
      );
    
    for (int f = 0; f < (P4EST_FACES); f++){

      int scalar_stride = ghost_mortar_side_data->mortar_quad_stride;
      int vector_stride = (P4EST_DIM)*ghost_mortar_side_data->boundary_quad_stride;
      int matrix_stride = (P4EST_DIM)*(P4EST_DIM)*ghost_mortar_side_data->boundary_quad_stride;
      int boundary_vector_stride = (P4EST_DIM)*ghost_mortar_side_data->boundary_quad_stride;
      
      double* xyz_on_f_m_quad [(P4EST_DIM)];
      double* xyz_on_f_m_lobatto [(P4EST_DIM)];
      double* n_m_mortar_quad [(P4EST_DIM)];
      double* drst_dxyz_m_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
      double* drst_dxyz_p_mortar_quad_porder [(P4EST_DIM)][(P4EST_DIM)];

      
      double* sj_m_mortar_quad = &schwarz_mortar_data->->sj_m_mortar_quad[scalar_stride];
      double* hm_mortar_quad = &schwarz_mortar_data->->hm_mortar_quad[scalar_stride];
      double* hp_mortar_quad = &schwarz_mortar_data->->hp_mortar_quad[scalar_stride];
  
      for (int i = 0; i < (P4EST_DIM); i++){
        n_m_mortar_quad[i] =
          &schwarz_mortar_data->->n_m_mortar_quad
          [vector_stride + i*total_mortar_nodes_quad];

        for (int j = 0; j < (P4EST_DIM); j++){
          drst_dxyz_m_mortar_quad[i][j] =
            &schwarz_mortar_data->->drst_dxyz_m_mortar_quad
            [matrix_stride + (i + j*(P4EST_DIM))*total_mortar_nodes_quad];
          
          drst_dxyz_p_mortar_quad_porder[i][j] =
            &schwarz_mortar_data->->drst_dxyz_p_mortar_quad_porder
            [matrix_stride + (i + j*(P4EST_DIM))*total_mortar_nodes_quad];
        }
      }

      if (ghost_mortar_side_data->faces_p == 0){

        d4est_element_data_t* e_m = &d4est_ghost->ghost_elements[gid];

        D4EST_ASSERT(faces_m == 1);
        D4EST_ASSERT(ghost_mortar_side_data->tree_m == e_m->tree);
        D4EST_ASSERT(ghost_mortar_side_data->tree_quadid_m[0] == e_m->tree_quadid);

        for (int i = 0; i < (P4EST_DIM); i++){
          xyz_on_f_m_quad [i] = &schwarz_mortar_data->xyz_on_f_m_quad[boundary_vector_stride];
          xyz_on_f_m_lobatto [i] = &schwarz_mortar_data->xyz_on_f_m_lobatto[boundary_vector_stride];
        }

        double* xyz_lobatto [(P4EST_DIM)];
        D4EST_ALLOC_DIM_VEC(xyz_lobatto, d4est_lgl_get_nodes((P4EST_DIM),ghost_mortar_side_data->deg_m[0]));
        
        d4est_geometry_compute_xyz
          (
           d4est_ops,
           d4est_geom,
           rst_points_lobatto,
           ghost_mortar_side_data->tree_m,
           ghost_mortar_side_data->deg_m[0],
           ghost_mortar_side_data->q_m[0],
           ghost_mortar_side_data->dq_m[0],
           xyz_lobatto
          );
  
        for (int d = 0; d < (P4EST_DIM); d++){

          d4est_operators_apply_slicer(d4est_ops,
                                       xyz_lobatto[d],
                                       (P4EST_DIM),
                                       ghost_mortar_side_data->f_m,
                                       ghost_mortar_side_data->deg_m[0],
                                       xyz_on_f_m_lobatto[d]);

        }

        P4EST_FREE(xyz_lobatto);
 
  
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
           &e_m->deg,
           &e_m->deg_quad,
           f_m,
           xyz_on_f_m_quad,
           drst_dxyz_m_mortar_quad,
           sj_m_mortar_quad,
           n_m_mortar_quad,
           NULL,
           j_div_sj_m_mortar_quad,
           COMPUTE_NORMAL_USING_JACOBIAN
          );

        d4est_mesh_calculate_mortar_h
          (
           p4est,
           &d4est_ghost->ghost_elements[gid],
           ghost_mortar_side_data->f_m,
           d4est_ops,
           d4est_geom,
           d4est_quad,
           face_h_type,
           j_div_sj_m_mortar_quad,
           hm_mortar_quad,
           1,
           1,
           &total_mortar_nodes_quad,
           d4est_mesh_get_size_parameters(d4est_factors)
          );

      }
      else {
        
        int faces_m = ghost_mortar_side_data->faces_m;
        int faces_p = ghost_mortar_side_data->faces_p;
        int orientation = ghost_mortar_side_data->orientation;
        int f_m = ghost_mortar_side_data->f_m;
        int f_p = ghost_mortar_side_data->f_p;
        int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;

        int deg_mortar_quad_porder [P4EST_HALF];
        int deg_mortar_lobatto_porder [P4EST_HALF];
        int deg_mortar_quad_morder [P4EST_HALF];
        int deg_mortar_lobatto_morder [P4EST_HALF];
        
        
        for (int i = 0; i < faces_m; i++){
          for (int j = 0; j < faces_p; j++){
            /* find max degree for each face pair of the two sides*/
            
            int jnew = j;
            if (faces_p == (P4EST_HALF)){
              jnew = d4est_reference_reorient_face_order((P4EST_DIM)-1,
                                                         f_m,
                                                         f_p,
                                                         orientation,
                                                         j);
            }
            deg_mortar_quad_morder[i + j] = d4est_util_max_int(ghost_mortar_side_data->deg_m_quad[i],
                                                               ghost_mortar_side_data->deg_p_quad[jnew]);
            deg_mortar_lobatto_morder[i + j] = d4est_util_max_int(ghost_mortar_side_data->deg_m[i],
                                                                  ghost_mortar_side_data->deg_p[jnew]);
            
          }
        }

        for (int i = 0; i < faces_mortar; i++){
          int inew = i;
          if (faces_mortar == P4EST_HALF){
            d4est_reference_reorient_face_order((P4EST_DIM)-1,
                                                f_m,
                                                f_p,
                                                orientation, i);
          }
          deg_mortar_quad_porder[inew] = deg_mortar_quad_morder[i];
          deg_mortar_lobatto_porder[inew] = deg_mortar_lobatto_morder[i];
        }
        
        d4est_mortars_compute_geometric_data_on_mortar
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           QUAD_INTEGRAND_UNKNOWN,
           ghost_mortar_side_data->tree_m,
           ghost_mortar_side_data->q_m[0],
           ghost_mortar_side_data->dq_m[0],
           ghost_mortar_side_data->mortar_side_id,
           ghost_mortar_side_data->faces_m,
           faces_mortar,
           &deg_mortar_lobatto[0],
           &deg_mortar_quad[0],
           f_m,
           NULL,
           drst_dxyz_m_mortar_quad,
           sj_m_mortar_quad,
           n_m_mortar_quad,
           NULL,
           j_div_sj_m_mortar_quad,
           COMPUTE_NORMAL_USING_JACOBIAN
          );

        d4est_mortars_compute_geometric_data_on_mortar
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           QUAD_INTEGRAND_UNKNOWN,
           ghost_mortar_side_data->tree_p,
           ghost_mortar_side_data->q_p[0],
           ghost_mortar_side_data->dq_p[0],
           ghost_mortar_side_data->mortar_side_id,
           ghost_mortar_side_data->faces_p,
           faces_mortar,
           &deg_mortar_lobatto_porder[0],
           &deg_mortar_quad_porder[0],
           ghost_mortar_side_data->f_p,
           NULL,
           drst_dxyz_p_mortar_quad_porder,
           NULL,
           NULL,
           NULL,
           j_div_sj_p_mortar_quad_porder,
           COMPUTE_NORMAL_USING_JACOBIAN
          );
  
        int face_mortar_stride = 0;
        for (int face = 0; face < faces_mortar; face++){
          int face_p = face;
          if (faces_mortar == (P4EST_HALF))
            face_p = d4est_reference_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, face);

          int oriented_face_mortar_stride = 0;
          for (int b = 0; b < face_p; b++){
            oriented_face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad_porder[b]);
          }

          d4est_operators_reorient_face_data
            (
             d4est_ops,
             &j_div_sj_p_mortar_quad_porder[oriented_face_mortar_stride],
             (P4EST_DIM)-1,
             deg_mortar_quad[face],
             ghost_mortar_side_data->orientation,
             ghost_mortar_side_data->f_m,
             ghost_mortar_side_data->f_p,
             &j_div_sj_p_mortar_quad[face_mortar_stride]
            );
    
          face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad[face]);
        }

        d4est_mesh_calculate_mortar_h
          (
           p4est,
           e_m,
           f_m,
           d4est_ops,
           d4est_geom,
           d4est_quad,
           *(d4est_mesh_face_h_t*)params,
           j_div_sj_m_mortar_quad,
           hm_mortar_quad,
           faces_mortar,
           faces_m,
           mortar_nodes_quad,
           d4est_mesh_get_size_parameters(d4est_factors)
          );
  
        d4est_mesh_calculate_mortar_h
          (
           p4est,
           &e_p_oriented[0],
           f_p,
           d4est_ops,
           d4est_geom,
           d4est_quad,
           *(d4est_mesh_face_h_t*)params,
           j_div_sj_p_mortar_quad,
           hp_mortar_quad,
           faces_mortar,
           faces_p,
           mortar_nodes_quad,
           d4est_mesh_get_size_parameters(d4est_factors)
          );
      }
      
    }
    
  }

}
