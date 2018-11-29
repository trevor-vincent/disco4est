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


static int
d4est_solver_schwarz_mortar_data_init
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost
)
{
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

  int mortar_quad_ghost_size = 0;




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
      for (int i = 0; i < faces_m; i++)
        for (int j = 0; j < faces_p; j++){
          /* find max degree for each face pair of the two sides*/
          deg_mortar_quad = d4est_util_max_int(ghost_mortar_side_data->deg_m_quad[i],
                                               ghost_mortar_side_data->deg_p_quad[j]);
          mortar_nodes_quad = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_quad);
          total_mortar_nodes_quad += mortar_nodes_quad;
        }
      
      /* compute mortar side and add to total */
      mortar_quad_ghost_size += total_mortar_nodes_quad;
      ghost_mortar_side_data++;
    }
    
  }

  /* allocate mortar shit */
  schwarz_mortar_data->sj_m_mortar_quad = P4EST_ALLOC(double, mortar_quad_ghost_size);
  schwarz_mortar_data->n_m_mortar_quad = P4EST_ALLOC(double, (P4EST_DIM)*mortar_quad_ghost_size);
  schwarz_mortar_data->hp_mortar_quad = P4EST_ALLOC(double, mortar_quad_ghost_size);
  schwarz_mortar_data->hm_mortar_quad = P4EST_ALLOC(double, mortar_quad_ghost_size);
  schwarz_mortar_data->drst_dxyz_mortar_quad = P4EST_ALLOC(double, (P4EST_DIM)*(P4EST_DIM)*mortar_quad_ghost_size);

    /* compute mortar shit */



  for (int gid = 0; gid < num_ghost_elements; gid++){

    d4est_mortar_side_data_t* ghost_mortar_side_data =
      d4est_ghost_data_ext_get_field_on_element
      (
       &d4est_ghost->ghost_elements[gid],
       0,
       *d4est_ghost_data_ext
      );
    
    for (int f = 0; f < (P4EST_FACES); f++){
  
      double* n_m_mortar_quad [(P4EST_DIM)];
      double* drst_dxyz_m_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
      double* drst_dxyz_p_mortar_quad_porder [(P4EST_DIM)][(P4EST_DIM)];
      double* sj_m_mortar_quad = &d4est_factors->sj_m_mortar_quad[scalar_stride];
  
      for (int i = 0; i < (P4EST_DIM); i++){
        n_m_mortar_quad[i] = &d4est_factors->n_m_mortar_quad[vector_stride + i*total_mortar_nodes_quad];
        for (int j = 0; j < (P4EST_DIM); j++){
          drst_dxyz_m_mortar_quad[i][j] = &d4est_factors->drst_dxyz_m_mortar_quad[matrix_stride + (i + j*(P4EST_DIM))*total_mortar_nodes_quad];
          drst_dxyz_p_mortar_quad_porder[i][j] = &d4est_factors->drst_dxyz_p_mortar_quad_porder[matrix_stride + (i + j*(P4EST_DIM))*total_mortar_nodes_quad];
        }
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
     e_p[0]->tree,
     e_p[0]->q,
     e_p[0]->dq,
     mortar_side_id_p,
     faces_p,
     faces_mortar,
     &deg_mortar_lobatto_porder[0],
     &deg_mortar_quad_porder[0],
     f_p,
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
       orientation,
       f_m,
       f_p,
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

  /* d4est_element_data_t* e_p_oriented [(P4EST_HALF)]; */
  /* d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented); */
  
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
