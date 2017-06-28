#include <pXest.h>
#include <util.h>
#include <d4est_linalg.h>
#include <d4est_mortars.h>
#include <d4est_geometry.h>
#include <d4est_quadrature.h>

void
d4est_mortars_compute_geometric_data_on_mortar_aux
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 p4est_topidx_t e0_tree,
 p4est_qcoord_t e0_q [(P4EST_DIM)], /* qcoord of first element of side */
 p4est_qcoord_t e0_dq, /* qcoord vector spanning first element of side */
 d4est_rst_t* rst_quad_mortar,
 int num_faces_side, 
 int num_faces_mortar,
 int* deg_mortar_quad,
 int face,
 double* drst_dxyz_on_mortar_quad [(P4EST_DIM)][(P4EST_DIM)],
 double* sj_on_mortar_quad,
 double* n_on_mortar_quad [(P4EST_DIM)],
 double* n_sj_on_mortar_quad [(P4EST_DIM)],
 double* j_div_sj_mortar_quad,
 normal_compute_method_t n_compute_method
)
{  
  double* dxyz_drst_on_face_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_on_face_quad[(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_times_jac_on_face_quad[(P4EST_DIM)][(P4EST_DIM)];
  
  int max_deg = 0;
  for (int i = 0; i < num_faces_mortar; i++){
    max_deg = (deg_mortar_quad[i] > max_deg) ? deg_mortar_quad[i] : max_deg;
  }
  int volume_nodes_max = d4est_lgl_get_nodes((P4EST_DIM), max_deg);
  int face_nodes_max = d4est_lgl_get_nodes((P4EST_DIM)-1, max_deg);
  for (int i = 0; i < (P4EST_DIM); i++)
    for (int j = 0; j < (P4EST_DIM); j++){
      /* dxyz_drst[i][j] = P4EST_ALLOC(double, volume_nodes_max); */
      dxyz_drst_on_face_quad[i][j] = P4EST_ALLOC(double, face_nodes_max);
      drst_dxyz_times_jac_on_face_quad[i][j] = P4EST_ALLOC(double, face_nodes_max);
      drst_dxyz_on_face_quad[i][j] = P4EST_ALLOC(double, face_nodes_max);
    }

  double* temp = P4EST_ALLOC(double, volume_nodes_max);
  double* J_on_face_quad = P4EST_ALLOC(double, face_nodes_max);

  p4est_qcoord_t q0 [(P4EST_HALF)][(P4EST_DIM)];
  p4est_qcoord_t q [(P4EST_DIM)];
  p4est_qcoord_t mortar_dq;
  
  d4est_geometry_compute_qcoords_on_mortar
    (
     e0_tree,
     e0_q, /* qcoord of first element of side */
     e0_dq, /* qcoord vector spanning first element of side */
     num_faces_side, 
     num_faces_mortar,
     face,
     q0,
     &mortar_dq
    );
 
  int face_mortar_quad_stride = 0;
  for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){
  
    int face_mortar_quad_nodes = d4est_lgl_get_nodes((P4EST_DIM) - 1, deg_mortar_quad[face_mortar]);
    /* double* xyz [(P4EST_DIM)]; */
    for (int d = 0; d < (P4EST_DIM); d++){
      q[d] = q0[face_mortar][d];
    }

    if (d4est_geom->DX_compute_method == GEOM_COMPUTE_ANALYTIC){
      d4est_geometry_compute_dxyz_drst_face_analytic
        (
         d4est_ops,
         d4est_geom,
         rst_quad_mortar[face_mortar],
         q,
         mortar_dq,
         e0_tree,
         face,
         deg_mortar_quad[face_mortar],
         dxyz_drst_on_face_quad
        );
    }
    else {
      mpi_abort("[D4EST_ERROR]: compute method type for DX not supported");
    }

    if (d4est_geom->JAC_compute_method == GEOM_COMPUTE_NUMERICAL){
      d4est_geometry_compute_jacobian
        (
         dxyz_drst_on_face_quad,
         J_on_face_quad,
         face_mortar_quad_nodes
        );
    }
    else {
      mpi_abort("[D4EST_ERROR]: compute method type for JAC not supported");
    }

    d4est_geometry_compute_drst_dxyz
      (
       dxyz_drst_on_face_quad,
       J_on_face_quad,
       drst_dxyz_on_face_quad,
       face_mortar_quad_nodes
      );

    d4est_geometry_compute_drst_dxyz_times_jacobian
      (
       dxyz_drst_on_face_quad,
       drst_dxyz_times_jac_on_face_quad,
       face_mortar_quad_nodes
      );
    
      int i0 = -1; 
      if (face == 0 || face == 1){
        i0 = 0;
      }
      else if (face == 2 || face == 3){
        i0 = 1;
      }
      else if (face == 4 || face == 5){
        i0 = 2;
      }
      else {
        mpi_abort("face must be < 6\n");
      }
      double sgn = (face == 0 || face == 2 || face == 4) ? -1. : 1.;
      d4est_geometry_face_info_t face_info;
      d4est_geometry_get_face_info(face, &face_info);

      for (int i = 0; i < face_mortar_quad_nodes; i++){
        int is = face_mortar_quad_stride + i;
        double n_is [] = {0.,0.,0.};
        double n_sj_is [] = {0.,0.,0.};
        double sj_is = 0.;

        if(n_compute_method == COMPUTE_NORMAL_USING_JACOBIAN){
          double sgn = (face == 0 || face == 2 || face == 4) ? -1. : 1.;
          for (int d = 0; d < (P4EST_DIM); d++){
            /* n_is[d] = sgn*drst_dxyz_on_face_quad[i0][d][i]*J_on_face_quad[i]; */
            n_sj_is[d] = sgn*drst_dxyz_times_jac_on_face_quad[i0][d][i];
            sj_is += n_sj_is[d]*n_sj_is[d];
          }

          
          sj_is = sqrt(sj_is);
          for (int d = 0; d < (P4EST_DIM); d++){
            n_is[d] = n_sj_is[d]/sj_is;
          }
        }
        
        else if (n_compute_method == COMPUTE_NORMAL_USING_CROSS_PRODUCT){
          double sgn = (face == 0 || face == 3 || face == 4) ? -1. : 1.;
          double vecs [2][3] = {{0.,0.,0.},{0.,0.,1.}};
          
          for (int d = 0; d < (P4EST_DIM); d++)
            for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
              vecs[dir][d] =  dxyz_drst_on_face_quad[d][dir == 0 ? face_info.a : face_info.b][i];
            }
          
          d4est_linalg_cross_prod
            (
             vecs[0][0],
             vecs[0][1],
             vecs[0][2],
             vecs[1][0],
             vecs[1][1],
             vecs[1][2],
             &(n_sj_is[0]),
             &(n_sj_is[1]),
             &(n_sj_is[2])
            );


          sj_is = 0.;
          for (int d = 0; d < (P4EST_DIM); d++){
            /* The normals are backwards for these 2(3) face_mortars in 2-d(3-d) */
            n_sj_is[d] *= sgn;
            sj_is += n_sj_is[d]*n_sj_is[d];
          }
          sj_is = sqrt(sj_is);
          
          for (int d = 0; d < (P4EST_DIM); d++)
            n_is[d] = n_sj_is[d]/sj_is;
          
        }
        else {
          mpi_abort("[D4EST_ERROR]: Only two ways to compute normal");
        }

        /* STORE COMPUTATIONS */
        if (drst_dxyz_on_mortar_quad != NULL){
          for (int d1 = 0; d1 < (P4EST_DIM); d1++){
            for (int d2 = 0; d2 < (P4EST_DIM); d2++){
              drst_dxyz_on_mortar_quad[d1][d2][is] = drst_dxyz_on_face_quad[d1][d2][i];
            }
          }
        }

        if (n_sj_on_mortar_quad != NULL){
          for (int d = 0; d < (P4EST_DIM); d++){
            n_sj_on_mortar_quad[d][is] = n_sj_is[d];
          }
        }
        
        if (n_on_mortar_quad != NULL){
          for (int d = 0; d < (P4EST_DIM); d++){
            n_on_mortar_quad[d][is] = n_is[d];
          }
        }
        
        if (sj_on_mortar_quad != NULL){
          sj_on_mortar_quad[is] = sj_is;
        }
        
        if (j_div_sj_mortar_quad != NULL){
          j_div_sj_mortar_quad[is] = J_on_face_quad[i]/sj_is;
        }
      }  

    
    face_mortar_quad_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad[face_mortar]);
  }

  P4EST_FREE(temp);
  P4EST_FREE(J_on_face_quad);

  for (int i = 0; i < (P4EST_DIM); i++)
    for (int j = 0; j < (P4EST_DIM); j++){
      P4EST_FREE(dxyz_drst_on_face_quad[i][j]);
      P4EST_FREE(drst_dxyz_on_face_quad[i][j]);
      P4EST_FREE(drst_dxyz_times_jac_on_face_quad[i][j]);
    }
}

void
d4est_mortars_compute_geometric_data_on_mortar
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_quadrature_integrand_type_t integrand_type,
 p4est_topidx_t e0_tree,
 p4est_qcoord_t e0_q [(P4EST_DIM)], /* qcoord of first element of side */
 p4est_qcoord_t e0_dq, /* qcoord vector spanning first element of side */
 int mortar_side_id,
 int num_faces_side, 
 int num_faces_mortar,
 int* deg_mortar_quad,
 int face,
 double* drst_dxyz_on_mortar_quad [(P4EST_DIM)][(P4EST_DIM)],
 double* sj_on_mortar_quad,
 double* n_on_mortar_quad [(P4EST_DIM)],
 double* n_sj_on_mortar_quad [(P4EST_DIM)],
 double* j_div_sj_mortar_quad,
 normal_compute_method_t n_compute_method
)
{
  d4est_rst_t rst_quad_mortar [(P4EST_HALF)];
  
  d4est_geometry_compute_rst_on_mortar
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     integrand_type,
     &rst_quad_mortar[0],
     e0_tree,
     e0_q,
     e0_dq,
     num_faces_side,
     num_faces_mortar,
     face,
     mortar_side_id,
     deg_mortar_quad
    );

    d4est_mortars_compute_geometric_data_on_mortar_aux
    (
     d4est_ops,
     d4est_geom,
     e0_tree,
     e0_q, /* qcoord of first element of side */
     e0_dq, /* qcoord vector spanning first element of side */
     rst_quad_mortar,
     num_faces_side, 
     num_faces_mortar,
     deg_mortar_quad,
     face,
     drst_dxyz_on_mortar_quad,
     sj_on_mortar_quad,
     n_on_mortar_quad,
     n_sj_on_mortar_quad,
     j_div_sj_mortar_quad,
     n_compute_method
    );
  
}

void
d4est_geometry_compute_rst_on_mortar
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_quadrature_integrand_type_t integrand_type,
 d4est_rst_t* d4est_rst,
 p4est_topidx_t e0_tree,
 p4est_qcoord_t e0_q [(P4EST_DIM)], /* qcoord of first element of side */
 p4est_qcoord_t e0_dq, /* qcoord vector spanning first element of side */
 int num_faces_side, 
 int num_faces_mortar,
 int face,
 int mortar_side_id,
 int* deg_mortar_quad
)
{
  p4est_qcoord_t mortar_q0 [(P4EST_HALF)][(P4EST_DIM)];
  p4est_qcoord_t q [(P4EST_DIM)];
  p4est_qcoord_t mortar_dq;

  d4est_geometry_compute_qcoords_on_mortar
    (
     e0_tree,
     e0_q, /* qcoord of first element of side */
     e0_dq, /* qcoord vector spanning first element of side */
     num_faces_side, 
     num_faces_mortar,
     face,
     mortar_q0,
     &mortar_dq
    );

  for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){
  
    d4est_quadrature_mortar_t face_object;
    face_object.tree = e0_tree;
    face_object.mortar_side_id = mortar_side_id;
    face_object.face = face;
    face_object.dq = mortar_dq;
    face_object.mortar_side_id = mortar_side_id;
    face_object.mortar_subface_id = face_mortar;
    for (int d = 0; d < (P4EST_DIM); d++){
      face_object.q[d] = mortar_q0[face_mortar][d];
    }
    
    d4est_rst[face_mortar] = d4est_quadrature_get_rst_points
                             (
                              d4est_ops,
                              d4est_quad,
                              d4est_geom,
                              &face_object,
                              QUAD_OBJECT_MORTAR,
                              integrand_type,
                              deg_mortar_quad[face_mortar]
                             );    
  }
  
}

void
d4est_geometry_compute_qcoords_on_mortar
(
 p4est_topidx_t e0_tree,
 p4est_qcoord_t e0_q [(P4EST_DIM)], /* qcoord of first element of side */
 p4est_qcoord_t e0_dq, /* qcoord vector spanning first element of side */
 int num_faces_side, 
 int num_faces_mortar,
 int face,
 p4est_qcoord_t mortar_q0 [(P4EST_HALF)][(P4EST_DIM)],
 p4est_qcoord_t* mortar_dq
)
{
  for (int j = 0; j < (P4EST_HALF); j++){
    int c = p4est_face_corners[face][j];
    for (int d = 0; d < (P4EST_DIM); d++){
      int cd = d4est_reference_is_child_left_or_right(c, d);
      mortar_q0[j][d] = e0_q[d] + cd*e0_dq;
    }
  }
  
  p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)];
  
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
      dqa[dir][d] = (mortar_q0[(dir+1)][d] - mortar_q0[0][d]);
      if (num_faces_side != num_faces_mortar)
        dqa[dir][d] /= 2;
    }
  }    

  p4est_qcoord_t dq0mf0 [(P4EST_DIM)];
  for (int d = 0; d < (P4EST_DIM); d++){
    dq0mf0[d] = (mortar_q0[0][d] - e0_q[d])/2;
  }
    
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int c = 0; c < (P4EST_HALF); c++){
      mortar_q0[c][d] = e0_q[d];
      if (num_faces_side != num_faces_mortar)
        mortar_q0[c][d] += dq0mf0[d];
      for (int dir = 0; dir < (P4EST_DIM) - 1; dir++){
        int cd = d4est_reference_is_child_left_or_right(c, dir);
        mortar_q0[c][d] += cd*dqa[dir][d];
      }
    }
  }

  *mortar_dq = (num_faces_side == num_faces_mortar) ? e0_dq : e0_dq/2;
}


void d4est_mortars_project_mortar_onto_side(d4est_operators_t* d4est_ops,
                                     double* in_mortar, int faces_mortar,
                                     int* deg_mortar, double* out_side,
                                     int faces_side, int* deg_side) {
  /* most-common so first */
  /* p-prolong to mortar space */
  if (faces_side == 1 && faces_mortar == 1) {
    d4est_operators_apply_p_restrict(d4est_ops, in_mortar, deg_mortar[0], (P4EST_DIM)-1,
                            deg_side[0], out_side);
  }

  /* hp-prolong to mortar space. */
  else if (faces_side == 1 && faces_mortar == (P4EST_HALF)) {
    d4est_operators_apply_hp_restrict(d4est_ops, in_mortar, &deg_mortar[0], (P4EST_DIM)-1,
                             deg_side[0], out_side);
  }

  /* if faces_side = (P4EST_HALF), then faces_mortar == (P4EST_HALF) */
  /* here we loop through faces and p-prolong */
  else if (faces_side == (P4EST_HALF) && faces_mortar == (P4EST_HALF)) {
    int i;
    int stride_side = 0;
    int stride_mortar = 0;
    for (i = 0; i < faces_side; i++) {
      d4est_operators_apply_p_restrict(d4est_ops, &in_mortar[stride_mortar], deg_mortar[i],
                              (P4EST_DIM)-1, deg_side[i],
                              &out_side[stride_side]);

      stride_side += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_side[i]);
      stride_mortar += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar[i]);
    }
  }

  else {
    mpi_abort("ERROR: d4est_mortars_project_side_onto_mortar_space");
  }
}


void d4est_mortars_project_mass_mortar_onto_side(d4est_operators_t* dgmath,
                                          double* in_mortar, int faces_mortar,
                                          int* deg_mortar, double* out_side,
                                          int faces_side, int* deg_side) {

  /* most-common so first */
  /* p-prolong to mortar space */
  if (faces_side == 1 && faces_mortar == 1) {
    d4est_operators_apply_p_prolong_transpose(dgmath, in_mortar, deg_mortar[0], (P4EST_DIM)-1,
                            deg_side[0], out_side);
  }

  /* hp-prolong to mortar space. */
  else if (faces_side == 1 && faces_mortar == (P4EST_HALF)) {
    d4est_operators_apply_hp_prolong_transpose(dgmath, in_mortar, &deg_mortar[0], (P4EST_DIM)-1,
                             deg_side[0], out_side);
  }

  /* if faces_side = (P4EST_HALF), then faces_mortar == (P4EST_HALF) */
  /* here we loop through faces and p-prolong */
  else if (faces_side == (P4EST_HALF) && faces_mortar == (P4EST_HALF)) {
    int i;
    int stride_side = 0;
    int stride_mortar = 0;
    for (i = 0; i < faces_side; i++) {
      d4est_operators_apply_p_prolong_transpose(dgmath, &in_mortar[stride_mortar], deg_mortar[i],
                              (P4EST_DIM)-1, deg_side[i],
                              &out_side[stride_side]);

      stride_side += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_side[i]);
      stride_mortar += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar[i]);
    }
  }

  else {
    mpi_abort("ERROR: d4est_mortars_project_mass_mortar_onto_side_space");
  }
}


void d4est_mortars_project_side_onto_mortar_space(d4est_operators_t* d4est_ops,
                                           double* in_side, int faces_side,
                                           int* deg_side, double* out_mortar,
                                           int faces_mortar, int* deg_mortar) {
  /* most-common so first */
  /* p-prolong to mortar space */
  if (faces_side == 1 && faces_mortar == 1) {
    d4est_operators_apply_p_prolong
      (
       d4est_ops,
       in_side,
       deg_side[0],
       (P4EST_DIM)-1,
       deg_mortar[0],
       out_mortar
      );
  }

  /* hp-prolong to mortar space. */
  else if (faces_side == 1 && faces_mortar == (P4EST_HALF)) {
    d4est_operators_apply_hp_prolong
      (
       d4est_ops,
       in_side,
       deg_side[0],
       (P4EST_DIM)-1,
       &deg_mortar[0],
       out_mortar
      );
  }

  /* if faces_side = (P4EST_HALF), then faces_mortar == (P4EST_HALF) */
  /* here we loop through faces and p-prolong */
  else if (faces_side == (P4EST_HALF) && faces_mortar == (P4EST_HALF)) {
    int i;
    int stride_side = 0;
    int stride_mortar = 0;
    for (i = 0; i < faces_side; i++) {
      d4est_operators_apply_p_prolong(d4est_ops, &in_side[stride_side], deg_side[i],
                             (P4EST_DIM)-1, deg_mortar[i],
                             &out_mortar[stride_mortar]);

      stride_side += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_side[i]);
      stride_mortar += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar[i]);
    }
  } else {
    mpi_abort("ERROR: d4est_mortars_project_side_onto_mortar_space");
  }
}
