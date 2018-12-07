#include <pXest.h>
#include <d4est_ghost.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_operators.h>
#include <d4est_mortars.h>
#include <d4est_quadrature.h>
#include <d4est_quadrature_lobatto.h>
#include <d4est_solver_schwarz_mortar_data.h>

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


void
d4est_solver_schwarz_mortar_data_boundary_callback
(
 p4est_t* p4est,
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* params
)
{
  d4est_solver_schwarz_mortar_data_t* schwarz_mortar_data = params;
  int total_mortar_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM)-1, e_m->deg_quad);
  /* printf("rank %d: mortar_side_id_m = %d, stride = %d\n", p4est->mpirank, mortar_side_id_m, schwarz_mortar_data->stride); */


  d4est_mortar_side_data_t* mortar_side_data = &schwarz_mortar_data->mortar_side_data[schwarz_mortar_data->stride];
  mortar_side_data->boundary_quad_vector_stride = e_m->boundary_quad_vector_stride[f_m];
  mortar_side_data->mortar_quad_scalar_stride = e_m->mortar_quad_scalar_stride[f_m];
  mortar_side_data->mortar_quad_vector_stride = e_m->mortar_quad_vector_stride[f_m];
  mortar_side_data->mortar_quad_matrix_stride = e_m->mortar_quad_matrix_stride[f_m];
  mortar_side_data->faces_m = 1;
  mortar_side_data->faces_p = 0;
  mortar_side_data->f_p = -1;
  mortar_side_data->f_m = f_m;
  mortar_side_data->tree_p = -1;
  mortar_side_data->tree_m = e_m->tree;
  mortar_side_data->orientation = -1;
  mortar_side_data->q_m[0][0] = e_m->q[0];
  mortar_side_data->q_m[0][1] = e_m->q[1];
#if (P4EST_DIM)==3
  mortar_side_data->q_m[0][2] = e_m->q[2];
#endif
  mortar_side_data->dq_m[0] = e_m->dq;
  mortar_side_data->tree_quadid_m[0] = e_m->tree_quadid;
  mortar_side_data->deg_m[0] = e_m->deg;
  mortar_side_data->deg_quad_m[0] = e_m->deg_quad;
  mortar_side_data->mortar_side_id = schwarz_mortar_data->stride;
  schwarz_mortar_data->stride += 1;
}

void
d4est_solver_schwarz_mortar_data_interface_callback
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
 void* params
){

  int stride;
  int deg_p_lobatto [(P4EST_HALF)];
  int face_nodes_p_lobatto [(P4EST_HALF)];
  int face_nodes_p_quad [(P4EST_HALF)];
  int deg_m_lobatto [(P4EST_HALF)];
  int deg_m_quad [(P4EST_HALF)];
  int deg_p_quad [(P4EST_HALF)];
  int face_nodes_m_lobatto [(P4EST_HALF)];
  int face_nodes_m_quad [(P4EST_HALF)];
  int mortar_nodes_quad [(P4EST_HALF)];
  int mortar_nodes_lobatto [(P4EST_HALF)];
  int deg_mortar_quad [(P4EST_HALF)];
  int deg_mortar_lobatto [(P4EST_HALF)];
  int deg_p_lobatto_porder [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;

  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);

  /* calculate degs and nodes of each face of (-) side */
  int total_side_nodes_m_lobatto = 0;
  int total_side_nodes_m_quad = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_lobatto[i] = e_m[i]->deg;
    deg_m_quad[i] = e_m[i]->deg_quad;
    
    face_nodes_m_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg);
    face_nodes_m_quad[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_m_quad[i]);
    
    total_side_nodes_m_lobatto += face_nodes_m_lobatto[i];
    total_side_nodes_m_quad += face_nodes_m_quad[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_lobatto = 0;
  int total_side_nodes_p_quad = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p_lobatto[i] = e_p_oriented[i]->deg;
    deg_p_quad[i] = e_p_oriented[i]->deg_quad;
    deg_p_lobatto_porder[i] = e_p[i]->deg;

    face_nodes_p_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    face_nodes_p_quad[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_p_quad[i]);
    
    total_side_nodes_p_lobatto += face_nodes_p_lobatto[i];
    total_side_nodes_p_quad += face_nodes_p_quad[i];
  }

  /* calculate degs and nodes of the mortar faces */
  int total_mortar_nodes_quad = 0;
  int total_mortar_nodes_lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_quad[i+j] = d4est_util_max_int(deg_m_quad[i],
                                          deg_p_quad[j]);
      deg_mortar_lobatto[i+j] = d4est_util_max_int( e_m[i]->deg,
                                              e_p_oriented[j]->deg );
      mortar_nodes_quad[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_quad[i+j] );
      mortar_nodes_lobatto[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_lobatto[i+j] );
      total_mortar_nodes_quad += mortar_nodes_quad[i+j];
      total_mortar_nodes_lobatto += mortar_nodes_lobatto[i+j];

    }

  int deg_mortar_quad_porder [(P4EST_HALF)];
  int deg_mortar_lobatto_porder [(P4EST_HALF)];
  int mortar_nodes_quad_porder [(P4EST_HALF)];
  
  for(int i = 0; i < faces_mortar; i++){
    int inew = i;
    if (faces_mortar == (P4EST_HALF)){
      inew = d4est_reference_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, i);
    }
    deg_mortar_quad_porder[inew] = deg_mortar_quad[i];
    deg_mortar_lobatto_porder[inew] = deg_mortar_lobatto[i];
    mortar_nodes_quad_porder[inew] = mortar_nodes_quad[i];
  }
  d4est_solver_schwarz_mortar_data_t* schwarz_mortar_data = params;
  /* printf("rank %d: mortar_side_id_m = %d, stride = %d\n", p4est->mpirank, mortar_side_id_m, schwarz_mortar_data->stride); */
  d4est_mortar_side_data_t* mortar_side_data = &schwarz_mortar_data->mortar_side_data[schwarz_mortar_data->stride];

  mortar_side_data->boundary_quad_vector_stride = -1;
  mortar_side_data->mortar_quad_scalar_stride = e_m[0]->mortar_quad_scalar_stride[f_m];
  mortar_side_data->mortar_quad_vector_stride = e_m[0]->mortar_quad_vector_stride[f_m];
  mortar_side_data->mortar_quad_matrix_stride = e_m[0]->mortar_quad_matrix_stride[f_m];
  
  mortar_side_data->faces_m = faces_m;
  mortar_side_data->faces_p = faces_p;
  mortar_side_data->f_p = f_p;
  mortar_side_data->f_m = f_m;
  mortar_side_data->tree_p = e_p[0]->tree;
  mortar_side_data->tree_m = e_m[0]->tree;
  mortar_side_data->orientation = orientation;

  for (int f = 0; f < faces_m; f++){
    mortar_side_data->q_m[f][0] = e_m[f]->q[0];
    mortar_side_data->q_m[f][1] = e_m[f]->q[1];
#if (P4EST_DIM)==3
    mortar_side_data->q_m[f][2] = e_m[f]->q[2];
#endif
    mortar_side_data->dq_m[f] = e_m[f]->dq;
    mortar_side_data->tree_quadid_m[f] = e_m[f]->tree_quadid;
    mortar_side_data->deg_m[f] = e_m[f]->deg;
    mortar_side_data->deg_quad_m[f] = e_m[f]->deg_quad;
  }

  for (int f = 0; f < faces_p; f++){
    mortar_side_data->q_p[f][0] = e_p[f]->q[0];
    mortar_side_data->q_p[f][1] = e_p[f]->q[1];
#if (P4EST_DIM)==3
    mortar_side_data->q_p[f][2] = e_p[f]->q[2];
#endif
    mortar_side_data->dq_p[f] = e_p[f]->dq;
    mortar_side_data->tree_quadid_p[f] = e_p[f]->tree_quadid;
    mortar_side_data->deg_p[f] = e_p[f]->deg;
    mortar_side_data->deg_quad_p[f] = e_p[f]->deg_quad;
  }

  mortar_side_data->mortar_side_id = stride;
  schwarz_mortar_data->stride += 1;
}

void
d4est_solver_schwarz_mortar_data_destroy
(
 d4est_solver_schwarz_mortar_data_t* schwarz_mortar_data
)
{

  d4est_ghost_data_ext_destroy(schwarz_mortar_data->mortar_side_ghost_data);

  P4EST_FREE(schwarz_mortar_data->drst_dxyz_m_mortar_quad);
  P4EST_FREE(schwarz_mortar_data->drst_dxyz_p_mortar_quad_porder);
  P4EST_FREE(schwarz_mortar_data->sj_m_mortar_quad);
  P4EST_FREE(schwarz_mortar_data->n_m_mortar_quad);
  P4EST_FREE(schwarz_mortar_data->hm_mortar_quad);
  P4EST_FREE(schwarz_mortar_data->hp_mortar_quad);
  P4EST_FREE(schwarz_mortar_data->xyz_on_f_m_quad);
  P4EST_FREE(schwarz_mortar_data->xyz_on_f_m_lobatto);
  
  P4EST_FREE(schwarz_mortar_data->mortar_side_data);
  P4EST_FREE(schwarz_mortar_data);
}

d4est_solver_schwarz_mortar_data_t*
d4est_solver_schwarz_mortar_data_init
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* d4est_ghost,
 d4est_mesh_face_h_t face_h_type
)
{
  d4est_solver_schwarz_mortar_data_t* schwarz_mortar_data =
    P4EST_ALLOC(d4est_solver_schwarz_mortar_data_t, 1);

  /* printf(" d4est_factors->local_sizes.local_mortar_sides = %d\n",  d4est_factors->local_sizes.local_mortar_sides); */
  schwarz_mortar_data->mortar_side_data = P4EST_ALLOC(d4est_mortar_side_data_t, d4est_factors->local_sizes.local_mortar_sides);


  schwarz_mortar_data->stride = 0;
  d4est_mortars_fcn_ptrs_t flux_fcns;
  flux_fcns.flux_interface_fcn = d4est_solver_schwarz_mortar_data_interface_callback;
  flux_fcns.flux_boundary_fcn = d4est_solver_schwarz_mortar_data_boundary_callback;
  flux_fcns.user_ctx = (void*)schwarz_mortar_data;

  d4est_mortars_compute_flux_on_local_elements
    (
     p4est,
     d4est_ghost,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &flux_fcns
    );

  

  d4est_mortar_side_data_t* mirror_data = P4EST_ALLOC(d4est_mortar_side_data_t, (P4EST_FACES)*d4est_ghost->ghost->mirrors.elem_count);

  for (int i = 0; i < d4est_ghost->ghost->mirrors.elem_count; i++){
    /* get mirror data */
    d4est_element_data_t* mirror_ed_ref = &d4est_ghost->mirror_elements[i];
    for (int f = 0; f < (P4EST_FACES); f++){


      d4est_element_data_t* mirror_ed = d4est_element_data_get_ptr
                                        (
                                         p4est,
                                         mirror_ed_ref->tree,
                                         mirror_ed_ref->tree_quadid
                                        );
      int mortar_side_id = mirror_ed->face_belongs_to_which_mortar[f];            
      mirror_data[i*(P4EST_FACES) + f] = schwarz_mortar_data->mortar_side_data[mortar_side_id];

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

  P4EST_FREE(mirror_data);
  
  /* loop through data and set strides */
  
  /* int mortar_quad_ghost_size = 0; */
  /* int boundary_quad_ghost_size = 0; */
  
  int boundary_quad_stride = 0;
  int mortar_quad_stride = 0;

  for (int gid = 0; gid < d4est_ghost->ghost->ghosts.elem_count; gid++){
    
    d4est_mortar_side_data_t* ghost_mortar_side_data =
      d4est_ghost_data_ext_get_field_on_element
      (
       &d4est_ghost->ghost_elements[gid],
       0,
       schwarz_mortar_data->mortar_side_ghost_data
      );

    for (int f = 0; f < (P4EST_FACES); f++){
      int total_mortar_nodes_quad = 0;
 
      if (ghost_mortar_side_data[f].faces_p == 0){
        ghost_mortar_side_data[f].boundary_quad_vector_stride = (P4EST_DIM)*boundary_quad_stride;
        ghost_mortar_side_data[f].mortar_quad_scalar_stride = mortar_quad_stride;
        ghost_mortar_side_data[f].mortar_quad_vector_stride = (P4EST_DIM)*mortar_quad_stride;
        ghost_mortar_side_data[f].mortar_quad_matrix_stride
          = (P4EST_DIM)*(P4EST_DIM)*mortar_quad_stride;
        total_mortar_nodes_quad += d4est_lgl_get_nodes((P4EST_DIM)-1, ghost_mortar_side_data[f].deg_quad_m[0]);
        boundary_quad_stride += total_mortar_nodes_quad;
      }
      else {
        for (int i = 0; i < ghost_mortar_side_data[f].faces_m; i++){
          for (int j = 0; j < ghost_mortar_side_data[f].faces_p; j++){
            /* find max degree for each face pair of the two sides*/
            
            int jnew = j;
            if (ghost_mortar_side_data[f].faces_p == (P4EST_HALF)){
              jnew = d4est_reference_reorient_face_order((P4EST_DIM)-1,
                                                         ghost_mortar_side_data[f].f_m,
                                                         ghost_mortar_side_data[f].f_p,
                                                         ghost_mortar_side_data[f].orientation, j);
            }
            int deg_mortar_quad = d4est_util_max_int(ghost_mortar_side_data[f].deg_quad_m[i],
                                                 ghost_mortar_side_data[f].deg_quad_p[jnew]);
            int mortar_nodes_quad = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_quad);
            total_mortar_nodes_quad += mortar_nodes_quad;
          }
        }
        ghost_mortar_side_data[f].mortar_quad_scalar_stride = mortar_quad_stride;
        ghost_mortar_side_data[f].mortar_quad_vector_stride = (P4EST_DIM)*mortar_quad_stride;
        ghost_mortar_side_data[f].mortar_quad_matrix_stride
          = (P4EST_DIM)*(P4EST_DIM)*mortar_quad_stride;
        ghost_mortar_side_data[f].boundary_quad_vector_stride = -1;
      }
      
      mortar_quad_stride += total_mortar_nodes_quad;

      printf("ghost_mortar_side_data[f].faces_m = %d\n", ghost_mortar_side_data[f].faces_m);
      printf("ghost_mortar_side_data[f].faces_p = %d\n", ghost_mortar_side_data[f].faces_p);
      printf("ghost_mortar_side_data[f].boundary_vector_stride = %d\n", ghost_mortar_side_data[f].boundary_quad_vector_stride);
      ghost_mortar_side_data[f].total_mortar_nodes_quad = total_mortar_nodes_quad;
      /* ghost_mortar_side_data++; */
    }
  }

  int mortar_quad_ghost_size = mortar_quad_stride;
  int boundary_quad_ghost_size = boundary_quad_stride;

  printf("mortar_quad_ghost_size = %d\n", mortar_quad_ghost_size);
  
  schwarz_mortar_data->sj_m_mortar_quad = P4EST_ALLOC_ZERO(double, mortar_quad_ghost_size);
  schwarz_mortar_data->n_m_mortar_quad = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*mortar_quad_ghost_size);
  schwarz_mortar_data->hp_mortar_quad = P4EST_ALLOC_ZERO(double, mortar_quad_ghost_size);
  schwarz_mortar_data->hm_mortar_quad = P4EST_ALLOC_ZERO(double, mortar_quad_ghost_size);
  schwarz_mortar_data->drst_dxyz_m_mortar_quad = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*(P4EST_DIM)*mortar_quad_ghost_size);
  schwarz_mortar_data->drst_dxyz_p_mortar_quad_porder = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*(P4EST_DIM)*mortar_quad_ghost_size);
  schwarz_mortar_data->xyz_on_f_m_quad = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*boundary_quad_ghost_size);
  schwarz_mortar_data->xyz_on_f_m_lobatto = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*boundary_quad_ghost_size);

  for (int gid = 0; gid < d4est_ghost->ghost->ghosts.elem_count; gid++){

    d4est_mortar_side_data_t* ghost_mortar_side_data =
      d4est_ghost_data_ext_get_field_on_element
      (
       &d4est_ghost->ghost_elements[gid],
       0,
       schwarz_mortar_data->mortar_side_ghost_data
      );
    
    for (int f = 0; f < (P4EST_FACES); f++){

      int scalar_stride = ghost_mortar_side_data[f].mortar_quad_scalar_stride;
      int vector_stride = ghost_mortar_side_data[f].mortar_quad_vector_stride;
      int matrix_stride = ghost_mortar_side_data[f].mortar_quad_matrix_stride;
      int boundary_vector_stride = ghost_mortar_side_data[f].boundary_quad_vector_stride;


      printf("boundary_vector_stride = %d\n", boundary_vector_stride);
      double* xyz_on_f_m_quad [(P4EST_DIM)];
      double* xyz_on_f_m_lobatto [(P4EST_DIM)];
      double* n_m_mortar_quad [(P4EST_DIM)];
      double* drst_dxyz_m_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
      double* drst_dxyz_p_mortar_quad_porder [(P4EST_DIM)][(P4EST_DIM)];

      printf("scalar stride = %d\n", scalar_stride);
      double* sj_m_mortar_quad = &schwarz_mortar_data->sj_m_mortar_quad[scalar_stride];
      double* hm_mortar_quad = &schwarz_mortar_data->hm_mortar_quad[scalar_stride];
      double* hp_mortar_quad = &schwarz_mortar_data->hp_mortar_quad[scalar_stride];
  
      for (int i = 0; i < (P4EST_DIM); i++){
        n_m_mortar_quad[i] =
          &schwarz_mortar_data->n_m_mortar_quad
          [vector_stride + i*ghost_mortar_side_data[f].total_mortar_nodes_quad];

        for (int j = 0; j < (P4EST_DIM); j++){
          drst_dxyz_m_mortar_quad[i][j] =
            &schwarz_mortar_data->drst_dxyz_m_mortar_quad
            [matrix_stride + (i + j*(P4EST_DIM))*ghost_mortar_side_data[f].total_mortar_nodes_quad];
          
          drst_dxyz_p_mortar_quad_porder[i][j] =
            &schwarz_mortar_data->drst_dxyz_p_mortar_quad_porder
            [matrix_stride + (i + j*(P4EST_DIM))*ghost_mortar_side_data[f].total_mortar_nodes_quad];
        }
      }

      if (ghost_mortar_side_data[f].faces_p == 0){

        d4est_element_data_t* e_m = &d4est_ghost->ghost_elements[gid];
        D4EST_ASSERT(ghost_mortar_side_data[f].faces_m == 1);
        D4EST_ASSERT(ghost_mortar_side_data[f].tree_m == e_m->tree);
        D4EST_ASSERT(ghost_mortar_side_data[f].tree_quadid_m[0] == e_m->tree_quadid);

        for (int i = 0; i < (P4EST_DIM); i++){
          xyz_on_f_m_quad [i] = &schwarz_mortar_data->xyz_on_f_m_quad[boundary_vector_stride];
          xyz_on_f_m_lobatto [i] = &schwarz_mortar_data->xyz_on_f_m_lobatto[boundary_vector_stride];
        }

        double* xyz_lobatto [(P4EST_DIM)];
        D4EST_ALLOC_DIM_VEC(xyz_lobatto, d4est_lgl_get_nodes((P4EST_DIM),ghost_mortar_side_data[f].deg_m[0]));


        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq =  e_m->dq;
        mesh_object.tree = e_m->tree;
        
        mesh_object.q[0] = e_m->q[0];
        mesh_object.q[1] = e_m->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = e_m->q[2];
#endif

        

        d4est_rst_t rst_points_lobatto;
        rst_points_lobatto.r = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, e_m->deg, 0);
        rst_points_lobatto.s = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, e_m->deg, 1);
        rst_points_lobatto.t = NULL;
#if (P4EST_DIM)==3
        rst_points_lobatto.t = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, e_m->deg, 2);
#endif
        
        
        d4est_geometry_compute_xyz
          (
           d4est_ops,
           d4est_geom,
           rst_points_lobatto,
           ghost_mortar_side_data[f].tree_m,
           ghost_mortar_side_data[f].deg_m[0],
           ghost_mortar_side_data[f].q_m[0],
           ghost_mortar_side_data[f].dq_m[0],
           xyz_lobatto
          );
  
        for (int d = 0; d < (P4EST_DIM); d++){

          d4est_operators_apply_slicer(d4est_ops,
                                       xyz_lobatto[d],
                                       (P4EST_DIM),
                                       ghost_mortar_side_data[f].f_m,
                                       ghost_mortar_side_data[f].deg_m[0],
                                       xyz_on_f_m_lobatto[d]);

        }

        D4EST_FREE_DIM_VEC(xyz_lobatto);
 
        double* j_div_sj_m_mortar_quad = P4EST_ALLOC(double, ghost_mortar_side_data[f].total_mortar_nodes_quad);
        d4est_mortars_compute_geometric_data_on_mortar
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           QUAD_INTEGRAND_UNKNOWN,
           e_m->tree,
           e_m->q,
           e_m->dq,
           ghost_mortar_side_data[f].mortar_side_id,
           1,
           1,
           &e_m->deg,
           &e_m->deg_quad,
           ghost_mortar_side_data[f].f_m,
           xyz_on_f_m_quad,
           drst_dxyz_m_mortar_quad,
           sj_m_mortar_quad,
           n_m_mortar_quad,
           NULL,
           j_div_sj_m_mortar_quad,
           COMPUTE_NORMAL_USING_JACOBIAN
          );


        DEBUG_PRINT_ARR_DBL(j_div_sj_m_mortar_quad, ghost_mortar_side_data[f].total_mortar_nodes_quad);
        d4est_element_data_t* ged = &d4est_ghost->ghost_elements[gid];
        
        d4est_mesh_calculate_mortar_h
          (
           p4est,
           &ged,
           ghost_mortar_side_data[f].f_m,
           d4est_ops,
           d4est_geom,
           d4est_quad,
           face_h_type,
           j_div_sj_m_mortar_quad,
           hm_mortar_quad,
           1,
           1,
           &ghost_mortar_side_data[f].total_mortar_nodes_quad,
           d4est_mesh_get_size_parameters(d4est_factors)
          );

        P4EST_FREE(j_div_sj_m_mortar_quad);
        
      }
      else {

        int faces_m = ghost_mortar_side_data[f].faces_m;
        int faces_p = ghost_mortar_side_data[f].faces_p;
        int orientation = ghost_mortar_side_data[f].orientation;
        int f_m = ghost_mortar_side_data[f].f_m;
        int f_p = ghost_mortar_side_data[f].f_p;
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
            deg_mortar_quad_morder[i + j] = d4est_util_max_int(ghost_mortar_side_data[f].deg_quad_m[i],
                                                               ghost_mortar_side_data[f].deg_quad_p[jnew]);
            deg_mortar_lobatto_morder[i + j] = d4est_util_max_int(ghost_mortar_side_data[f].deg_m[i],
                                                                  ghost_mortar_side_data[f].deg_p[jnew]);
            
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
        double* j_div_sj_m_mortar_quad = P4EST_ALLOC(double, ghost_mortar_side_data[f].total_mortar_nodes_quad);
        double* j_div_sj_p_mortar_quad = P4EST_ALLOC(double, ghost_mortar_side_data[f].total_mortar_nodes_quad);
        double* j_div_sj_p_mortar_quad_porder = P4EST_ALLOC(double, ghost_mortar_side_data[f].total_mortar_nodes_quad);
        
        d4est_mortars_compute_geometric_data_on_mortar
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           QUAD_INTEGRAND_UNKNOWN,
           ghost_mortar_side_data[f].tree_m,
           ghost_mortar_side_data[f].q_m[0],
           ghost_mortar_side_data[f].dq_m[0],
           ghost_mortar_side_data[f].mortar_side_id,
           ghost_mortar_side_data[f].faces_m,
           faces_mortar,
           &deg_mortar_lobatto_morder[0],
           &deg_mortar_quad_morder[0],
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
           ghost_mortar_side_data[f].tree_p,
           ghost_mortar_side_data[f].q_p[0],
           ghost_mortar_side_data[f].dq_p[0],
           ghost_mortar_side_data[f].mortar_side_id,
           ghost_mortar_side_data[f].faces_p,
           faces_mortar,
           &deg_mortar_lobatto_porder[0],
           &deg_mortar_quad_porder[0],
           ghost_mortar_side_data[f].f_p,
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
             deg_mortar_quad_morder[face],
             ghost_mortar_side_data[f].orientation,
             ghost_mortar_side_data[f].f_m,
             ghost_mortar_side_data[f].f_p,
             &j_div_sj_p_mortar_quad[face_mortar_stride]
            );
    
          face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad_morder[face]);
        }
        
        /* d4est_mesh_calculate_mortar_h */
        /*   ( */
        /*    p4est, */
        /*    e_m, */
        /*    f_m, */
        /*    d4est_ops, */
        /*    d4est_geom, */
        /*    d4est_quad, */
        /*    face_h_type, */
        /*    j_div_sj_m_mortar_quad, */
        /*    hm_mortar_quad, */
        /*    faces_mortar, */
        /*    faces_m, */
        /*    mortar_nodes_quad, */
        /*    d4est_mesh_get_size_parameters(d4est_factors) */
        /*   ); */
  
      /*   d4est_mesh_calculate_mortar_h */
      /*     ( */
      /*      p4est, */
      /*      &e_p_oriented[0], */
      /*      f_p, */
      /*      d4est_ops, */
      /*      d4est_geom, */
      /*      d4est_quad, */
      /*      face_h_type, */
      /*      j_div_sj_p_mortar_quad, */
      /*      hp_mortar_quad, */
      /*      faces_mortar, */
      /*      faces_p, */
      /*      mortar_nodes_quad, */
      /*      d4est_mesh_get_size_parameters(d4est_factors) */
      /*     ); */
      /* } */

        if (face_h_type == FACE_H_EQ_J_DIV_SJ_QUAD){
          /* D4EST_ABORT("H_EQ_J_DIV_SJ_QUAD is no longer supported"); */
          d4est_mesh_calculate_mortar_h_eq_j_div_sj_quad
            (
             j_div_sj_m_quad_mortar,
             faces_mortar,
             nodes_mortar_quad,
             hm_mortar_quad
            );

          d4est_mesh_calculate_mortar_h_eq_j_div_sj_quad
            (
             j_div_sj_m_quad_mortar,
             faces_mortar,
             nodes_mortar_quad,
             hm_mortar_quad
            );
        }
        else if (face_h_type == FACE_H_EQ_TREE_H){

        }
        else {
          D4EST_ABORT("face_h_type");
        }

        P4EST_FREE(j_div_sj_m_mortar_quad);
        P4EST_FREE(j_div_sj_p_mortar_quad);
        P4EST_FREE(j_div_sj_p_mortar_quad_porder);
        
      }
    }
  }

  return schwarz_mortar_data;
}


/* TODO: to check schwarz mortar data do a sum over 
sj mortar data in serial and parallel over all schwarz subdomains */
/* This file was automatically generated.  Do not edit! */
