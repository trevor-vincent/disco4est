#include <pXest.h>
#include <d4est_solver_schwarz_laplacian.h>
#include <d4est_solver_schwarz_helpers.h>
#include <d4est_laplacian.h>
#include <d4est_laplacian_flux.h>

static
int
convert_mortar_side_data_to_elements
(
 p4est_t* p4est,
 d4est_solver_schwarz_subdomain_metadata_t* sub_data,
 d4est_mortar_side_data_t* mortar_data,
 d4est_element_data_t** e_m,
 d4est_element_data_t** e_p,
 int* zero_and_skip_m,
 int* zero_and_skip_p,
 int* element_m_sub_id
)
{
  int skip_p_sum = 0;
  
  for (int i = 0; i < mortar_data->faces_m; i++){
    int mortar_side_id = mortar_data->mortar_side_id;
    int f_m = mortar_data->f_m;

    e_m[0]->boundary_quad_vector_stride[f_m] = -1;
    e_m[0]->mortar_quad_vector_stride[f_m] = -1;
    e_m[0]->mortar_quad_scalar_stride[f_m] = -1;
    e_m[0]->mortar_quad_matrix_stride[f_m] = -1;
    
    e_m[i]->q[0] = mortar_data->q_m[i][0];
    e_m[i]->q[1] = mortar_data->q_m[i][1];
#if (P4EST_DIM)==3
    e_m[i]->q[2] = mortar_data->q_m[i][2];
#endif
    e_m[i]->dq = mortar_data->dq_m[i];
    e_m[i]->deg = mortar_data->deg_m[i];
    e_m[i]->deg_quad = mortar_data->deg_quad_m[i];
    e_m[i]->tree = mortar_data->tree_m;
    e_m[i]->tree_quadid = mortar_data->tree_quadid_m[i];
    e_m[i]->face_belongs_to_which_mortar[f_m] = mortar_side_id;
    e_m[i]->mpirank = p4est->mpirank;
    e_m[i]->id = -1;
    e_m[i]->sqr_nodal_stride = -1;
    e_m[i]->quad_stride = -1;
    int id = d4est_solver_schwarz_is_element_in_subdomain
             (
              sub_data,
              e_m[i]->tree,
              e_m[i]->tree_quadid
             );
    if (id >= 0){
      element_m_sub_id[i] = id;
      zero_and_skip_m[i] = 0;
      e_m[i]->nodal_stride = sub_data->element_metadata[id].nodal_stride;
    }
    else {
      element_m_sub_id[i] = -1;
      zero_and_skip_m[i] = 1;
      e_m[i]->nodal_stride = -1;
    }
  }
  for (int i = 0; i < mortar_data->faces_p; i++){
    int mortar_side_id = mortar_data->mortar_side_id;
    int f_p = mortar_data->f_p;
    e_p[0]->boundary_quad_vector_stride[f_p] = -1;
    e_p[0]->mortar_quad_vector_stride[f_p] = -1;
    e_p[0]->mortar_quad_scalar_stride[f_p] = -1;
    e_p[0]->mortar_quad_matrix_stride[f_p] = -1;
    e_p[i]->q[0] = mortar_data->q_p[i][0];
    e_p[i]->q[1] = mortar_data->q_p[i][1];
#if (P4EST_DIM)==3
    e_p[i]->q[2] = mortar_data->q_p[i][2];
#endif
    e_p[i]->dq = mortar_data->dq_p[i];
    e_p[i]->deg = mortar_data->deg_p[i];
    e_p[i]->deg_quad = mortar_data->deg_quad_p[i];
    e_p[i]->tree = mortar_data->tree_p;
    e_p[i]->tree_quadid = mortar_data->tree_quadid_p[i];
    e_p[i]->face_belongs_to_which_mortar[f_p] = -1;
    e_p[i]->mpirank = p4est->mpirank;
    e_p[i]->id = -1;
    e_p[i]->sqr_nodal_stride = -1;
    e_p[i]->quad_stride = -1;
    int id = d4est_solver_schwarz_is_element_in_subdomain
             (
              sub_data,
              e_p[i]->tree,
              e_p[i]->tree_quadid
             );

    if (id >= 0){
      zero_and_skip_p[i] = 0;
      e_p[i]->nodal_stride = sub_data->element_metadata[id].nodal_stride;
    }
    else {
      zero_and_skip_p[i] = 1;
      e_p[i]->nodal_stride = -1;
      skip_p_sum++;
    }
  }

  e_m[0]->boundary_quad_vector_stride[mortar_data->f_m] = mortar_data->boundary_quad_vector_stride;
  e_m[0]->mortar_quad_vector_stride[mortar_data->f_m] = mortar_data->mortar_quad_vector_stride;
  e_m[0]->mortar_quad_scalar_stride[mortar_data->f_m] = mortar_data->mortar_quad_scalar_stride;
  e_m[0]->mortar_quad_matrix_stride[mortar_data->f_m] = mortar_data->mortar_quad_matrix_stride;

  return (skip_p_sum == mortar_data->faces_p);
}

static void
d4est_solver_schwarz_laplacian_compute_dudr
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_solver_schwarz_subdomain_metadata_t* sub_data,
 d4est_mesh_data_t* d4est_factors,
 double * restrict  restrict_transpose_u_over_subdomain,
 double* dudr_over_subdomain [(P4EST_DIM)]
){
 int stride = 0;
 for (int i = 0; i < sub_data->num_elements; i++){
    d4est_element_data_t* ed = d4est_element_data_get_ptr(
                                                          p4est,
                                                          sub_data->element_metadata[i].tree,
                                                          sub_data->element_metadata[i].tree_quadid
                                                         );

   int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM),ed->deg);
   for (int d = 0; d < (P4EST_DIM); d++){
     d4est_operators_apply_dij(
                               d4est_ops,
                               &restrict_transpose_u_over_subdomain[stride],
                               (P4EST_DIM),
                               ed->deg,
                               d,
                               &dudr_over_subdomain[d][stride]
                              );
   }
   stride += volume_nodes_lobatto;
 }
}

static void
d4est_solver_schwarz_laplacian_compute_stiffness
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_subdomain_metadata_t* sub_data,
 double* dudr_over_subdomain [(P4EST_DIM)],
 double * restrict  restrict_transpose_u_over_subdomain,
 double * restrict  Au
)
{
 int stride = 0;
 for (int i = 0; i < sub_data->num_elements; i++){

    d4est_element_data_t* ed = d4est_element_data_get_ptr(
                                                          p4est,
                                                          sub_data->element_metadata[i].tree,
                                                          sub_data->element_metadata[i].tree_quadid
                                                         );

   
   int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM),ed->deg);
   d4est_laplacian_apply_stiffness_matrix_on_element
     (
      d4est_ops,
      d4est_geom,
      d4est_quad,
      d4est_factors,
      ed,
      &restrict_transpose_u_over_subdomain[stride],
      &Au[stride]
     );
   stride += volume_nodes_lobatto;
 }
}

void
d4est_solver_schwarz_laplacian_apply_over_subdomain
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* d4est_ghost,
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_laplacian_flux_data_t* flux_fcn_data,
 double* u_restricted_field_over_subdomain,
 double* Au_restricted_field_over_subdomain,
 int subdomain
)
{
  D4EST_ASSERT(p4est->mpisize == 1);
  d4est_solver_schwarz_subdomain_metadata_t* sub_data = &schwarz_data->subdomain_metadata[subdomain];

  double* Au = P4EST_ALLOC(double, sub_data->nodal_size);
  double* transpose_restricted_u = P4EST_ALLOC(double, sub_data->nodal_size);
  double* dudr [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(dudr, sub_data->nodal_size);
  
  d4est_solver_schwarz_apply_restrict_transpose_to_restricted_field_over_subdomain
    (
     schwarz_data,
     schwarz_ops,
     u_restricted_field_over_subdomain,
     transpose_restricted_u,
     subdomain
    );

  /* restrict transpose subdomain */
  d4est_element_data_t** e_m = P4EST_ALLOC(d4est_element_data_t*, (P4EST_HALF));
  d4est_element_data_t** e_p = P4EST_ALLOC(d4est_element_data_t*, (P4EST_HALF));

  for (int i = 0; i < (P4EST_HALF); i++){
    e_m[i] = P4EST_ALLOC(d4est_element_data_t, 1);
    e_p[i] = P4EST_ALLOC(d4est_element_data_t, 1);
  }
  
  /* Compute dudr */  
  d4est_solver_schwarz_laplacian_compute_dudr
    (
     p4est,
     d4est_ops,
     sub_data,
     d4est_factors,
     transpose_restricted_u,
     dudr
    );

  /* Compute stiffness */
  d4est_solver_schwarz_laplacian_compute_stiffness
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     sub_data,
     dudr,
     transpose_restricted_u,
     Au
    );

  
  /* Compute mortar integrals */
  flux_fcn_data->d4est_ghost_data = NULL;
  flux_fcn_data->d4est_ghost = d4est_ghost;
  flux_fcn_data->p4est = p4est;
  flux_fcn_data->which_field = 0;
  flux_fcn_data->local_nodes = sub_data->nodal_size;
  flux_fcn_data->u = transpose_restricted_u;
  flux_fcn_data->Au = Au;
  for (int i = 0; i < (P4EST_DIM); i++){
    flux_fcn_data->dudr_local[i] = dudr[i];
    flux_fcn_data->dudr_ghost[i] = NULL;
  }

  int* mortar_sides_done = P4EST_ALLOC_ZERO(int, sub_data->num_elements*(P4EST_FACES));
  /* int num_mortar_sides_done = 0; */
  
  for (int i = 0; i < sub_data->num_elements; i++){
    
    d4est_element_data_t* ed = d4est_element_data_get_ptr(
                                                          p4est,
                                                          sub_data->element_metadata[i].tree,
                                                          sub_data->element_metadata[i].tree_quadid
                                                         );
    
    for (int f = 0; f < (P4EST_FACES); f++){

      int mortar_side_id_m = ed->face_belongs_to_which_mortar[f];
      d4est_mortar_side_data_t* mortar_side_data = &d4est_factors->mortar_side_data[mortar_side_id_m];
                             
      
      int compute = 1;
      /* check if we have already computed this mortar */;
      if (mortar_sides_done[i*(P4EST_FACES) + f] == 1){
        compute = 0;
      }
        
      /* for (int c = 0; c < num_mortar_sides_done; i++){ */
      /*   if (mortar_sides_done[c] == mortar_side_data->mortar_side_id){ */
      /*     compute = 0; */
      /*     break; */
      /*   } */
      /* } */

      int zero_and_skip_m [] = {0,0,0,0};
      int zero_and_skip_p [] = {0,0,0,0};
      int e_m_is_ghost [] = {0,0,0,0};
      int element_m_sub_id [] = {-1,-1,-1,-1};
      
      if (compute){
        int skip_if_overlap_less_than_element_nodal_size
          = convert_mortar_side_data_to_elements
          (
           p4est,
           sub_data,
           mortar_side_data,
           e_m,
           e_p,
           zero_and_skip_m,
           zero_and_skip_p,
           element_m_sub_id
          );
        if (mortar_side_data->faces_p != 0 &&
            skip_if_overlap_less_than_element_nodal_size == 1
            && schwarz_data->num_nodes_overlap < ed->deg + 1){
          compute = 0;
        }            
      }

      /* printf("i e_m[0]->deg, f, mortar = %d %d %d %d\n", i, e_m[0]->deg, f, mortar_side_id_m); */
      
      /* optimization, don't compute if we skip all p elements and overlap size is not an entire element, because then none of these sides contribute anyway (e.g. flux is zero) */
  
      if (compute){
          d4est_laplacian_flux_schwarz
          (
           p4est,
           e_m,
           mortar_side_data->faces_m,
           mortar_side_data->f_m,
           mortar_side_data->mortar_side_id,
           e_p,
           mortar_side_data->faces_p,
           mortar_side_data->f_p,
           -1,
           e_m_is_ghost,
           zero_and_skip_m,
           zero_and_skip_p,
           mortar_side_data->orientation,
           d4est_ops,
           d4est_geom,
           d4est_quad,
           d4est_factors,
           flux_fcn_data
          );
          for (int m = 0; m < mortar_side_data->faces_m; m++){
            if (!zero_and_skip_m[m]){
              mortar_sides_done[element_m_sub_id[m]*(P4EST_FACES)
                                +  mortar_side_data->f_m] = 1;
            }
          }
              
        /* mortar_sides_done[num_mortar_sides_done] = mortar_side_data->mortar_side_id; */
        /* num_mortar_sides_done += 1; */
      }
    }   
  }

  d4est_solver_schwarz_convert_field_over_subdomain_to_restricted_field_over_subdomain
    (
     schwarz_data,
     schwarz_ops,
     Au,
     Au_restricted_field_over_subdomain,
     subdomain
    );
  
 
  for (int i = 0; i < (P4EST_HALF); i++){
    P4EST_FREE(e_m[i]);
    P4EST_FREE(e_p[i]);
  }
  P4EST_FREE(e_m);
  P4EST_FREE(e_p);

  P4EST_FREE(mortar_sides_done);
  P4EST_FREE(Au);
  P4EST_FREE(transpose_restricted_u);
  D4EST_FREE_DIM_VEC(dudr);
}
