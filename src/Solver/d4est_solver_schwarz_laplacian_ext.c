#include <pXest.h>
#include <d4est_solver_schwarz_laplacian_ext.h>
#include <d4est_solver_schwarz_helpers.h>
#include <d4est_solver_schwarz_geometric_data.h>
#include <d4est_laplacian.h>
#include <d4est_laplacian_flux.h>

static void
d4est_solver_schwarz_laplacian_ext_compute_dudr
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
   d4est_solver_schwarz_element_metadata_t* schwarz_ed = &sub_data->element_metadata[i];

   int volume_nodes_lobatto
     = d4est_lgl_get_nodes(
                           (P4EST_DIM),
                           schwarz_ed->deg
                          );
   
   for (int d = 0; d < (P4EST_DIM); d++){
     d4est_operators_apply_dij
       (
        d4est_ops,
        &restrict_transpose_u_over_subdomain[stride],
        (P4EST_DIM),
        schwarz_ed->deg,
        d,
        &dudr_over_subdomain[d][stride]
       );
   }
   stride += volume_nodes_lobatto;
 }
}

static void
d4est_solver_schwarz_laplacian_ext_compute_stiffness
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_ghost_t* d4est_ghost,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_subdomain_metadata_t* sub_data,
 d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data,
 double * restrict  restrict_transpose_u_over_subdomain,
 double * restrict  Au
)
{
  /* printf("\n\n mpirank %d, subdomain id %d\n", p4est->mpirank, sub_data->subdomain_id); */
 int stride = 0;
 for (int el = 0; el < sub_data->num_elements; el++){
   /* printf("\n element %d\n", i); */
   d4est_solver_schwarz_element_metadata_t* schwarz_ed = &sub_data->element_metadata[el];

   d4est_element_data_t* mesh_ed = NULL;
   double* J_quad = NULL;

   int is_ghost = schwarz_ed->mpirank != p4est->mpirank;
   double* rst_xyz_quad [P4EST_DIM][P4EST_DIM];

   if (!is_ghost){
     mesh_ed = d4est_element_data_get_ptr
               (
                p4est,
                schwarz_ed->tree,
                schwarz_ed->tree_quadid
               );

     J_quad = d4est_mesh_get_jacobian_on_quadrature_points
              (
               d4est_factors,
               mesh_ed
              );

     for (int i = 0; i < (P4EST_DIM); i++){
       for (int j = 0; j < (P4EST_DIM); j++){
         rst_xyz_quad[i][j] = &d4est_factors->rst_xyz_quad
                              [(i*(P4EST_DIM) + j)*d4est_factors->local_sizes.local_nodes_quad
                               + mesh_ed->quad_stride];
       }
     }     
   }
   else {
     mesh_ed = &d4est_ghost->ghost_elements[schwarz_ed->id];
     int ghost_quad_stride = schwarz_geometric_data->volume_quad_strides_per_ghost[schwarz_ed->id];
     int total_ghost_size = schwarz_geometric_data->total_ghost_volume_quad_size;
     J_quad = &schwarz_geometric_data->J_quad_ghost[ghost_quad_stride];
     for (int i = 0; i < (P4EST_DIM); i++){
       for (int j = 0; j < (P4EST_DIM); j++){
         rst_xyz_quad[i][j] = &schwarz_geometric_data->rst_xyz_quad_ghost
                              [(i*(P4EST_DIM) + j)*total_ghost_size
                               + ghost_quad_stride];
       }
     }     
     /* printf("element is ghost\n"); */
   }
   int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM),
                                                  schwarz_ed->deg);



  d4est_quadrature_volume_t mesh_vol = {.dq = mesh_ed->dq,
                                        .tree = mesh_ed->tree,
                                        .q[0] = mesh_ed->q[0],
                                        .q[1] = mesh_ed->q[1],
#if (P4EST_DIM)==3                                              
                                        .q[2] = mesh_ed->q[2],
#endif
                                        .element_id = mesh_ed->id
                                       };



  
  d4est_quadrature_apply_stiffness_matrix
    (
     d4est_ops,
     d4est_quad,
     d4est_geom,
     &mesh_vol,
     QUAD_OBJECT_VOLUME,
     QUAD_INTEGRAND_UNKNOWN,
     &restrict_transpose_u_over_subdomain[stride],
     mesh_ed->deg,
     J_quad,
     rst_xyz_quad,
     mesh_ed->deg_quad,
     &Au[stride]
    );

  /* double Au_sum = 0; */
  /* double J_sum = 0; */
  /* double rst_xyz_sum = 0; */

  /* for (int k = 0; k < volume_nodes_lobatto; k++){ */
  /*   J_sum += J_quad[k]; */
  /*   Au_sum += Au[stride + k]; */
  /*    for (int i = 0; i < (P4EST_DIM); i++){ */
  /*      for (int j = 0; j < (P4EST_DIM); j++){ */
  /*        rst_xyz_sum += rst_xyz_quad[i][j][k]; */
  /*      } */
  /*    } */
  /* } */
  
  /* if (!is_ghost){ */
  /* printf("core_tree, elem,  J, rstxyz, Au  %d %d %d %.15f %.15f %.15f\n", sub_data->core_tree, el, is_ghost, J_quad[1], rst_xyz_quad[1][0][0], Au[stride]); */
  /* } */
  /* else { */
     /* int ghost_quad_stride = schwarz_geometric_data->volume_quad_strides_per_ghost[schwarz_ed->id]; */
     /* printf("core_tree, elem,  gstride, J, rstxyz, Au  %d %d %d %.15f %.15f %.15f\n", sub_data->core_tree, el, ghost_quad_stride, J_quad[1], rst_xyz_quad[1][0][0], Au[stride]); */
  /* } */

  stride += volume_nodes_lobatto;
 }
}

void
d4est_solver_schwarz_laplacian_ext_apply_over_subdomain
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* d4est_ghost,
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data,
 d4est_laplacian_flux_data_t* flux_fcn_data,
 double* u_restricted_field_over_subdomain,
 double* Au_restricted_field_over_subdomain,
 int subdomain
)
{
  /* D4EST_ASSERT(p4est->mpisize == 1); */
  d4est_solver_schwarz_subdomain_metadata_t* sub_data = &schwarz_data->subdomain_metadata[subdomain];

  double* Au = P4EST_ALLOC_ZERO(double, sub_data->nodal_size);
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

  /* Compute dudr */ 
  d4est_solver_schwarz_laplacian_ext_compute_dudr
    (
     p4est,
     d4est_ops,
     sub_data,
     d4est_factors,
     transpose_restricted_u,
     dudr
    );

  /* Compute stiffness */
  d4est_solver_schwarz_laplacian_ext_compute_stiffness
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_ghost,
     d4est_factors,
     sub_data,
     schwarz_geometric_data,
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
  flux_fcn_data->subdomain = subdomain;
  
  for (int i = 0; i < (P4EST_DIM); i++){
    flux_fcn_data->dudr_local[i] = dudr[i];
    flux_fcn_data->dudr_ghost[i] = NULL;
  }

   int subdomain_mortar_stride =
     schwarz_geometric_data->mortar_strides_per_subdomain[subdomain];
   int subdomain_num_mortars =
     schwarz_geometric_data->num_of_mortars_per_subdomain[subdomain];


   for (int i = 0; i < subdomain_num_mortars; i++){

    d4est_mortar_side_data_t* mortar_side_data =
      schwarz_geometric_data->subdomain_mortars[subdomain_mortar_stride + i];

    
    if(!mortar_side_data->is_ghost){
      flux_fcn_data->using_provided_mesh_data = 0;
    }
    else{

      int scalar_stride = mortar_side_data->mortar_quad_stride;
      int vector_stride = mortar_side_data->mortar_quad_stride*(P4EST_DIM);
      int matrix_stride = mortar_side_data->mortar_quad_stride*(P4EST_DIM)*(P4EST_DIM);
      int boundary_vector_stride = mortar_side_data->boundary_quad_stride*(P4EST_DIM);
      
      flux_fcn_data->using_provided_mesh_data = 1;
      flux_fcn_data->hm_mortar_quad = &schwarz_geometric_data->hm_mortar_quad[scalar_stride];
      flux_fcn_data->hp_mortar_quad = &schwarz_geometric_data->hp_mortar_quad[scalar_stride];
      flux_fcn_data->sj_m_mortar_quad = &schwarz_geometric_data->sj_m_mortar_quad[scalar_stride];
      flux_fcn_data->n_m_mortar_quad = &schwarz_geometric_data->n_m_mortar_quad[vector_stride];
      flux_fcn_data->xyz_m_mortar_quad = &schwarz_geometric_data->xyz_on_f_m_quad[boundary_vector_stride];
      flux_fcn_data->xyz_m_mortar_lobatto = &schwarz_geometric_data->xyz_on_f_m_lobatto[boundary_vector_stride];
      flux_fcn_data->drst_dxyz_m_mortar_quad = &schwarz_geometric_data->drst_dxyz_m_mortar_quad[matrix_stride];
      flux_fcn_data->drst_dxyz_p_mortar_quad_porder =
        &schwarz_geometric_data->drst_dxyz_p_mortar_quad_porder[matrix_stride];
    }
    
    int skip_stride = (P4EST_HALF)*(subdomain_mortar_stride + i);
    int* zero_and_skip_m
      = &schwarz_geometric_data->zero_and_skip_m[skip_stride];
    int* zero_and_skip_p
      = &schwarz_geometric_data->zero_and_skip_p[skip_stride];
    int* nodal_stride_m
      = &schwarz_geometric_data->nodal_stride_m[skip_stride];
    int* nodal_stride_p
      = &schwarz_geometric_data->nodal_stride_p[skip_stride];
    
    int faces_mortar = (mortar_side_data->faces_p >
                        mortar_side_data->faces_m) ?
                       mortar_side_data->faces_p :
                       mortar_side_data->faces_m;

      d4est_element_data_t* e_m [P4EST_HALF];
      d4est_element_data_t* e_p [P4EST_HALF];
      for (int e = 0; e < mortar_side_data->faces_m; e++){
        e_m[e] = &mortar_side_data->e_m[e];
        e_m[e]->nodal_stride = nodal_stride_m[e];
        e_m[e]->mpirank = p4est->mpirank;
      }
      for (int e = 0; e < mortar_side_data->faces_p; e++){
        e_p[e] = &mortar_side_data->e_p[e];
        e_p[e]->nodal_stride = nodal_stride_p[e];
        e_p[e]->mpirank = p4est->mpirank;
      }
      
      int e_m_is_ghost [] = {0,0,0,0};
      /* if (!do_not_compute){ */
      d4est_laplacian_flux_schwarz
        (
         p4est,
         &e_m[0],
         mortar_side_data->faces_m,
         mortar_side_data->f_m,
         mortar_side_data->mortar_side_id,
         &e_p[0],
         mortar_side_data->faces_p,
         mortar_side_data->f_p,
         -1,
         &e_m_is_ghost[0],
         zero_and_skip_m,
         zero_and_skip_p,
         mortar_side_data->orientation,
         d4est_ops,
         d4est_geom,
         d4est_quad,
         d4est_factors,
         flux_fcn_data
        );
      /* } */
  }

  
  d4est_solver_schwarz_convert_field_over_subdomain_to_restricted_field_over_subdomain
    (
     schwarz_data,
     schwarz_ops,
     Au,
     Au_restricted_field_over_subdomain,
     subdomain
    );


  /* for (int i = 0; i < sub_data->nodal_size; i++){ */
  /*   if (Au[i] != Au[i]){ */
  /*     printf("mpirank = %d, subdomain id = %d, Au[i] = %.15f\n", p4est->mpirank, subdomain, Au[i]); */
  /*     D4EST_ABORT(""); */
  /*   } */
  /* } */
  
  /* if(p4est->mpirank == 3 && subdomain == 1){ */
  /*   printf("p4est->mpirank = %d, subdomain id = %d\n", p4est->mpirank, subdomain); */
    /* DEBUG_PRINT_ARR_DBL(Au_restricted_field_over_subdomain, */
                        /* sub_data->restricted_nodal_size); */
  /*   D4EST_ABORT(""); */
  /* } */

  /* P4EST_FREE(mortar_sides_done); */
  P4EST_FREE(Au);
  P4EST_FREE(transpose_restricted_u);
  D4EST_FREE_DIM_VEC(dudr);
}
