#include <d4est_ip_energy_norm.h>
#include <pXest.h>
#include <d4est_util.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_element_data.h>
#include <d4est_linalg.h>
#include <d4est_mesh.h>
#include <d4est_estimator_bi.h>
#include <d4est_gradient.h>
#include <d4est_poisson_flux.h>
#include <d4est_mortars.h>
#include <d4est_poisson.h>

static void
d4est_ip_energy_norm_boundary
(
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_poisson_flux_boundary_data_t* boundary_data,
 void* bc_params,
 void* params
)
{
  d4est_ip_energy_norm_data_t* data = params;
  d4est_quadrature_mortar_t* face_object = boundary_data->face_object;
  int face_nodes_m_lobatto = boundary_data->face_nodes_m_lobatto;
  int face_nodes_m_quad = boundary_data->face_nodes_m_quad;
  int deg_mortar_quad = boundary_data->deg_mortar_quad;

  double* u_m_on_f_m_quad = boundary_data->u_m_on_f_m_quad;
  double* sj_on_f_m_quad = boundary_data->sj_on_f_m_quad;
  double* j_div_sj_quad = boundary_data->j_div_sj_quad;
  double* xyz_on_f_m_lobatto [(P4EST_DIM)]; 
  double* drst_dxyz_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* dudx_m_on_f_m_quad [(P4EST_DIM)];  
  double* n_on_f_m_quad [(P4EST_DIM)];
  D4EST_COPY_DBYD_MAT(boundary_data->drst_dxyz_quad, drst_dxyz_quad);
  D4EST_COPY_DIM_VEC(boundary_data->dudx_m_on_f_m_quad, dudx_m_on_f_m_quad);
  D4EST_COPY_DIM_VEC(boundary_data->n_on_f_m_quad, n_on_f_m_quad);
  D4EST_COPY_DIM_VEC(boundary_data->xyz_on_f_m_lobatto, xyz_on_f_m_lobatto); 
  double* ip_energy_norm_prefactor = P4EST_ALLOC(double, face_nodes_m_quad);
  double* ip_energy_norm = P4EST_ALLOC(double, face_nodes_m_quad);

  double* h_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  d4est_poisson_flux_sipg_calculate_h
    (
     &e_m,
     f_m,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     data->sipg_flux_h,
     j_div_sj_quad,
     h_quad,
     1,
     1,
     &face_nodes_m_quad
    );

  for (int i = 0; i < face_nodes_m_quad; i++){
    ip_energy_norm_prefactor[i] = data->u_penalty_fcn
                       (
                        e_m->deg,
                        h_quad[i],
                        e_m->deg,
                        h_quad[i],
                        data->penalty_prefactor
                       );
  }  
  P4EST_FREE(h_quad);

    for(int i = 0; i < face_nodes_m_quad; i++){
      ip_energy_norm[i] = 0.;
      for (int dim = 0; dim < (P4EST_DIM); dim++){
        ip_energy_norm[i] += n_on_f_m_quad[dim][i]*(u_m_on_f_m_quad[i])*n_on_f_m_quad[dim][i]*(u_m_on_f_m_quad[i]);
      }
    }

    data->ip_energy_norm_sqr_boundary_term
      += d4est_quadrature_innerproduct
      (
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &face_object,
       QUAD_OBJECT_MORTAR,
       QUAD_INTEGRAND_UNKNOWN,
       ip_energy_norm_prefactor,
       ip_energy_norm,
       sj_on_f_m_quad,
       deg_mortar_quad
      );
    

  P4EST_FREE(ip_energy_norm_prefactor);
  P4EST_FREE(ip_energy_norm);
}

static void
d4est_ip_energy_norm_interface
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
 d4est_poisson_flux_interface_data_t* interface_data,
 void* params
)
{
  d4est_quadrature_mortar_t* mortar_face_object = interface_data->mortar_face_object;
  
  int faces_mortar = interface_data->faces_mortar;
  int total_side_nodes_m_lobatto = interface_data->total_side_nodes_m_lobatto;
  int total_nodes_mortar_lobatto = interface_data->total_nodes_mortar_lobatto;
  int total_nodes_mortar_quad = interface_data->total_nodes_mortar_quad;
  
  double* u_m_on_f_m_mortar_quad = interface_data->u_m_on_f_m_mortar_quad;
  double* sj_on_f_m_mortar_quad = interface_data->sj_on_f_m_mortar_quad;
  double* j_div_sj_on_f_m_mortar_quad = interface_data->j_div_sj_on_f_m_mortar_quad;
  double* u_p_on_f_p_mortar_quad = interface_data->u_p_on_f_p_mortar_quad;
  double* j_div_sj_on_f_p_mortar_quad = interface_data->j_div_sj_on_f_p_mortar_quad;
  
  double* drst_dxyz_m_on_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_quad [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_quad [(P4EST_DIM)];
  double* n_on_f_m_mortar_quad [(P4EST_DIM)];

  D4EST_COPY_DBYD_MAT(interface_data->drst_dxyz_m_on_mortar_quad, drst_dxyz_m_on_mortar_quad);
  D4EST_COPY_DIM_VEC(interface_data->dudx_m_on_f_m_mortar_quad, dudx_m_on_f_m_mortar_quad);
  D4EST_COPY_DIM_VEC(interface_data->dudx_p_on_f_p_mortar_quad, dudx_p_on_f_p_mortar_quad);
  D4EST_COPY_DIM_VEC(interface_data->n_on_f_m_mortar_quad, n_on_f_m_mortar_quad);
  
  int* deg_mortar_quad = interface_data->deg_mortar_quad;
  int* nodes_mortar_quad = interface_data->nodes_mortar_quad;
  int* nodes_mortar_lobatto = interface_data->nodes_mortar_lobatto;
  int* deg_mortar_lobatto = interface_data->deg_mortar_lobatto;
  int* face_nodes_m_lobatto = interface_data->deg_mortar_lobatto;
  int* deg_m_lobatto = interface_data->deg_m_lobatto;
  int* deg_p_lobatto = interface_data->deg_p_lobatto;
  
  d4est_ip_energy_norm_data_t* data = params;
  

  double* ip_energy_norm_prefactor = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* ip_energy_norm = P4EST_ALLOC(double, total_nodes_mortar_quad);


  double* hm_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* hp_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);

  d4est_poisson_flux_sipg_calculate_h
    (
     e_m,
     f_m,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     data->sipg_flux_h,
     j_div_sj_on_f_m_mortar_quad,
     hm_mortar_quad,
    interface_data->faces_mortar,
     faces_m,
     nodes_mortar_quad
    );


  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
 d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);
  
    d4est_poisson_flux_sipg_calculate_h
    (
     &e_p_oriented[0],
     f_p,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     data->sipg_flux_h,
     j_div_sj_on_f_p_mortar_quad,
     hp_mortar_quad,
     interface_data->faces_mortar,
     faces_p,
     nodes_mortar_quad
    );
  
  int stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_quad[f]; k++){
      int ks = k + stride;
      ip_energy_norm_prefactor[ks] = data->u_penalty_fcn
                          (
                           (faces_m == faces_mortar) ? deg_m_lobatto[f] : deg_m_lobatto[0],
                           hm_mortar_quad[ks],
                           (faces_p == faces_mortar) ? deg_p_lobatto[f] : deg_p_lobatto[0],
                           hp_mortar_quad[ks],
                           data->penalty_prefactor
                          );    
      

    }
    stride += nodes_mortar_quad[f];
  }

  P4EST_FREE(hm_mortar_quad);
  P4EST_FREE(hp_mortar_quad);
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
      ip_energy_norm[ks] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        n_ks = n_on_f_m_mortar_quad[d][ks];        
        ip_energy_norm[ks] += pow(n_ks*u_m_on_f_m_mortar_quad[ks] - n_ks*u_p_on_f_p_mortar_quad[ks],2);
      }
    }
    stride += nodes_mortar_quad[f];
  }

  /* the contribution in every direction must be added up due to it being a vector norm */
  stride = 0;
  for (f = 0; f < faces_mortar; f++){
    for (int d = 0; d < (P4EST_DIM); d++){

      data->ip_energy_norm_sqr_interface_term +=
        d4est_quadrature_innerproduct
        (
         d4est_ops,
         d4est_geom,
         d4est_quad,
         mortar_face_object,
         QUAD_OBJECT_MORTAR,
         QUAD_INTEGRAND_UNKNOWN,
         &ip_energy_norm_prefactor[stride],
         &ip_energy_norm[stride],
         &sj_on_f_m_mortar_quad[stride],
         deg_mortar_quad[f]
        );
    }
    stride += nodes_mortar_quad[f];
  }

  P4EST_FREE(ip_energy_norm);
  P4EST_FREE(ip_energy_norm_prefactor);
}

static int
d4est_ip_energy_norm_get_deg_mortar_quad
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  return elem_data->deg_quad;
}

double
d4est_ip_energy_norm_compute
(
 p4est_t* p4est,
 double* u,
 d4est_ip_energy_norm_data_t* energy_norm_data,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors
)
{
  energy_norm_data->ip_energy_norm_sqr_volume_term = 0.;
  energy_norm_data->ip_energy_norm_sqr_interface_term = 0.;
  energy_norm_data->ip_energy_norm_sqr_boundary_term = 0.;

  
  int nodes = 0;
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
        nodes += volume_nodes_lobatto;
        
        d4est_util_copy_1st_to_2nd
          (
           &(u[ed->nodal_stride]),
           &(ed->u_elem)[0],
           volume_nodes_lobatto
          );

        d4est_quadrature_volume_t volume_object = {.dq = ed->dq,
                                              .tree = ed->tree,
                                              .q[0] = ed->q[0],
                                              .q[1] = ed->q[1],
#if (P4EST_DIM)==3                                              
                                              .q[2] = ed->q[2],
#endif
                                              .element_id = ed->id
                                             };

        double* J_quad = d4est_mesh_get_jacobian_on_quadrature_points(d4est_factors,
                                                                      ed);

        
        energy_norm_data->ip_energy_norm_sqr_volume_term += d4est_gradient_l2_norm
                                                           (
                                                            d4est_ops,
                                                            d4est_geom,
                                                            d4est_quad,
                                                            &volume_object,
                                                            QUAD_OBJECT_VOLUME,
                                                            QUAD_INTEGRAND_UNKNOWN,
                                                            &(u[ed->nodal_stride]),
                                                            ed->deg,
                                                            ed->deg_quad,
                                                            ed->rst_xyz_quad,
                                                            J_quad
                                                           );
    
      }

    }


  p4est_ghost_exchange_data(p4est,ghost,ghost_data);

  int ghost_nodes = d4est_mesh_get_ghost_nodes(ghost, ghost_data);
  double* dudr [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(dudr,
                                                  nodes
                                                  + ghost_nodes
                                                 );


  d4est_poisson_compute_dudr
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     dudr
    );
  
  
  d4est_poisson_flux_data_t flux_data;
  flux_data.interface_fcn = d4est_ip_energy_norm_interface;
  flux_data.boundary_fcn = d4est_ip_energy_norm_boundary;
  flux_data.bc_type = BC_NOT_SET;
  flux_data.flux_data = energy_norm_data;
  flux_data.get_deg_mortar_quad = d4est_ip_energy_norm_get_deg_mortar_quad;
  flux_data.get_deg_mortar_quad_ctx = NULL;
  
  d4est_mortars_fcn_ptrs_t flux_fcns = d4est_poisson_flux_fetch_fcns(&flux_data);
  
  d4est_mortars_compute_flux_on_local_elements
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &flux_fcns,
     DO_NOT_EXCHANGE_GHOST_DATA
    );


  double ip_energy_norm_sqr = 0.;
  ip_energy_norm_sqr += energy_norm_data->ip_energy_norm_sqr_volume_term;
  ip_energy_norm_sqr += energy_norm_data->ip_energy_norm_sqr_boundary_term;
  ip_energy_norm_sqr += energy_norm_data->ip_energy_norm_sqr_interface_term;
  return ip_energy_norm_sqr;

  D4EST_FREE_DIM_VEC(dudr);
  
}
