#include <pXest.h>
#include <d4est_element_data.h>
#include <d4est_elliptic_data.h>
#include <d4est_operators.h>
#include <d4est_linalg.h>
#include <d4est_poisson.h>
#include <d4est_poisson_flux.h>
#include <d4est_xyz_functions.h>
#include <d4est_quadrature.h>
#include <d4est_mortars.h>
#include <util.h>

void
d4est_poisson_apply_aij
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_poisson_flux_data_t* flux_fcn_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad
)
{
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
  
        d4est_linalg_copy_1st_to_2nd
          (
           &(prob_vecs->u[ed->nodal_stride]),
           &(ed->u_elem)[0],
           volume_nodes_lobatto
          );

        for (int i = 0; i < (P4EST_DIM); i++){
          d4est_operators_apply_dij(d4est_ops, &(prob_vecs->u[ed->nodal_stride]), (P4EST_DIM), ed->deg, i, &ed->dudr_elem[i][0]);
        }

        double* stiff_u = P4EST_ALLOC(double, volume_nodes_lobatto);
        d4est_quadrature_volume_t mesh_vol = {.dq = ed->dq,
                                              .tree = ed->tree,
                                              .element_id = ed->id
                                             };

        for (int i = 0; i < (P4EST_DIM); i++)
          mesh_vol.q[i] = ed->q[i];
  
        d4est_quadrature_apply_stiffness_matrix
          (
           d4est_ops,
           d4est_quad,
           d4est_geom,
           &mesh_vol,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &ed->u_elem[0],
           ed->deg,
           ed->J_quad,
           ed->rst_xyz_quad,
           ed->deg_quad,
           ed->Au_elem
          );
  
        P4EST_FREE(stiff_u);
      }

    }


  d4est_mortar_fcn_ptrs_t flux_fcns = d4est_poisson_flux_fetch_fcns(flux_fcn_data);
  d4est_mortar_compute_flux_on_local_elements
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
    
}
