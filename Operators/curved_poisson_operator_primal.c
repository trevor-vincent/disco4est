#include <pXest.h>
#include <d4est_element_data.h>
#include <problem_data.h>
#include <d4est_operators.h>
#include <d4est_linalg.h>
#include <curved_poisson_operator_primal.h>
#include <curved_poisson_debug_vecs.h>
#include <grid_functions.h>
#include <d4est_quadrature.h>
#include <util.h>

/* #define NASTY_DEBUG */

typedef struct {

  d4est_geometry_t* d4est_geom;
  d4est_quadrature_t* d4est_quad;
  d4est_operators_t* d4est_ops;
  problem_data_t* problem_data;
#ifndef NDEBUG
  curved_poisson_debug_vecs_t* debug_vecs;
#endif
} curved_poisson_operator_primal_user_data_t;

void curved_poisson_operator_primal_init_vecs
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  d4est_element_data_t* elem_data = (d4est_element_data_t *) q->p.user_data;

  curved_poisson_operator_primal_user_data_t* curved_poisson_operator_primal_user_data = (curved_poisson_operator_primal_user_data_t*) user_data;
  problem_data_t* problem_data = (problem_data_t*) curved_poisson_operator_primal_user_data->problem_data;
  d4est_operators_t* d4est_ops = (d4est_operators_t*) curved_poisson_operator_primal_user_data->d4est_ops;
  
  int dim = (P4EST_DIM);
  int deg = elem_data->deg;
  int volume_nodes_Lobatto = d4est_operators_get_nodes(dim,deg);
  /* int face_nodes_Lobatto = d4est_operators_get_nodes(dim-1,deg); */
  /* int volume_nodes_Gauss = d4est_operators_get_nodes(dim, elem_data->deg_quad); */
  
  elem_data->Au_elem = &(problem_data->Au[elem_data->nodal_stride]);
  
  d4est_linalg_copy_1st_to_2nd
    (
     &(problem_data->u[elem_data->nodal_stride]),
     &(elem_data->u_elem)[0],
     volume_nodes_Lobatto
    );

  d4est_linalg_fill_vec
    (
     elem_data->Au_elem,
     0.0,
     volume_nodes_Lobatto
    );
  
  for (int i = 0; i < (P4EST_DIM); i++){
    d4est_operators_apply_Dij(d4est_ops, &(problem_data->u[elem_data->nodal_stride]), dim, elem_data->deg, i, &elem_data->dudr_elem[i][0]);
  }

}


void curved_poisson_operator_primal_compute_stiffmatrixterm
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t *q = info->quad;
  d4est_element_data_t* element_data = q->p.user_data;

  curved_poisson_operator_primal_user_data_t* curved_poisson_operator_primal_user_data =  user_data;
  d4est_operators_t* d4est_ops = curved_poisson_operator_primal_user_data->d4est_ops;
  d4est_geometry_t* d4est_geom = curved_poisson_operator_primal_user_data->d4est_geom;
  d4est_quadrature_t* d4est_quad =
    curved_poisson_operator_primal_user_data->d4est_quad;
  
  
  int dim = (P4EST_DIM);
  int volume_nodes_lobatto = d4est_operators_get_nodes(dim,element_data->deg);

  double* stiff_u = P4EST_ALLOC(double, volume_nodes_lobatto);
  d4est_quadrature_volume_t mesh_vol = {.dq = element_data->dq,
                                        .tree = element_data->tree,
                                        .element_id = element_data->id
                                       };

  for (int i = 0; i < (P4EST_DIM); i++)
    mesh_vol.q[i] = element_data->q[i];
  
  d4est_quadrature_apply_stiffness_matrix
    (
     d4est_ops,
     d4est_quad,
     d4est_geom,
     &mesh_vol,
     QUAD_VOLUME,
     QUAD_DRDX_JAC_DRDX_TIMES_POLY_INTEGRAND,
     &element_data->u_elem[0],
     element_data->deg,
     element_data->J_quad,
     element_data->rst_xyz_quad,
     element_data->deg_quad,
     stiff_u
    );
  
  for (int i = 0; i < volume_nodes_lobatto; i++){
    element_data->Au_elem[i] += stiff_u[i];
  }
  
#ifdef NASTY_DEBUG
  printf("Stiffness Matrix Element id %d\n", element_data->id);
  DEBUG_PRINT_ARR_DBL_SUM(stiff_u, volume_nodes_lobatto);
  DEBUG_PRINT_ARR_DBL_SUM(element_data->Au_elem, volume_nodes_lobatto);
#endif
  
  P4EST_FREE(stiff_u);
}



void
curved_poisson_operator_primal_apply_aij
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad
)
{
  curved_poisson_operator_primal_user_data_t curved_poisson_operator_primal_user_data;
  curved_poisson_operator_primal_user_data.d4est_ops = d4est_ops;
  curved_poisson_operator_primal_user_data.problem_data = prob_vecs;
  curved_poisson_operator_primal_user_data.d4est_geom = d4est_geom;
  curved_poisson_operator_primal_user_data.d4est_quad = d4est_quad;
#ifndef NDEBUG
  curved_poisson_operator_primal_user_data.debug_vecs = NULL;
#endif
  
  p4est_iterate(p4est,
		NULL,
		(void *) &curved_poisson_operator_primal_user_data,
		curved_poisson_operator_primal_init_vecs,
		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);

 
  curved_compute_flux_on_local_elements
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &prob_vecs->curved_scalar_flux_fcn_data,
     EXCHANGE_GHOST_DATA
    );
     
  p4est_iterate (p4est,
  		 NULL,
                 (void*)&curved_poisson_operator_primal_user_data,
  		 curved_poisson_operator_primal_compute_stiffmatrixterm,
  		 NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  


}
