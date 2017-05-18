#include <pXest.h>
#include <curved_element_data.h>
#include <problem_data.h>
#include <dgmath.h>
#include <linalg.h>
#include <curved_poisson_operator_primal.h>
#include <curved_poisson_debug_vecs.h>
#include <grid_functions.h>
#include <util.h>

typedef struct {

  d4est_geometry_t* d4est_geom;
  dgmath_jit_dbase_t* dgmath_jit_dbase;
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
  curved_element_data_t* elem_data = (curved_element_data_t *) q->p.user_data;

  curved_poisson_operator_primal_user_data_t* curved_poisson_operator_primal_user_data = (curved_poisson_operator_primal_user_data_t*) user_data;
  problem_data_t* problem_data = (problem_data_t*) curved_poisson_operator_primal_user_data->problem_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) curved_poisson_operator_primal_user_data->dgmath_jit_dbase;
  
  int dim = (P4EST_DIM);
  int deg = elem_data->deg;
  int volume_nodes_Lobatto = dgmath_get_nodes(dim,deg);
  /* int face_nodes_Lobatto = dgmath_get_nodes(dim-1,deg); */
  /* int volume_nodes_Gauss = dgmath_get_nodes(dim, elem_data->deg_integ); */
  
  elem_data->Au_elem = &(problem_data->Au[elem_data->nodal_stride]);
  
  linalg_copy_1st_to_2nd
    (
     &(problem_data->u[elem_data->nodal_stride]),
     &(elem_data->u_storage)[0],
     volume_nodes_Lobatto
    );

  linalg_fill_vec
    (
     elem_data->Au_elem,
     0.0,
     volume_nodes_Lobatto
    );
  
  for (int i = 0; i < (P4EST_DIM); i++){
    dgmath_apply_Dij(dgmath_jit_dbase, &(problem_data->u[elem_data->nodal_stride]), dim, elem_data->deg, i, &elem_data->dudr_elem[i][0]);
  }

}


void curved_poisson_operator_primal_compute_stiffmatrixterm
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t *q = info->quad;
  curved_element_data_t* element_data = q->p.user_data;

  curved_poisson_operator_primal_user_data_t* curved_poisson_operator_primal_user_data =  user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = curved_poisson_operator_primal_user_data->dgmath_jit_dbase;
  d4est_geometry_t* d4est_geom = curved_poisson_operator_primal_user_data->d4est_geom;

  
  int dim = (P4EST_DIM);
  int volume_nodes_lobatto = dgmath_get_nodes(dim,element_data->deg);

  double* stiff_u = P4EST_ALLOC(double, volume_nodes_lobatto);

  curved_element_data_apply_curved_stiffness_matrix
    (
     dgmath_jit_dbase,
     curved_poisson_operator_primal_user_data->d4est_geom,
     element_data,
     d4est_geom->geom_quad_type,
     &element_data->u_storage[0],
     stiff_u
    );

#ifdef NASTY_DEBUG
  printf("Stiffness Matrix Element id %d\n", element_data->id);
  DEBUG_PRINT_ARR_DBL_SUM(stiff_u, volume_nodes_lobatto);
#endif
  
  
  for (int i = 0; i < volume_nodes_lobatto; i++){
    element_data->Au_elem[i] += stiff_u[i];
  }

 
  

  P4EST_FREE(stiff_u);
}



void
curved_poisson_operator_primal_apply_aij
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 problem_data_t* prob_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* geom
)
{
  curved_poisson_operator_primal_user_data_t curved_poisson_operator_primal_user_data;
  curved_poisson_operator_primal_user_data.dgmath_jit_dbase = dgmath_jit_dbase;
  curved_poisson_operator_primal_user_data.problem_data = prob_vecs;
  curved_poisson_operator_primal_user_data.d4est_geom = geom;
#ifndef NDEBUG
  curved_poisson_operator_primal_user_data.debug_vecs = NULL;
#endif
  
  curved_compute_flux_user_data_t curved_compute_flux_user_data;
  curved_compute_flux_user_data.dgmath_jit_dbase = dgmath_jit_dbase;
  curved_compute_flux_user_data.geom = geom;
  
  void* tmpptr = p4est->user_pointer;
  
  p4est_iterate(p4est,
		NULL,
		(void *) &curved_poisson_operator_primal_user_data,
		curved_poisson_operator_primal_init_vecs,
		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);

  p4est_ghost_exchange_data(p4est,ghost,ghost_data);
 
  curved_compute_flux_user_data.flux_fcn_ptrs = &prob_vecs->curved_scalar_flux_fcn_data;
  p4est->user_pointer = &curved_compute_flux_user_data;

  p4est_iterate(p4est,
		ghost,
		(void*) ghost_data,
		NULL,
		curved_compute_flux_on_local_elements,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL);
 
  
  p4est_iterate (p4est,
  		 NULL,
                 (void*)&curved_poisson_operator_primal_user_data,
  		 curved_poisson_operator_primal_compute_stiffmatrixterm,
  		 NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  
  p4est->user_pointer = tmpptr;
}
