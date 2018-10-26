#include <pXest.h>
#include <d4est_util.h>
#include <d4est_ghost.h>
#include <d4est_operators.h>
#include <d4est_mesh.h>
#include <d4est_solver_schwarz.h>
#include <d4est_solver_schwarz_operators.h>
#include <d4est_solver_schwarz_helpers.h>
#include <d4est_field.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_linalg.h>

/** 
 * Iterate over a field and restrict it to each
 * subdomain
 * 
 * @param [in] schwarz_data 
 * @param [in] field 
 * @param [in] ghost_data 
 * @param [in] field_num 
 * @param [out] field_over_subdomains 
 */
void
d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomain
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_data_t* ghost_data,
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* field,
 int ghost_data_num_of_field,
 double* restricted_field_over_subdomain,
 int subdomain
)
{
  d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[subdomain];
  for (int j = 0; j < sub_data.num_elements; j++){
    d4est_solver_schwarz_element_data_t ed = sub_data.element_data[j];

    double* field_ed = NULL;      
    if (ed.mpirank == p4est->mpirank){
      field_ed = &field[d4est_factors->element_data[ed.id]->nodal_stride];
    }
    else {
      field_ed =
        d4est_ghost_data_get_field_on_element
        (
         d4est_factors->element_data[ed.id + p4est->local_num_quadrants],
         ghost_data_num_of_field,
         ghost_data
        );
    }
    d4est_solver_schwarz_operators_apply_schwarz_restrictor
      (
       schwarz_ops,
       field_ed,
       (P4EST_DIM),
       &(ed.faces[0]),
       ed.deg,
       schwarz_data->num_nodes_overlap,
       D4OPS_NO_TRANSPOSE,
       &restricted_field_over_subdomain[ed.restricted_nodal_stride]
      );
  }
}

void
d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomains
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_data_t* ghost_data,
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* field,
 int ghost_data_num_of_field,
 double* restricted_field_over_subdomains
)
{
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[i];
    d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomain
      (
       p4est,
       d4est_ops,
       d4est_factors,
       ghost_data,
       schwarz_data,
       schwarz_ops,
       field,
       ghost_data_num_of_field,
       &restricted_field_over_subdomains[sub_data.restricted_nodal_stride],
       i
      );
  }
}

void
d4est_solver_schwarz_convert_field_over_subdomain_to_restricted_field_over_subdomain
(
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* field_over_subdomains,
 double* restricted_field_over_subdomains,
 int subdomain
)
{
  d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[subdomain];
  for (int i = 0; i < sub_data.num_elements; i++){
    d4est_solver_schwarz_element_data_t ed = sub_data.element_data[i];
    d4est_solver_schwarz_operators_apply_schwarz_restrictor
      (
       schwarz_ops,
       &field_over_subdomains[ed.nodal_stride],
       (P4EST_DIM),
       &(ed.faces[0]),
       ed.deg,
       schwarz_data->num_nodes_overlap,
       D4OPS_NO_TRANSPOSE,
       &restricted_field_over_subdomains[ed.restricted_nodal_stride]      
      );
  }
}

void
d4est_solver_schwarz_convert_field_over_subdomains_to_restricted_field_over_subdomains
(
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* field_over_subdomains,
 double* restricted_field_over_subdomains
)
{
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[i];
    for (int j = 0; j < sub_data.num_elements; j++){
      d4est_solver_schwarz_element_data_t ed = sub_data.element_data[j];

      d4est_solver_schwarz_convert_field_over_subdomain_to_restricted_field_over_subdomain
        (
         schwarz_data,
         schwarz_ops,
         &field_over_subdomains[sub_data.nodal_stride],
         &restricted_field_over_subdomains[sub_data.restricted_nodal_stride],
         i
        );
    }
  }
}

void d4est_solver_schwarz_apply_restrict_transpose_to_restricted_field_over_subdomain
(
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* restricted_field_over_subdomain,
 double* transpose_restricted_field_over_subdomain,
 int subdomain
)
{
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[subdomain];
    for (int j = 0; j < sub_data.num_elements; j++){
      d4est_solver_schwarz_element_data_t ed = sub_data.element_data[j]; 
      d4est_solver_schwarz_operators_apply_schwarz_restrictor
        (
         schwarz_ops,
         &restricted_field_over_subdomain[ed.restricted_nodal_stride],
         (P4EST_DIM),
         &(ed.faces[0]),
         ed.deg,
         schwarz_data->num_nodes_overlap,
         D4OPS_TRANSPOSE,
         &transpose_restricted_field_over_subdomain[ed.nodal_stride]
        );
    }
}
/** 
 * 
 * 
 * @param schwarz_data 
 * @param field_over 
 */
void
d4est_solver_schwarz_apply_restrict_transpose_to_restricted_field_over_subdomains
(
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* restricted_field_over_subdomains,
 double* transpose_restricted_field_over_subdomains
)
{
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[i];
    d4est_solver_schwarz_apply_restrict_transpose_to_restricted_field_over_subdomain
      (
       schwarz_data,
       schwarz_ops,
       &restricted_field_over_subdomains[sub_data.restricted_nodal_stride],
       &transpose_restricted_field_over_subdomains[sub_data.nodal_stride],
       i
      );
  }
}


/** 
 * Take a single subdomain field on a certain process
 * and turn it into a global nodal field (e.g. the parts
 * on ghost elements will be returned to their given 
 * processes). This is useful for visualization.
 * This uses communication.
 * 
 * @param p4est 
 * @param d4est_factors 
 * @param schwarz_data 
 * @param [in] field_on_subdomain restricted field only on a single subdomain
 * @param [out] field 
 * @param mpirank mpirank of core id of subdomain, may be different from p4est mpirank if using multiple cores 
 * @param sub_id id of subdomain
 * @param local_nodes number of nodes for the output field
 * @param is_it_zeroed is the output field zeroed, this is needed because elements not in this subdomain should have zero fields 
 */
void
d4est_solver_schwarz_convert_restricted_subdomain_field_to_global_nodal_field
(
 p4est_t* p4est,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* field_on_subdomain, 
 double* field,
 int mpirank,
 int sub_id,
 int local_nodes,
 d4est_field_init_type_t is_it_zeroed 
)
{
  /* transfer subdomains if multi-process */
  /* copy over relevant data if process matches and element matches, otherwise zero */
  if (is_it_zeroed == FIELD_NOT_ZEROED){
    for (int i = 0; i < local_nodes; i++){
      field[i] = 0.;
    }
  }

  int sub_nodes = schwarz_data->subdomain_data[sub_id].nodal_size;
  double* transpose_restricted_field_on_subdomain = P4EST_ALLOC(double, sub_nodes);
  
  d4est_solver_schwarz_apply_restrict_transpose_to_restricted_field_over_subdomain
    (
     schwarz_data,
     schwarz_ops,
     field_on_subdomain,
     transpose_restricted_field_on_subdomain,
     sub_id
    );
  
  /* Don't need to check if p4est->mpirank == mpirank in this case */
  if (p4est->mpisize == 1){    
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[sub_id];
    for (int j = 0; j < sub_data.num_elements; j++){
      d4est_solver_schwarz_element_data_t ed = sub_data.element_data[j];
      int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed.deg);
      int local_nodal_stride = d4est_factors->element_data[ed.id]->nodal_stride;
      int sub_nodal_stride = ed.nodal_stride;

      d4est_util_copy_1st_to_2nd(&transpose_restricted_field_on_subdomain[sub_nodal_stride],
                                 &field[local_nodal_stride],
                                 volume_nodes);
      
    }
  }
  else {
    /* TRANSER STUFF HERE BEFORE IF */
    if (p4est->mpirank == mpirank){
      /* GO THROUGH LOCAL SUBDOMAINS AND COPY DATA IF NEEDED */
    }
    else {
      /* GO THROUGH GHOST SUBDOMAIN DATA AND COPY DATA */
    }
  }
  
}

void d4est_solver_schwarz_apply_weights_over_subdomain
(
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_solver_schwarz_subdomain_data_t sub_data,
 double* restricted_field_over_subdomain,
 double* weighted_restricted_field_over_subdomain
)
{
  for (int j = 0; j < sub_data.num_elements; j++){
    d4est_solver_schwarz_element_data_t ed = sub_data.element_data[j];

    double* in = &restricted_field_over_subdomain[ed.restricted_nodal_stride];
    double* out = &weighted_restricted_field_over_subdomain[ed.restricted_nodal_stride];
                    
    d4est_solver_schwarz_operators_apply_schwarz_weights
      (
       schwarz_ops,
       in,
       (P4EST_DIM),
       &ed.core_faces[0],
       ed.deg,
       schwarz_data->num_nodes_overlap,
       out
      );
  }
}

void d4est_solver_schwarz_apply_weights_over_all_subdomains
(
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* restricted_field_over_subdomains,
 double* weighted_restricted_field_over_subdomains
)
{

  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[i];
    double* in = &restricted_field_over_subdomains[sub_data.restricted_nodal_stride];
    double* out = &weighted_restricted_field_over_subdomains
                  [sub_data.restricted_nodal_stride];
    d4est_solver_schwarz_apply_weights_over_subdomain
      (
       schwarz_data,
       schwarz_ops,
       sub_data,
       in,
       out
      );
  }
}


void
d4est_solver_schwarz_zero_field_over_subdomain_single_core
(
 p4est_t* p4est,
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_subdomain_data_t sub_data,
 double* field
)
{
  D4EST_ASSERT(p4est->mpisize == 1);
  for (int j = 0; j < sub_data.num_elements; j++){
    d4est_solver_schwarz_element_data_t ed = sub_data.element_data[j];
    int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed.deg);
    int local_nodal_stride = d4est_factors->element_data[ed.id]->nodal_stride;
    for (int i = 0; i < volume_nodes; i++){
      field[local_nodal_stride + i] = 0.;
    }
  }

}

void
d4est_solver_schwarz_compute_correction
(
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* du_restricted_field_over_subdomains,
 double* correction_field_over_subdomains
)
{
  double* weighted_du_restricted_field_over_subdomains
    = P4EST_ALLOC(double, schwarz_data->restricted_nodal_size);
  
  d4est_solver_schwarz_apply_weights_over_all_subdomains
    (
     schwarz_data,
     schwarz_ops,
     du_restricted_field_over_subdomains,
     weighted_du_restricted_field_over_subdomains
    );

  d4est_solver_schwarz_apply_restrict_transpose_to_restricted_field_over_subdomains
    (
     schwarz_data,
     schwarz_ops,
     weighted_du_restricted_field_over_subdomains,
     correction_field_over_subdomains
    );
  
  P4EST_FREE(weighted_du_restricted_field_over_subdomains);
}


void d4est_solver_schwarz_apply_lhs_single_core
(
 p4est_t* p4est,
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_elliptic_eqns_t* elliptic_eqns,
 d4est_elliptic_data_t* elliptic_data,
 double* u_restricted_field_over_subdomains,
 double* Au_restricted_field_over_subdomains,
 double* zeroed_u_field_over_mesh,
 double* zeroed_Au_field_over_mesh,
 int subdomain
)
{
  D4EST_ASSERT(p4est->mpisize == 1);
  int local_nodes = elliptic_data->local_nodes;
  d4est_elliptic_data_t elliptic_data_for_schwarz;

  d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[subdomain];
      
  /* for (int i = 0; i < schwarz_data->num_subdomains; i++){ */
  d4est_solver_schwarz_convert_restricted_subdomain_field_to_global_nodal_field
    (
     p4est,
     d4est_factors,
     schwarz_data,
     schwarz_ops,
     &u_restricted_field_over_subdomains[sub_data.restricted_nodal_stride],
     zeroed_u_field_over_mesh,
     p4est->mpirank,
     subdomain,
     local_nodes,
     FIELD_ZEROED
    );

    d4est_elliptic_data_copy_ptrs
      (
       elliptic_data,
       &elliptic_data_for_schwarz
      );
      
    elliptic_data_for_schwarz.Au = zeroed_Au_field_over_mesh;
    elliptic_data_for_schwarz.u = zeroed_u_field_over_mesh;

    d4est_elliptic_eqns_apply_lhs
      (
       p4est,
       ghost,
       ghost_data,
       elliptic_eqns,
       &elliptic_data_for_schwarz,
       schwarz_ops->d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors
      );

    for (int i = 0; i < local_nodes; i++){
      zeroed_u_field_over_mesh[i] = 0.;
      zeroed_Au_field_over_mesh[i] = 0.;
    }
    
    d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomain
      (
       p4est,
       schwarz_ops->d4est_ops,
       d4est_factors,
       ghost_data,
       schwarz_data,
       schwarz_ops,
       elliptic_data_for_schwarz.Au,
       0,
       &Au_restricted_field_over_subdomains[sub_data.restricted_nodal_stride],
       subdomain
      );
    
  P4EST_FREE(zeroed_u_field_over_mesh);
  P4EST_FREE(zeroed_Au_field_over_mesh);
}

double
d4est_solver_schwarz_visualize_subdomain_helper_single_core
(
 p4est_t* p4est,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 int subdomain,
 int* which_elements,
 double* subdomain_field_for_viz,
 double* weighted_subdomain_field_for_viz
)
{
  d4est_solver_schwarz_subdomain_data_t sub_data
    = schwarz_data->subdomain_data[subdomain];
  double* restricted_subdomain_field = P4EST_ALLOC(double, sub_data.restricted_nodal_size);
  double* weighted_restricted_subdomain_field = P4EST_ALLOC(double, sub_data.restricted_nodal_size);
  for (int i = 0; i < sub_data.restricted_nodal_size;i++){
    restricted_subdomain_field[i] = 1.;
  }
  int local_nodes = d4est_factors->local_sizes.local_nodes;

  d4est_solver_schwarz_apply_weights_over_subdomain
    (
     schwarz_data,
     schwarz_ops,
     sub_data,
     restricted_subdomain_field,
     weighted_restricted_subdomain_field
    );  

  d4est_solver_schwarz_convert_restricted_subdomain_field_to_global_nodal_field
    (
     p4est,
     d4est_factors,
     schwarz_data,
     schwarz_ops,
     restricted_subdomain_field,
     subdomain_field_for_viz,
     p4est->mpirank,
     subdomain,
     local_nodes,
     FIELD_NOT_ZEROED
    );


  /* printf("sub_data.restricted_nodal_size = %d\n", sub_data.restricted_nodal_size); */
  
  /* DEBUG_PRINT_ARR_DBL(restricted_subdomain_field, sub_data.restricted_nodal_size); */

  /* DEBUG_PRINT_ARR_DBL(subdomain_field_for_viz, local_nodes); */
  
d4est_solver_schwarz_convert_restricted_subdomain_field_to_global_nodal_field
    (
     p4est,
     d4est_factors,
     schwarz_data,
     schwarz_ops,
     weighted_restricted_subdomain_field,
     weighted_subdomain_field_for_viz,
     p4est->mpirank,
     subdomain,
     local_nodes,
     FIELD_NOT_ZEROED
    );
    
  for (int i = 0; i < p4est->local_num_quadrants; i++){
    which_elements[i] = 0;
  }
  for (int j = 0; j < sub_data.num_elements; j++){
    d4est_solver_schwarz_element_data_t ed = sub_data.element_data[j];
    which_elements[ed.id] = 1;
  }
  P4EST_FREE(restricted_subdomain_field);
  P4EST_FREE(weighted_restricted_subdomain_field);
}

double
d4est_solver_schwarz_cg_solve_subdomain_single_core
(
 p4est_t* p4est,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_elliptic_data_t* elliptic_data,
 d4est_elliptic_eqns_t* elliptic_eqns,
 double* du_restricted_field_over_subdomain,
 double* rhs_restricted_field_over_subdomain,
 double* zeroed_u_field_over_mesh,
 double* zeroed_Au_field_over_mesh,
 int nodes,
 int iter,
 int subdomain
)
{
  double delta_new, delta_old, temp_max, temp_min, d_dot_Ad;
  double spectral_bound = -1;
  double alpha = -1.;
  double beta = -1.;
  double alpha_old = -1;
  double beta_old = -1;

  double* d = P4EST_ALLOC(double, nodes); 
  double* Ad = P4EST_ALLOC(double, nodes); 
  double* r = P4EST_ALLOC(double, nodes);
  
  d4est_solver_schwarz_apply_lhs_single_core
    (
     p4est,
     schwarz_data,
     schwarz_ops,
     ghost,
     ghost_data,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     elliptic_eqns,
     elliptic_data,
     du_restricted_field_over_subdomain,
     Ad,
     zeroed_u_field_over_mesh,
     zeroed_Au_field_over_mesh,
     subdomain
    );
  
  d4est_util_copy_1st_to_2nd(Ad, r, nodes);
  d4est_linalg_vec_xpby(rhs_restricted_field_over_subdomain, -1., r, nodes);
  d4est_util_copy_1st_to_2nd(r, d, nodes);
  delta_new = d4est_linalg_vec_dot(r,r,nodes);
  
  for (int i = 0; i < iter; i++){

    d4est_solver_schwarz_apply_lhs_single_core
      (
       p4est,
       schwarz_data,
       schwarz_ops,
       ghost,
       ghost_data,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       elliptic_eqns,
       elliptic_data,
       d,
       Ad,
       zeroed_u_field_over_mesh,
       zeroed_Au_field_over_mesh,       
       subdomain
      );
    
    d_dot_Ad = d4est_linalg_vec_dot(d, Ad, nodes);

    alpha_old = alpha;
    alpha = delta_new/d_dot_Ad;
        
    d4est_linalg_vec_axpy(alpha, d, du_restricted_field_over_subdomain, nodes);
    d4est_linalg_vec_axpy(-alpha, Ad, r, nodes);

    delta_old = delta_new;
    delta_new = d4est_linalg_vec_dot(r, r, nodes);
    
    beta_old = beta;
    beta = delta_new/delta_old;
    d4est_linalg_vec_xpby(r, beta, d, nodes);    
    
  }
  
  free(Ad);
  free(d);
  free(r);
  return spectral_bound;
}

void
d4est_solver_schwarz_compute_and_add_correction
(
 p4est_t* p4est,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_elliptic_eqns_t* elliptic_eqns,
 d4est_elliptic_data_t* elliptic_data,
 double* u,
 double* r
)
{

  D4EST_ASSERT(p4est->mpisize == 1);
  int local_nodes = elliptic_data->local_nodes;
  double* du = P4EST_ALLOC_ZERO(double, schwarz_data->restricted_nodal_size); 
  double* restricted_r = P4EST_ALLOC_ZERO(double, schwarz_data->restricted_nodal_size);
  double* correction_field_over_subdomains = P4EST_ALLOC_ZERO(double, schwarz_data->nodal_size);
  
  d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomains
    (
     p4est,
     schwarz_ops->d4est_ops,
     d4est_factors,
     ghost_data,
     schwarz_data,
     schwarz_ops,
     r,
     -1, /* ghost data num of field not known in single core */
     restricted_r
    );

  double* zeroed_u_field_over_mesh = P4EST_ALLOC(double, local_nodes);
  double* zeroed_Au_field_over_mesh = P4EST_ALLOC(double, local_nodes);
  
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[i];
d4est_solver_schwarz_cg_solve_subdomain_single_core
      (
       p4est,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       ghost,
       ghost_data,
       schwarz_data,
       schwarz_ops,
       elliptic_data,
       elliptic_eqns,
       du,
       restricted_r,
       zeroed_u_field_over_mesh,
       zeroed_Au_field_over_mesh,
       sub_data.restricted_nodal_size,
       100,
       i
      );


  }
  
  d4est_solver_schwarz_compute_correction
    (
     schwarz_data,
     schwarz_ops,
     du,
     correction_field_over_subdomains
    );

  /* add correction du to u */
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[i];
    for (int j = 0; j < sub_data.num_elements; j++){
      d4est_solver_schwarz_element_data_t ed = sub_data.element_data[j];

      int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed.deg);
      int local_nodal_stride = d4est_factors->element_data[ed.id]->nodal_stride;
      int sub_nodal_stride = sub_data.nodal_stride + ed.nodal_stride;
      for (int k = 0; k < volume_nodes; k++){
        u[local_nodal_stride + k] += correction_field_over_subdomains[sub_nodal_stride + k];
      }
    }
  }

  P4EST_FREE(zeroed_u_field_over_mesh);
  P4EST_FREE(zeroed_Au_field_over_mesh);
  P4EST_FREE(restricted_r);
  P4EST_FREE(correction_field_over_subdomains);
  P4EST_FREE(du);
}


