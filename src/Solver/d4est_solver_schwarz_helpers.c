#include <pXest.h>
#include <d4est_util.h>
#include <d4est_ghost.h>
#include <d4est_operators.h>
#include <d4est_mesh.h>
#include <d4est_solver_schwarz.h>
#include <d4est_solver_schwarz_operators.h>
#include <d4est_solver_schwarz_helpers.h>

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
d4est_solver_schwarz_restrict_nodal_field_to_subdomain
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_data_t* ghost_data,
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* field,
 int ghost_data_num_of_field,
 double* restricted_field_over_subdomains,
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
    d4est_operators_apply_schwarz_restrictor
      (
       schwarz_ops,
       field_ed,
       (P4EST_DIM),
       &(ed.faces[0]),
       ed.deg,
       schwarz_data->num_nodes_overlap,
       D4OPS_NO_TRANSPOSE,
       &restricted_field_over_subdomains[sub_data.restricted_nodal_stride
                                         + ed.restricted_nodal_stride]
      );
  }

}

void
d4est_solver_schwarz_restrict_nodal_field_to_subdomains
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

    d4est_solver_schwarz_restrict_nodal_field_to_subdomain
      (
       p4est,
       d4est_ops,
       d4est_factors,
       ghost_data,
       schwarz_data,
       schwarz_ops,
       field,
       ghost_data_num_of_field,
       restricted_field_over_subdomains,
       i
      );
  }
}

inline
void d4est_solver_schwarz_restrict_transpose_restricted_field_over_single_subdomain
(
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* restricted_field,
 double* transpose_restricted_field,
 int subdomain
)
{
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[subdomain];
    for (int j = 0; j < sub_data.num_elements; j++){
      d4est_solver_schwarz_element_data_t ed = sub_data.element_data[j]; 
      d4est_operators_apply_schwarz_restrictor
        (
         schwarz_ops,
         &restricted_field_over_subdomains[ed.restricted_nodal_stride],
         (P4EST_DIM),
         &(ed.faces[0]),
         ed.deg,
         schwarz_data->num_nodes_overlap,
         D4OPS_TRANSPOSE,
         &field_over_subdomains[ed.nodal_stride]
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
d4est_solver_schwarz_restrict_transpose_restricted_field_over_subdomains
(
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* restricted_field_over_subdomains,
 double* field_over_subdomains
)
{
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[i];
    for (int j = 0; j < sub_data.num_elements; j++){
      d4est_solver_schwarz_element_data_t ed = sub_data.element_data[j]; 
      d4est_operators_apply_schwarz_restrictor
        (
         schwarz_ops,
         &restricted_field_over_subdomains[sub_data.restricted_nodal_stride + ed.restricted_nodal_stride],
         (P4EST_DIM),
         &(ed.faces[0]),
         ed.deg,
         schwarz_data->num_nodes_overlap,
         D4OPS_TRANSPOSE,
         &field_over_subdomains[sub_data.nodal_stride
                                + ed.nodal_stride]
        );
    }
  }
}

void
d4est_solver_schwarz_restrict_field_over_subdomains
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
      
      d4est_operators_apply_schwarz_restrictor
        (
         schwarz_ops,
         &field_over_subdomains[sub_data.nodal_stride
                                + ed.nodal_stride],
         (P4EST_DIM),
         &(ed.faces[0]),
         ed.deg,
         schwarz_data->num_nodes_overlap,
         D4OPS_NO_TRANSPOSE,
         &restricted_field_over_subdomains[sub_data.restricted_nodal_stride
                                           + ed.restricted_nodal_stride]
         
        );
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

  d4est_solver_schwarz_restrict_transpose_restricted_field_over_subdomains
    (
     schwarz_data,
     schwarz_ops,
     weighted_du_restricted_field_over_subdomains,
     correction_field_over_subdomains
    );
  
  P4EST_FREE(weighted_du_restricted_field_over_subdomains);
}

double
d4est_solver_schwarz_cg_solve_single_core
(
 double* A,
 double* x0,
 double* b,
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

  /* double* Ax = D4EST_TEST_ALLOC(double, nodes);  */
  double* d = D4EST_TEST_ALLOC(double, nodes); 
  double* r = D4EST_TEST_ALLOC(double, nodes);

  /* apply_mat(A, x0, Ax, nodes); */

  int sub_res_nodal_stride = schwarz_data->sub_data[subdomain].restricted_nodal_stride;
  
  d4est_solver_schwarz_apply_lhs_single_core
    (
     schwarz_data,
     schwarz_ops,
     elliptic_fcns,
     elliptic_data_for_cg,
     u_restricted_field_over_subdomains,
     Au_restricted_field_over_subdomains,
     subdomain
    )
  
          
  d4est_util_copy_1st_to_2nd(Ax, r, nodes);
  d4est_linalg_vec_xpby(b, -1., r, nodes);
  d4est_util_copy_1st_to_2nd(r, d, nodes);
  delta_new = d4est_linalg_vec_dot(r,r,nodes);
  
  for (int i = 0; i < iter; i++){

    /* apply_mat(A, d, Ax, nodes); */    
    d4est_solver_schwarz_apply_lhs_single_core
      (
       schwarz_data,
       schwarz_ops,
       elliptic_fcns,
       elliptic_data,
       u_restricted_field_over_subdomains,
       Au_restricted_field_over_subdomains,
       subdomain
      )
    
    d_dot_Ad = d4est_linalg_vec_dot(d, Ax, nodes);

    alpha_old = alpha;
    alpha = delta_new/d_dot_Ad;
        
    d4est_linalg_vec_axpy(alpha, d, x0, nodes);
    d4est_linalg_vec_axpy(-alpha, Ax, r, nodes);

    delta_old = delta_new;
    delta_new = d4est_linalg_vec_dot(r, r, nodes);
    
    beta_old = beta;
    beta = delta_new/delta_old;
    d4est_linalg_vec_xpby(r, beta, d, nodes);    
    
  }
  
  free(Ax);
  free(d);
  free(r);
  return spectral_bound;
}




void
d4est_solver_schwarz_compute_and_add_correction
(
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* u,
 double* r
)
{

  D4EST_ASSERT(p4est->mpisize == 1);
  double* du = P4EST_ALLOC_ZERO(double, schwarz_data->restricted_nodal_size); 

  double* restricted_r = P4EST_ALLOC_ZERO(double, schwarz_data->restricted_nodal_size);
  
  d4est_solver_schwarz_restrict_nodal_field_to_subdomains
    (
     p4est,
     d4est_ops,
     d4est_factors,
     ghost_data,
     schwarz_data,
     schwarz_ops,
     r,
     -1, /* ghost data num of field not known in single core */
     restricted_r
    );

  d4est_solver_schwarz_cg_solve_single_core
    (
     

    );
  
  d4est_solver_schwarz_compute_correction
    (
     schwarz_data,
     schwarz_ops,
     du_restricted_field_over_subdomains,
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
  
  P4EST_FREE(weighted_du_restricted_field_over_subdomains);
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

  

d4est_solver_schwarz_restrict_transpose_restricted_field_over_single_subdomain
(
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 field_on_subdomain,
 double* transpose_restricted_field,
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

      d4est_util_copy_1st_to_2nd(&field_over_subdomains[sub_nodal_stride],
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

void
d4est_solver_schwarz_take_step_single_core
(


)
{


}

void d4est_solver_schwarz_apply_lhs_single_core
(
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_elliptic_fcns_t* elliptic_fcns,
 d4est_elliptic_data_t* elliptic_data,
 double* u_restricted_field_over_subdomains,
 double* Au_restricted_field_over_subdomains,
 int subdomain
)
{
  D4EST_ASSERT(p4est->mpisize == 1);
  int local_nodes = elliptic_data->local_nodes;
  double* field = P4EST_ALLOC_ZERO(double, local_nodes);
  double* A_field = P4EST_ALLOC_ZERO(double, local_nodes);
  d4est_elliptic_data_t elliptic_data_for_schwarz;
  
  /* for (int i = 0; i < schwarz_data->num_subdomains; i++){ */
    d4est_solver_schwarz_convert_restricted_subdomain_field_to_global_nodal_field
      (
       p4est,
       d4est_factors,
       schwarz_data,
       &u_restricted_field_over_subdomains[],
       field,
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
      
    elliptic_data_for_schwarz.Au = A_field;
    elliptic_data_for_schwarz.u = field;

    d4est_elliptic_eqns_apply_lhs
      (
       p4est,
       ghost,
       ghost_data,
       elliptic_eqns,
       &elliptic_data_for_schwarz,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors
      );
    
    d4est_solver_restrict_nodal_field_to_subdomain
      (
       p4est,
       d4est_ops,
       d4est_factors,
       ghost_data,
       schwarz_data,
       schwarz_ops,
       A_field,
       0,
       Au_restricted_field_over_subdomains,
       subdomain
      );

  P4EST_FREE(field);
  P4EST_FREE(A_field);
}

 
void d4est_solver_schwarz_apply_weights_over_all_subdomains
(
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* du_restricted_field_over_subdomains,
 double* weighted_du_restricted_field_over_subdomains
)
{

  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[i];
    for (int j = 0; j < sub_data.num_elements; j++){
      d4est_solver_schwarz_element_data_t ed = sub_data.element_data[j];

      double* in = &du_restricted_field_over_subdomains[sub_data.restricted_nodal_stride + ed.restricted_nodal_stride];
      double* out = &weighted_du_restricted_field_over_subdomains[sub_data.restricted_nodal_stride + ed.restricted_nodal_stride];
                    
      d4est_operators_apply_schwarz_weights(
                                            schwarz_ops,
                                            in,
                                            (P4EST_DIM),
                                            &ed.core_faces[0],
                                            ed.deg,
                                            schwarz_data->num_nodes_overlap,
                                            out);
    }
  }
}




d4est_solver_schwarz_zero_field_over_subdomain_single_core
(
 p4est_t* p4est,
 d4est_solver_schwarz_data_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* field,
)
{
  D4EST_ASSERT(p4est->mpisize == 1);
  d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[subdomain];
  for (int j = 0; j < sub_data.num_elements; j++){
    d4est_solver_schwarz_element_data_t ed = sub_data.element_data[j];
    int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed.deg);
    int local_nodal_stride = d4est_factors->element_data[ed.id]->nodal_stride;
    for (int i = 0; i < volume_nodes; i++){
      field[local_nodal_stride + i] = 0.;
    }
  }

}
