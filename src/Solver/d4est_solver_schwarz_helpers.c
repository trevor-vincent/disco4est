#include <pXest.h>
#include <d4est_util.h>
#include <d4est_ghost.h>
#include <d4est_operators.h>
#include <d4est_mesh.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_solver_schwarz_operators.h>
#include <d4est_solver_schwarz_helpers.h>
#include <d4est_field.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_linalg.h>
#include <d4est_laplacian_flux.h>
#include <d4est_solver_schwarz_laplacian.h>

/** 
 * Returns the element id within the schwarz subdomain
 * if the element exists in the subdomain otherwise 
 * it returns zero.
 * 
 * @param sub_data 
 * @param element_tree 
 * @param element_tree_quadid 
 * 
 * @return subdomain element id
 */
int 
d4est_solver_schwarz_is_element_in_subdomain
(
 d4est_solver_schwarz_subdomain_metadata_t* sub_data,
 int element_mpirank,
 int element_tree,
 int element_tree_quadid
)
{
  for (int i = 0; i < sub_data->num_elements; i++){
    int mpirank = sub_data->element_metadata[i].mpirank;
    int tree = sub_data->element_metadata[i].tree;
    int tree_quadid = sub_data->element_metadata[i].tree_quadid;
    if (tree == element_tree
        &&
        tree_quadid == element_tree_quadid
        &&
        mpirank == element_mpirank
       ){
      return i;
    }
  }
  return -1;
}

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
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* field,
 int ghost_data_num_of_field,
 double* restricted_field_over_subdomain,
 int subdomain
)
{
  d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[subdomain];
  for (int j = 0; j < sub_data.num_elements; j++){
    d4est_solver_schwarz_element_metadata_t schwarz_ed = sub_data.element_metadata[j];
    double* field_ed = NULL;      
    if (schwarz_ed.mpirank == p4est->mpirank){
      d4est_element_data_t* mesh_ed = d4est_element_data_get_ptr(
                                                                 p4est,
                                                                 schwarz_ed.tree,
                                                                 schwarz_ed.tree_quadid
                                                                );
      field_ed = &field[mesh_ed->nodal_stride];
    }
    else {
      field_ed =
        d4est_ghost_data_get_field_on_element
        (
         &ghost->ghost_elements[schwarz_ed.id],/* d4est_factors->element_data[ed.id + p4est->local_num_quadrants], */
         /* d4est_element_data_get_ptr(p4est,ed.tree,schwarz_ed.tree_quadid), */
         ghost_data_num_of_field,
         ghost_data
        );

      /* int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), */
                                             /* ghost->ghost_elements[schwarz_ed.id].deg); */
      /* DEBUG_PRINT_MPI_ARR_DBL(p4est->mpirank, field_ed, volume_nodes); */
    }
    d4est_solver_schwarz_operators_apply_schwarz_restrictor
      (
       schwarz_ops,
       field_ed,
       (P4EST_DIM),
       &(schwarz_ed.faces[0]),
       schwarz_ed.deg,
       schwarz_data->num_nodes_overlap,
       D4OPS_NO_TRANSPOSE,
       &restricted_field_over_subdomain[schwarz_ed.restricted_nodal_stride]
      );
    /* double* tmp =&restricted_field_over_subdomain[schwarz_ed.restricted_nodal_stride]; */
    /* DEBUG_PRINT_MPI_ARR_DBL(p4est->mpirank, */
                            /* tmp, */
                            /* schwarz_ed.restricted_nodal_size); */
                            
  }
}

void
d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomains
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* field,
 int ghost_data_num_of_field,
 double* restricted_field_over_subdomains
)
{
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[i];
    d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomain
      (
       p4est,
       d4est_ops,
       d4est_factors,
       ghost,
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
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* field_over_subdomains,
 double* restricted_field_over_subdomains,
 int subdomain
)
{
  d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[subdomain];
  for (int i = 0; i < sub_data.num_elements; i++){
    d4est_solver_schwarz_element_metadata_t ed = sub_data.element_metadata[i];
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
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* field_over_subdomains,
 double* restricted_field_over_subdomains
)
{
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[i];
    for (int j = 0; j < sub_data.num_elements; j++){
      d4est_solver_schwarz_element_metadata_t ed = sub_data.element_metadata[j];

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
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* restricted_field_over_subdomain,
 double* transpose_restricted_field_over_subdomain,
 int subdomain
)
{
    d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[subdomain];
    for (int j = 0; j < sub_data.num_elements; j++){
      d4est_solver_schwarz_element_metadata_t ed = sub_data.element_metadata[j]; 
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
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* restricted_field_over_subdomains,
 double* transpose_restricted_field_over_subdomains
)
{
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[i];
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
 d4est_solver_schwarz_metadata_t* schwarz_data,
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

  int sub_nodes = schwarz_data->subdomain_metadata[sub_id].nodal_size;
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
    d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[sub_id];
    for (int j = 0; j < sub_data.num_elements; j++){
      d4est_solver_schwarz_element_metadata_t schwarz_ed = sub_data.element_metadata[j];
      d4est_element_data_t* mesh_ed = d4est_element_data_get_ptr(p4est, schwarz_ed.tree, schwarz_ed.tree_quadid);
      
      int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), schwarz_ed.deg);
      int local_nodal_stride = mesh_ed->nodal_stride;
      int sub_nodal_stride = schwarz_ed.nodal_stride;

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
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_solver_schwarz_subdomain_metadata_t sub_data,
 double* restricted_field_over_subdomain,
 double* weighted_restricted_field_over_subdomain
)
{
  for (int j = 0; j < sub_data.num_elements; j++){
    d4est_solver_schwarz_element_metadata_t ed = sub_data.element_metadata[j];

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
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 double* restricted_field_over_subdomains,
 double* weighted_restricted_field_over_subdomains
)
{
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[i];
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
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_subdomain_metadata_t sub_data,
 double* field
)
{
  D4EST_ASSERT(p4est->mpisize == 1);
  for (int j = 0; j < sub_data.num_elements; j++){
    d4est_solver_schwarz_element_metadata_t schwarz_ed = sub_data.element_metadata[j];
    d4est_element_data_t* mesh_ed = d4est_element_data_get_ptr(p4est, schwarz_ed.tree, schwarz_ed.tree_quadid);
    
    int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), mesh_ed->deg);
    int local_nodal_stride = mesh_ed->nodal_stride;
    for (int i = 0; i < volume_nodes; i++){
      field[local_nodal_stride + i] = 0.;
    }
  }
}

void
d4est_solver_schwarz_compute_correction
(
 d4est_solver_schwarz_metadata_t* schwarz_data,
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
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_elliptic_eqns_t* elliptic_eqns,
 d4est_elliptic_data_t* elliptic_data,
 double* u_restricted_field_over_subdomain,
 double* Au_restricted_field_over_subdomain,
 double* zeroed_u_field_over_mesh,
 double* zeroed_Au_field_over_mesh,
 int subdomain,
 d4est_vtk_helper_array_t* helper_array
)
{
  D4EST_ASSERT(p4est->mpisize == 1);
  int local_nodes = elliptic_data->local_nodes;
  d4est_elliptic_data_t elliptic_data_for_schwarz;

  d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[subdomain];
  
  d4est_solver_schwarz_convert_restricted_subdomain_field_to_global_nodal_field
    (
     p4est,
     d4est_factors,
     schwarz_data,
     schwarz_ops,
     u_restricted_field_over_subdomain,
     zeroed_u_field_over_mesh,
     p4est->mpirank,
     subdomain,
     local_nodes,
     FIELD_NOT_ZEROED
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

    if(helper_array != NULL){
      double* Au_sub = d4est_vtk_helper_array_alloc_and_add_nodal_dbl_field
                       (
                        helper_array,
                        "Au_fcn_sub",
                        subdomain
                       );

      d4est_util_copy_1st_to_2nd(elliptic_data_for_schwarz.Au, Au_sub, elliptic_data_for_schwarz.local_nodes);
    }
    
    d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomain
      (
       p4est,
       schwarz_ops->d4est_ops,
       d4est_factors,
       ghost,
       ghost_data,
       schwarz_data,
       schwarz_ops,
       elliptic_data_for_schwarz.Au,
       0,
       Au_restricted_field_over_subdomain,
       subdomain
      );
}

double
d4est_solver_schwarz_visualize_subdomain_helper_single_core
(
 p4est_t* p4est,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 int subdomain,
 int* which_elements,
 double* subdomain_field_for_viz,
 double* weighted_subdomain_field_for_viz
 /* double* weighted_subdomain_field_for_viz_test */
)
{
  d4est_solver_schwarz_subdomain_metadata_t sub_data
    = schwarz_data->subdomain_metadata[subdomain];
  double* restricted_subdomain_field = P4EST_ALLOC(double, sub_data.restricted_nodal_size);
  double* weighted_restricted_subdomain_field = P4EST_ALLOC(double, sub_data.restricted_nodal_size);
  /* double* weighted_restricted_subdomain_field_test = P4EST_ALLOC(double, sub_data.restricted_nodal_size); */
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
    d4est_solver_schwarz_element_metadata_t ed = sub_data.element_metadata[j];
    which_elements[ed.id] = 1;
  }
  P4EST_FREE(restricted_subdomain_field);
  P4EST_FREE(weighted_restricted_subdomain_field);
  /* P4EST_FREE(weighted_restricted_subdomain_field_test); */
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
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_elliptic_data_t* elliptic_data,
 d4est_elliptic_eqns_t* elliptic_eqns,
 double* du_restricted_field_over_subdomain,
 double* rhs_restricted_field_over_subdomain,
 double* zeroed_u_field_over_mesh,
 double* zeroed_Au_field_over_mesh,
 int iter,
 double atol,
 double rtol,
 int subdomain
)
{


  int nodes = schwarz_data->subdomain_metadata[subdomain].restricted_nodal_size;
  
  double delta_new, delta_old, temp_max, temp_min, d_dot_Ad;
  double spectral_bound = -1;
  double alpha = -1.;
  double beta = -1.;
  double alpha_old = -1;
  double beta_old = -1;

  /* printf("nodes = %d\n", nodes); */
  
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
     subdomain,
     NULL
    );
  
  d4est_util_copy_1st_to_2nd(Ad, r, nodes);
  d4est_linalg_vec_xpby(rhs_restricted_field_over_subdomain, -1., r, nodes);
  d4est_util_copy_1st_to_2nd(r, d, nodes);
  delta_new = d4est_linalg_vec_dot(r,r,nodes);
  double delta_0 = delta_new;
  
  int i;
  for (i = 0; i < iter && (delta_new > atol*atol + delta_0 * rtol*rtol); i++){

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
       subdomain,
       NULL
      );
    
    d_dot_Ad = d4est_linalg_vec_dot(d, Ad, nodes);

    alpha_old = alpha;
    alpha = delta_new/d_dot_Ad;
        
    d4est_linalg_vec_axpy(alpha, d, du_restricted_field_over_subdomain, nodes);
    d4est_linalg_vec_axpy(-alpha, Ad, r, nodes);

    delta_old = delta_new;
    delta_new = d4est_linalg_vec_dot(r, r, nodes);
    /* printf("delta_new = %.25f %.25f %d\n", sqrt(delta_new), sqrt(atol*atol + delta_0 * rtol*rtol), (delta_new > atol*atol + delta_0 * rtol*rtol)); */
    /* printf("rtol, atol = %.25f %.25f\n", rtol, atol); */
    beta_old = beta;
    beta = delta_new/delta_old;
    d4est_linalg_vec_xpby(r, beta, d, nodes);    
    
  }

  printf("SOLVE FOR SUBDOMAIN %d, iter %d, r %.15f\n", subdomain,
         i, sqrt(delta_new));
  
  P4EST_FREE(Ad);
  P4EST_FREE(d);
  P4EST_FREE(r);
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
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_elliptic_eqns_t* elliptic_eqns,
 d4est_elliptic_data_t* elliptic_data,
 double* u,
 double* r,
 double iter,
 double atol,
 double rtol,
 d4est_vtk_helper_array_t* helper_array,
 int suffix_id
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
     ghost,
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
    d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[i];
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
       &du[sub_data.restricted_nodal_stride],
       &restricted_r[sub_data.restricted_nodal_stride],
       zeroed_u_field_over_mesh,
       zeroed_Au_field_over_mesh,
       iter,
       atol,
       rtol,
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

  if (helper_array != NULL){

    for (int i = 0; i < schwarz_data->num_subdomains; i++){
      d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[i];
      char string_vtk [25];
      sprintf(string_vtk, "%s_%d", "du_sub", i);
      double* du_sub = d4est_vtk_helper_array_alloc_and_add_nodal_dbl_field
                   (
                    helper_array,
                    &string_vtk[0],
                    suffix_id
                   );

      d4est_solver_schwarz_convert_restricted_subdomain_field_to_global_nodal_field
        (
         p4est,
         d4est_factors,
         schwarz_data,
         schwarz_ops,
         &du[sub_data.restricted_nodal_stride],
         du_sub,
         p4est->mpirank,
         i,
         local_nodes,
         FIELD_NOT_ZEROED
        );
    }
  }

  /* add correction du to u */
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[i];
    for (int j = 0; j < sub_data.num_elements; j++){
      d4est_solver_schwarz_element_metadata_t schwarz_ed = sub_data.element_metadata[j];
      d4est_element_data_t* mesh_ed = d4est_element_data_get_ptr(p4est, schwarz_ed.tree, schwarz_ed.tree_quadid);

      int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), mesh_ed->deg);
      int local_nodal_stride = mesh_ed->nodal_stride;
      int sub_nodal_stride = sub_data.nodal_stride + schwarz_ed.nodal_stride;
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




double
d4est_solver_schwarz_cg_solve_subdomain_single_core_version2
(
 p4est_t* p4est,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_laplacian_flux_data_t* flux_fcn_data,
 d4est_solver_schwarz_laplacian_mortar_data_t* laplacian_mortar_data,
 double* du_restricted_field_over_subdomain,
 double* rhs_restricted_field_over_subdomain,
 int iter,
 double atol,
 double rtol,
 int subdomain
)
{
  int nodes = schwarz_data->subdomain_metadata[subdomain].restricted_nodal_size;
  
  double delta_new, delta_old, temp_max, temp_min, d_dot_Ad;
  double spectral_bound = -1;
  double alpha = -1.;
  double beta = -1.;
  double alpha_old = -1;
  double beta_old = -1;

  /* printf("nodes = %d\n", nodes); */
  
  double* d = P4EST_ALLOC(double, nodes); 
  double* Ad = P4EST_ALLOC(double, nodes); 
  double* r = P4EST_ALLOC(double, nodes);
  
  D4EST_ABORT("Need to add support for this again");
  /* d4est_solver_schwarz_laplacian_apply_over_subdomain */
  /*   ( */
  /*    p4est, */
  /*    schwarz_ops->d4est_ops, */
  /*    d4est_geom, */
  /*    d4est_quad, */
  /*    d4est_factors, */
  /*    ghost, */
  /*    schwarz_data, */
  /*    schwarz_ops, */
  /*    flux_fcn_data, */
  /*    laplacian_mortar_data, */
  /*    du_restricted_field_over_subdomain, */
  /*    Ad, */
  /*    subdomain */
  /*   ); */
  
  
  d4est_util_copy_1st_to_2nd(Ad, r, nodes);
  d4est_linalg_vec_xpby(rhs_restricted_field_over_subdomain, -1., r, nodes);
  d4est_util_copy_1st_to_2nd(r, d, nodes);
  delta_new = d4est_linalg_vec_dot(r,r,nodes);
  double delta_0 = delta_new;
  
  
  for (int i = 0; i < iter && (delta_new > atol*atol + delta_0 * rtol*rtol); i++){

    D4EST_ABORT("Add support for this again");
  /* d4est_solver_schwarz_laplacian_apply_over_subdomain */
  /*   ( */
  /*    p4est, */
  /*    schwarz_ops->d4est_ops, */
  /*    d4est_geom, */
  /*    d4est_quad, */
  /*    d4est_factors, */
  /*    ghost, */
  /*    schwarz_data, */
  /*    schwarz_ops, */
  /*    flux_fcn_data, */
  /*    laplacian_mortar_data, */
  /*    d, */
  /*    Ad, */
  /*    subdomain */
  /*   ); */
    
    d_dot_Ad = d4est_linalg_vec_dot(d, Ad, nodes);

    alpha_old = alpha;
    alpha = delta_new/d_dot_Ad;
        
    d4est_linalg_vec_axpy(alpha, d, du_restricted_field_over_subdomain, nodes);
    d4est_linalg_vec_axpy(-alpha, Ad, r, nodes);

    delta_old = delta_new;
    delta_new = d4est_linalg_vec_dot(r, r, nodes);
    /* printf("delta_new = %.25f %.25f %d\n", sqrt(delta_new), sqrt(atol*atol + delta_0 * rtol*rtol), (delta_new > atol*atol + delta_0 * rtol*rtol)); */
    /* printf("rtol, atol = %.25f %.25f\n", rtol, atol); */
    beta_old = beta;
    beta = delta_new/delta_old;
    d4est_linalg_vec_xpby(r, beta, d, nodes);    
    
  }
  
  P4EST_FREE(Ad);
  P4EST_FREE(d);
  P4EST_FREE(r);
  return spectral_bound;
}

void
d4est_solver_schwarz_compute_and_add_correction_version2
(
 p4est_t* p4est,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_laplacian_flux_data_t* flux_fcn_data,
 d4est_solver_schwarz_laplacian_mortar_data_t* laplacian_mortar_data,
 double* u,
 double* r,
 int subdomain_iter,
 double subdomain_atol,
 double subdomain_rtol
)
{

  D4EST_ASSERT(p4est->mpisize == 1);
  int local_nodes = d4est_factors->local_sizes.local_nodes;
  double* du = P4EST_ALLOC_ZERO(double, schwarz_data->restricted_nodal_size); 
  double* restricted_r = P4EST_ALLOC_ZERO(double, schwarz_data->restricted_nodal_size);
  double* correction_field_over_subdomains = P4EST_ALLOC_ZERO(double, schwarz_data->nodal_size);
  
  d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomains
    (
     p4est,
     schwarz_ops->d4est_ops,
     d4est_factors,
     ghost,
     ghost_data,
     schwarz_data,
     schwarz_ops,
     r,
     -1, /* ghost data num of field not known in single core */
     restricted_r
    );

  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[i];
d4est_solver_schwarz_cg_solve_subdomain_single_core_version2
      (
       p4est,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       ghost,
       ghost_data,
       schwarz_data,
       schwarz_ops,
       flux_fcn_data,
       laplacian_mortar_data,
       &du[sub_data.restricted_nodal_stride],
       &restricted_r[sub_data.restricted_nodal_stride],
       subdomain_iter,
       subdomain_atol,
       subdomain_rtol,
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
    d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[i];
    for (int j = 0; j < sub_data.num_elements; j++){
      d4est_solver_schwarz_element_metadata_t schwarz_ed = sub_data.element_metadata[j];
      d4est_element_data_t* mesh_ed = d4est_element_data_get_ptr(
                                                                 p4est,
                                                                 schwarz_ed.tree,
                                                                 schwarz_ed.tree_quadid
                                                                );
      int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), mesh_ed->deg);
      int local_nodal_stride = mesh_ed->nodal_stride;
      int sub_nodal_stride = sub_data.nodal_stride + schwarz_ed.nodal_stride;
      for (int k = 0; k < volume_nodes; k++){
        u[local_nodal_stride + k] += correction_field_over_subdomains[sub_nodal_stride + k];
      }
    }
  }

  P4EST_FREE(restricted_r);
  P4EST_FREE(correction_field_over_subdomains);
  P4EST_FREE(du);
}
