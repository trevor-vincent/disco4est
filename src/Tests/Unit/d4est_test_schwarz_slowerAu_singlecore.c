#include <pXest.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_operators.h>
#include <d4est_amr.h>
#include <d4est_amr_random.h>
#include <d4est_vtk.h>
#include <d4est_h5.h>
#include <d4est_checkpoint.h>
#include <d4est_element_data.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_laplacian.h>
#include <d4est_laplacian_flux.h>
#include <d4est_laplacian_flux_sipg.h>
#include <d4est_solver_schwarz_helpers.h>
#include <d4est_solver_schwarz_laplacian.h>
#include <petscsnes.h>
#include <zlog.h>

static void
d4est_test_apply_lhs
(
 p4est_t* p4est,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* user
)
{
  d4est_laplacian_flux_data_t* flux_data_for_apply_lhs = user;  
  d4est_laplacian_apply_aij(p4est, ghost, ghost_data, prob_vecs, flux_data_for_apply_lhs, d4est_ops, d4est_geom, d4est_quad, d4est_factors, 0);
}

double
poly_vec_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
){
#if (P4EST_DIM)==3
  return (x-1.)*x*(y-1.)*y*(z-1)*z;
#else
  return (x-1.)*x*(y-1.)*y;
#endif
}


double
boundary_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void *user
)
{
  return poly_vec_fcn(x,
                      y,
#if(P4EST_DIM)==3
                      z,
#endif
                      user);
}

double
neg_laplacian_poly_vec_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
){
#if (P4EST_DIM)==3
  return - (2 * ((-1 + x) * x * (-1 + y) * y + (-1 + x) * x * (-1 + z) * z + (-1 + y) * y * (-1 + z)*z));
#else
  return - (2.*(x - 1.)*x + 2.*(y - 1.)*y);
#endif
}

static void 
d4est_test_iterate_refine_only_one_element
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  d4est_amr_t* d4est_amr = (d4est_amr_t*) info->p4est->user_pointer;
  int* element_to_refine = (int*) (d4est_amr->scheme->amr_scheme_data);
  d4est_element_data_t* elem_data = (d4est_element_data_t*) info->quad->p.user_data;
  if (elem_data->id == *element_to_refine){
    d4est_amr->refinement_log[elem_data->id] = -elem_data->deg;
  }
  else {
    d4est_amr->refinement_log[elem_data->id] = elem_data->deg;
  }
}


int main(int argc, char *argv[])
{
  sc_MPI_Comm mpicomm;
  PetscInitialize(&argc,&argv,(char*)0,NULL);
  mpicomm = PETSC_COMM_WORLD;
  
  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  
#ifndef NDEBUG
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DEBUG MODE ON\n");
  p4est_init(NULL, SC_LP_ERROR);
  /* p4est_init(NULL, SC_LP_ALWAYS); */
#else
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DEBUG MODE OFF\n");
  p4est_init(NULL, SC_LP_ERROR);
#endif
  
#if (P4EST_DIM)==3
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DIM = 3\n");
#else
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DIM = 2\n");
#endif


  char* input_file = P4EST_ALLOC(char, 100);
  sprintf(input_file, "%s", (argc == 2) ? argv[1] : "d4est_test_schwarz_fasterAu_singlecore.input");

  
  if (proc_rank == 0)
    printf("[D4EST_INFO]: options file = %s\n", input_file);
 
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    input_file,
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse(input_file, d4est_geom);

  p4est_t* p4est;
  p4est = p4est_new_ext
          (
           mpicomm,
           d4est_geom->p4est_conn,
           initial_grid_input->min_quadrants,
           initial_grid_input->min_level,
           initial_grid_input->fill_uniform,
           sizeof(d4est_element_data_t),
           NULL,
           NULL
          );
  

  printf("num_quadrants = %d\n", p4est->local_num_quadrants);

  p4est_partition(p4est, 1, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  
  d4est_ghost_t* d4est_ghost = NULL;
   

  if (proc_rank == 0){
    printf("[D4EST_INFO]: mpisize = %d\n", proc_size);
  }
  if (proc_rank == 0 && initial_grid_input->load_from_checkpoint == 0){
    printf("[D4EST_INFO]: min_quadrants = %d\n", initial_grid_input->min_quadrants);
    printf("[D4EST_INFO]: min_level = %d\n", initial_grid_input->min_level);
    printf("[D4EST_INFO]: fill_uniform = %d\n", initial_grid_input->fill_uniform);
  }
  
  sc_MPI_Barrier(mpicomm);
  printf("[D4EST_INFO]: elements on proc %d = %d\n", proc_rank, p4est->local_num_quadrants);
  sc_MPI_Barrier(mpicomm);
  
  /* start just-in-time dg-math */
  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_data_t* d4est_factors = d4est_mesh_data_init(p4est);
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, input_file, "quadrature");
  

  d4est_mesh_local_sizes_t local_sizes= d4est_mesh_update
                                        (
                                         p4est,
                                         &d4est_ghost,
                                         d4est_ops,
                                         d4est_geom,
                                         d4est_quad,
                                         d4est_factors,
                                         initial_grid_input,
                                         INITIALIZE_GHOST,
                                         INITIALIZE_QUADRATURE_DATA,
                                         INITIALIZE_GEOMETRY_DATA,
                                         INITIALIZE_GEOMETRY_ALIASES,
                                         d4est_mesh_set_initial_extents,
                                         (void*)initial_grid_input
                                        );



  /* create amr scheme */
  int add_amr_in_one_element = 0;
  if (add_amr_in_one_element){
    int* refinement_log = P4EST_ALLOC(int, p4est->local_num_quadrants);
    int element_to_refine = 0;
    d4est_amr_t* d4est_amr = P4EST_ALLOC(d4est_amr_t, 1);
    d4est_amr_scheme_t* scheme = P4EST_ALLOC(d4est_amr_scheme_t, 1);
    scheme->post_balance_callback = NULL;
    scheme->pre_refine_callback = NULL;
    scheme->refine_replace_callback_fcn_ptr = NULL;
    scheme->balance_replace_callback_fcn_ptr = NULL;
    scheme->mark_elements = d4est_test_iterate_refine_only_one_element;
    scheme->amr_scheme_data = &element_to_refine;
    scheme->destroy = NULL;
    
    d4est_amr->mpirank = p4est->mpirank;
    d4est_amr->scheme = scheme;
    d4est_amr->balance_log = NULL;
    d4est_amr->refinement_log = NULL;
    d4est_amr->initial_log = NULL;
    d4est_amr->max_degree = 1000;
    
    d4est_amr_step
      (
       p4est,
       NULL,
       d4est_amr,
       NULL,
       NULL,
       NULL
      );

    P4EST_FREE(refinement_log);
    P4EST_FREE(d4est_amr);
    P4EST_FREE(scheme);

    local_sizes = d4est_mesh_update
                  (
                   p4est,
                   &d4est_ghost,
                   d4est_ops,
                   d4est_geom,
                   d4est_quad,
                   d4est_factors,
                   initial_grid_input,
                   INITIALIZE_GHOST,
                   INITIALIZE_QUADRATURE_DATA,
                   INITIALIZE_GEOMETRY_DATA,
                   INITIALIZE_GEOMETRY_ALIASES,
                   d4est_mesh_set_quadratures_after_amr,
                   (void*)initial_grid_input
                  );
  
    initial_grid_input->initial_nodes = local_sizes.local_nodes;
  }


       
    dirichlet_bndry_eval_method_t eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;

    /* / Setup boundary conditions */
    d4est_laplacian_dirichlet_bc_t bc_data_for_lhs;
    bc_data_for_lhs.dirichlet_fcn = poly_vec_fcn;
    bc_data_for_lhs.eval_method = eval_method;  

    
    d4est_laplacian_flux_data_t* flux_data_for_apply_lhs = d4est_laplacian_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_for_lhs);

    d4est_elliptic_eqns_t prob_fcns;
    prob_fcns.build_residual = NULL;
    prob_fcns.apply_lhs = d4est_test_apply_lhs;
    prob_fcns.user = flux_data_for_apply_lhs;
    d4est_field_type_t field_type = NODAL;
    
    double* poly_vec = P4EST_ALLOC(double, local_sizes.local_nodes);
    double* Apoly_vec = P4EST_ALLOC(double, local_sizes.local_nodes);
    double* f = P4EST_ALLOC(double, local_sizes.local_nodes);
    double* Apoly_vec_compare = P4EST_ALLOC(double, local_sizes.local_nodes);
    d4est_elliptic_data_t elliptic_data;
    elliptic_data.u = poly_vec;
    elliptic_data.Au = Apoly_vec;
    elliptic_data.local_nodes = local_sizes.local_nodes;
    elliptic_data.field_types = &field_type;
    elliptic_data.num_of_fields = 1;

    d4est_ghost_data_t* d4est_ghost_data = d4est_ghost_data_init(p4est,
                                                                 d4est_ghost,
                                                                 &field_type,
                                                                 1);

    
    
    d4est_mesh_init_field
      (
       p4est,
       poly_vec,
       poly_vec_fcn,
       d4est_ops, // unnecessary?
       d4est_geom, // unnecessary?
       d4est_factors,
       INIT_FIELD_ON_LOBATTO,
       NULL
      );

    d4est_mesh_init_field
      (
       p4est,
       f,
       neg_laplacian_poly_vec_fcn,
       d4est_ops, // unnecessary?
       d4est_geom, // unnecessary?
       d4est_factors,
       INIT_FIELD_ON_LOBATTO,
       NULL
      );
      
    int* element_id = P4EST_ALLOC(int, p4est->local_num_quadrants);
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
          d4est_quadrature_volume_t mesh_object;
          mesh_object.dq = ed->dq;
          mesh_object.tree = ed->tree;
          mesh_object.element_id = ed->id;
          mesh_object.q[0] = ed->q[0];
          mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
          mesh_object.q[2] = ed->q[2];
#endif
          element_id[ed->id] = ed->id;
          
          double* J_quad = d4est_mesh_get_jacobian_on_quadrature_points(d4est_factors,
                                                                        ed);

          /* printf("ed->deg = %d\n", ed->deg); */
          /* D4EST_ASSERT(ed->deg >= 2); */
        
          d4est_quadrature_apply_mass_matrix
            (
             d4est_ops,
             d4est_geom,
             d4est_quad,
             &mesh_object,
             QUAD_OBJECT_VOLUME,
             QUAD_INTEGRAND_UNKNOWN,
             &f[ed->nodal_stride],
             ed->deg,
             J_quad,
             ed->deg_quad,
             &Apoly_vec_compare[ed->nodal_stride]
            );
        }
      }
    P4EST_FREE(f);
      
    
  int num_nodes_overlap = 2;
  d4est_solver_schwarz_metadata_t* schwarz_data
    = d4est_solver_schwarz_metadata_init
    (
     p4est,
     d4est_ghost,
     input_file,
     "d4est_solver_schwarz"
    );

  d4est_solver_schwarz_operators_t* schwarz_ops
    = d4est_solver_schwarz_operators_init
    (d4est_ops);


  /* d4est_solver_schwarz_metadata_print */
  /*   ( */
  /*    p4est, */
  /*    schwarz_data */
  /*   ); */
  
  /* d4est_mesh_print_out_mortar_data */
    /* ( */
     /* d4est_factors */
    /* ); */

  
  double* poly_vec_over_subdomains = P4EST_ALLOC(double,
                                                 schwarz_data->restricted_nodal_size);

  double* Apoly_vec_over_subdomains = P4EST_ALLOC(double,
                                                    schwarz_data->restricted_nodal_size);

  double* Apoly_vec_over_subdomains_2 = P4EST_ALLOC(double,
                                                    schwarz_data->restricted_nodal_size);
  
  double* Apoly_vec_compare_over_subdomains = P4EST_ALLOC(double,
                                                    schwarz_data->restricted_nodal_size);
    

  printf("restricted_nodal_size = %d\n", schwarz_data->restricted_nodal_size);

  double* u = P4EST_ALLOC_ZERO(double, local_sizes.local_nodes);
  double* Au = P4EST_ALLOC_ZERO(double, local_sizes.local_nodes);

  d4est_elliptic_eqns_apply_lhs
    (
     p4est,
     d4est_ghost,
     d4est_ghost_data,
     &prob_fcns,
     &elliptic_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors
    );

  d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomains
    (
     p4est,
     d4est_ops,
     d4est_factors,
     d4est_ghost,
     d4est_ghost_data,
     schwarz_data,
     schwarz_ops,
     Apoly_vec_compare,
     -1,
     Apoly_vec_compare_over_subdomains
    );
    
  d4est_solver_schwarz_convert_nodal_field_to_restricted_field_over_subdomains
    (
     p4est,
     d4est_ops,
     d4est_factors,
     d4est_ghost,
     d4est_ghost_data,
     schwarz_data,
     schwarz_ops,
     poly_vec,
     -1,
     poly_vec_over_subdomains
    );


  /* DEBUG_PRINT_ARR_DBL(poly_vec_over_subdomains,schwarz_data->restricted_nodal_size); */
  
  for (int i = 0; i < schwarz_data->num_subdomains; i++){

    d4est_solver_schwarz_subdomain_metadata_t* sub_data = &schwarz_data->subdomain_metadata[i];
    
    d4est_solver_schwarz_apply_lhs_single_core
      (
       p4est,
       schwarz_data,
       schwarz_ops,
       d4est_ghost,
       d4est_ghost_data,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       &prob_fcns,
       &elliptic_data,
       &poly_vec_over_subdomains[sub_data->restricted_nodal_stride],
       &Apoly_vec_over_subdomains[sub_data->restricted_nodal_stride],
       u,
       Au,
       i,
       NULL
      );

    /* d4est_solver_schwarz_laplacian_apply_over_subdomain */
    /*   ( */
    /*    p4est, */
    /*    d4est_ops, */
    /*    d4est_geom, */
    /*    d4est_quad, */
    /*    d4est_factors, */
    /*    d4est_ghost, */
    /*    schwarz_data, */
    /*    schwarz_ops, */
    /*    flux_data_for_apply_lhs, */
    /*    &poly_vec_over_subdomains[sub_data->restricted_nodal_stride], */
    /*    &Apoly_vec_over_subdomains_2[sub_data->restricted_nodal_stride], */
    /*    i */
    /*   ); */


    /* double* temp2 = &Apoly_vec_over_subdomains_2[sub_data->restricted_nodal_stride]; */
    /* double* temp1 = &Apoly_vec_over_subdomains[sub_data->restricted_nodal_stride]; */

    /* DEBUG_PRINT_2ARR_DBL(temp1, temp2, sub_data->restricted_nodal_size); */
    /* DEBUG_PRINT_ARR_DBL(temp1, sub_data->restricted_nodal_size); */
    
    
      /* DEBUG_PRINT_2ARR_DBL(Apoly_vec_over_subdomains_2, */
                           /* Apoly_vec_over_subdomains_2, */
                           /* sub_data->restricted_nodal_size); */
  }

  printf("\n******RESULT******\n");
  /* DEBUG_PRINT_2ARR_DBL(poly_vec_over_subdomains, */
  /*                      Apoly_vec_over_subdomains, */
  /*                      schwarz_data->restricted_nodal_size); */

  
  DEBUG_PRINT_ARR_DBL_SUM(Apoly_vec_over_subdomains,
                          schwarz_data->restricted_nodal_size);

  if (d4est_ghost_data != NULL){
    d4est_ghost_data_destroy(d4est_ghost_data);
    d4est_ghost_data = NULL;
  } 
  
  P4EST_FREE(u);
  P4EST_FREE(poly_vec);
  P4EST_FREE(Apoly_vec);
  P4EST_FREE(Apoly_vec_compare);
  P4EST_FREE(Au);
  P4EST_FREE(poly_vec_over_subdomains);
  P4EST_FREE(Apoly_vec_over_subdomains);
  P4EST_FREE(Apoly_vec_over_subdomains_2);
  P4EST_FREE(element_id);
  P4EST_FREE(Apoly_vec_compare_over_subdomains);
  d4est_laplacian_flux_destroy(flux_data_for_apply_lhs);
  
  d4est_solver_schwarz_metadata_destroy
    (
     schwarz_data
    );

  d4est_solver_schwarz_operators_destroy
    (
     schwarz_ops
    );

  d4est_mesh_initial_extents_destroy(initial_grid_input);
  d4est_mesh_data_destroy(d4est_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  
  if (d4est_ghost) {
    d4est_ghost_destroy(d4est_ghost);
  }
    
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);

  P4EST_FREE(input_file);
  
  PetscFinalize();
  return 0;
}
