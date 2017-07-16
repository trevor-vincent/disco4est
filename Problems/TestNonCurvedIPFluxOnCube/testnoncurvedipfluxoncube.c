#include <sc_reduce.h>
#include <pXest.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <element_data.h>
#include <sipg_flux_scalar_fcns.h>
#include <sipg_flux_vector_fcns.h>
#include <problem.h>
#include <problem_data.h>
#include <d4est_elliptic_eqns.h>
#include <central_flux_params.h>
#include <poisson_operator.h>
#include <krylov_petsc.h>
#include <matrix_sym_tester.h>
#include <dg_norm.h>
#include <d4est_geometry.h>
#include <d4est_geometry_sphere.h>
#include <curved_poisson_debug_vecs.h>
#include <bi_estimator.h>
#include <bi_estimator_flux_fcns.h>
#include <d4est_vtk.h>
#include <ini.h>
#include <cg.h>
#include "time.h"
#include "d4est_util.h"

double pi = 3.1415926535897932384626433832795;

static
double analytic_fcn
(
 double x,
 double y,
 double z
)
{
  return sin(pi*x)*sin(pi*y)*sin(pi*z);
}

static
double boundary_fcn
(
 double x,
 double y,
 double z
)
{
  return analytic_fcn(x,y,z);
}

static
double f_fcn
(
 double x,
 double y,
 double z
)
{
  return -3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z);
}

static int
random_h_refine(p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * quadrant)
{
  /* d4est_element_data_t* data = (d4est_element_data_t*) quadrant->p.user_data; */
  return rand()%2;
}

static int
refine_function
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t *quadrant
)
{
  return 1;
}
/*  */
void
problem_help()
{
}

p4est_connectivity_t*
problem_build_conn()
{
#if (P4EST_DIM)==3
  return p8est_connectivity_new_unitcube();
#elif (P4EST_DIM)==2
  return p4est_connectivity_new_unitsquare();
#else
  D4EST_ABORT("[D4EST_ERROR]: Dim not supported");
#endif
}

typedef struct {
 
  int num_unifrefs;
  int num_randrefs;
  int deg;
  int deg_gauss;
  
  double ip_flux_penalty;  
  int count;
  
} problem_input_t;


static
int problem_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  problem_input_t* pconfig = (problem_input_t*)user;
  if (d4est_util_match_couple(section,"amr",name,"num_unifrefs")) {
    D4EST_ASSERT(pconfig->num_unifrefs == -1);
    pconfig->num_unifrefs = atoi(value);
    pconfig->count += 1;
  }
  else if (d4est_util_match_couple(section,"amr",name,"num_randrefs")) {
    D4EST_ASSERT(pconfig->num_randrefs == -1);
    pconfig->num_randrefs = atoi(value);
    pconfig->count += 1;
  }
  else if (d4est_util_match_couple(section,"amr",name,"deg_gauss")){
    D4EST_ASSERT(pconfig->deg_gauss == -1);
    pconfig->deg_gauss = atoi(value);
    pconfig->count += 1;    
  }
  else if (d4est_util_match_couple(section,"amr",name,"deg")){
    D4EST_ASSERT(pconfig->deg == -1);
    pconfig->deg = atoi(value);
    pconfig->count += 1;    
  }
  else if (d4est_util_match_couple(section,"flux",name,"ip_flux_penalty")) {
    D4EST_ASSERT(pconfig->ip_flux_penalty == -1);
    pconfig->ip_flux_penalty = atof(value);
    pconfig->count += 1;
  } 
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


static
void
build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops
)
{
  prob_vecs->scalar_flux_fcn_data.bndry_fcn = zero_fcn;
  prob_vecs->vector_flux_fcn_data.bndry_fcn = zero_fcn;
  
  poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops);
  d4est_linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);

  prob_vecs->scalar_flux_fcn_data.bndry_fcn = boundary_fcn;
  prob_vecs->vector_flux_fcn_data.bndry_fcn = boundary_fcn;
}

static
void problem_build_rhs
(
 p4est_t* p4est,
 d4est_elliptic_data_t* prob_vecs,
 d4est_elliptic_eqns_t* prob_fcns,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 void* user
)
{
  d4est_poisson_flux_sipg_params_t* ip_flux_params = user;
  double* f = P4EST_ALLOC(double, prob_vecs->local_nodes);
  
  element_data_init_node_vec
    (
     p4est,
     f,
     f_fcn,
     d4est_ops
    );

  /* DEBUG_PRINT_ARR_DBL(f, prob_vecs->local_nodes); */
  
  prob_vecs->vector_flux_fcn_data = sipg_flux_vector_dirichlet_fetch_fcns
                                           (
                                            boundary_fcn,
                                            ip_flux_params
                                           );
  prob_vecs->scalar_flux_fcn_data = sipg_flux_scalar_dirichlet_fetch_fcns
                                           (boundary_fcn);


    
  element_data_apply_mij_on_vec
    (
     p4est,
     f,
     prob_vecs->rhs,
     d4est_ops
    );
  
  d4est_linalg_vec_scale(-1., prob_vecs->rhs, prob_vecs->local_nodes);

  int local_nodes = prob_vecs->local_nodes;
  double* u_eq_0 = P4EST_ALLOC_ZERO(double, local_nodes);
  double* tmp = prob_vecs->u;
  
  prob_vecs->u = u_eq_0; 
  poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops);
  d4est_linalg_vec_axpy(-1., prob_vecs->Au, prob_vecs->rhs, local_nodes);

  prob_vecs->u = tmp;
  P4EST_FREE(u_eq_0);

  /* DEBUG_PRINT_2ARR_DBL(f,prob_vecs->rhs, prob_vecs->local_nodes); */
  prob_vecs->vector_flux_fcn_data = sipg_flux_vector_dirichlet_fetch_fcns
                                           (
                                            zero_fcn,
                                            ip_flux_params
                                           );
  prob_vecs->scalar_flux_fcn_data = sipg_flux_scalar_dirichlet_fetch_fcns
                                           (zero_fcn);

  P4EST_FREE(f);
}


static
problem_input_t
problem_input
(
 const char* input_file
)
{
  int num_of_options = 5;
  
  problem_input_t input;
  input.num_unifrefs = -1;
  input.num_randrefs = -1;
  input.deg_gauss = -1;
  input.deg = -1;
  input.ip_flux_penalty = -1;
  input.count = 0;
  
  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  if (input.count != num_of_options){
    printf("input.count = %d\n", input.count);
    printf("num_of_options = %d\n", num_of_options); 
    D4EST_ABORT("input.count != num_of_options");
  }
  
  /* D4EST_ASSERT(input.count == num_of_options); */
  return input;
}

p4est_geometry_t*
problem_build_geom
(
 p4est_connectivity_t* conn
)
{
  /* D4EST_ASSERT((P4EST_DIM)==3); */
  /* p4est_geometry_t* geom; */
  /* problem_input_t input = problem_input("options.input"); */

  /* geom = p4est_geometry_new_connectivity(conn); */
  
  return NULL;  
}

p4est_t*
problem_build_p4est
(
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t* conn,
 p4est_locidx_t min_quadrants,
 int min_level, 
 int fill_uniform
)
{
  return p4est_new_ext
    (
     mpicomm,
     conn,
     min_quadrants,
     min_level,
     fill_uniform,
     sizeof(element_data_t),
     NULL,
     NULL
    );
}

p4est_t*
problem_load_p4est_from_checkpoint
(
 const char* filename,
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t** conn
){
  int autopartition = 1;
  int load_data = 1;
  int broadcasthead = 0;
  
  return p4est_load_ext (filename,
                         mpicomm,
                         sizeof(element_data_t),
                         load_data,
                         autopartition,
                         broadcasthead,
                         NULL,
                         conn);
}

void
problem_set_degrees
(
 element_data_t* elem_data,
 void* user_ctx
)
{
  problem_input_t* input = user_ctx;
  elem_data->deg = rand()%5 + 1;
  /* elem_data->deg = 4; */
  elem_data->deg_quad = elem_data->deg;
}


void
problem_init
(
 int argc,
 char* argv [],
 p4est_t* p4est,
 p4est_geometry_t* p4est_geom,
 d4est_operators_t* d4est_ops,
 int proc_size,
 sc_MPI_Comm mpicomm,
 int load_from_checkpoint
)
{

  problem_input_t input = problem_input("options.input");
  
  int level;
  /* int degree = input.degree;          */
  /* global_Rinf = input.Rinf; */
  
  D4EST_ASSERT((P4EST_DIM) == 2 || (P4EST_DIM) == 3);
  int world_rank, world_size;
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &world_rank);
  sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &world_size);

  double* Au = P4EST_ALLOC_ZERO(double, 1);
  double* rhs = P4EST_ALLOC_ZERO(double, 1);
  double* u = P4EST_ALLOC_ZERO(double, 1);
  double* u_analytic = P4EST_ALLOC_ZERO(double, 1);
  int local_nodes = 1;

  
  d4est_poisson_flux_sipg_params_t ip_flux_params;
  ip_flux_params.sipg_penalty_prefactor = input.ip_flux_penalty;
  ip_flux_params.sipg_penalty_fcn = sipg_flux_vector_calc_penalty_maxp2_over_minh;
      
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  /* create space for storing the ghost data */
  element_data_t* ghost_data = P4EST_ALLOC (element_data_t,
                                                   ghost->ghosts.elem_count);

  p4est_partition(p4est, 0, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FACE, NULL);
  /* d4est_mesh_geometry_storage_t* geometric_factors = geometric_factors_init(p4est); */


  /* d4est_xyz_fcn_t boundary_flux_fcn = zero_fcn; */
  
  d4est_elliptic_data_t prob_vecs;
  prob_vecs.rhs = rhs;
  prob_vecs.Au = Au;
  prob_vecs.u = u;
  prob_vecs.local_nodes = local_nodes;

  for (level = 0; level < input.num_unifrefs; ++level){

    if (level != 0){
      p4est_refine_ext(p4est,
                       0,
                       -1,
                       refine_function,
                       NULL,
                       NULL
                      );

      p4est_partition(p4est, 0, NULL);
      p4est_balance_ext
        (
         p4est,
         P4EST_CONNECT_FULL,
         NULL,
         NULL
        );

      p4est_ghost_destroy(ghost);
      P4EST_FREE(ghost_data);

      ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
      ghost_data = P4EST_ALLOC(element_data_t, ghost->ghosts.elem_count);      
    }

  }


    
  srand(102230213);
  int num_randrefs = input.num_randrefs;  
  for (int randrefs = 0; randrefs < num_randrefs; ++randrefs){
      p4est_refine_ext(p4est,
                       0,
                       -1,
                       random_h_refine,
                       NULL,
                       NULL
                      );

      p4est_balance_ext
        (
         p4est,
         P4EST_CONNECT_FULL,
         NULL,
         NULL
        );

      p4est_ghost_destroy(ghost);
      P4EST_FREE(ghost_data);

      ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
      ghost_data = P4EST_ALLOC(element_data_t, ghost->ghosts.elem_count);

    /* d4est_element_data_init_new(p4est, geometric_factors, d4est_ops, &dgeom, degree, input.gauss_quad_deg); */

      
  }
  element_data_init_new(p4est,
                        problem_set_degrees,
                        NULL);
  
  
  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.apply_lhs = poisson_apply_aij;
  prob_fcns.build_residual = build_residual;

     
  /* d4est_mesh_geometry_storage_t* geometric_factors = geometric_factors_init(p4est); */


  /* d4est_geometry_t dgeom; */
  /* dgeom.p4est_geom = p4est_geom; */
  /* dgeom.interp_to_gauss = 1; */
  /* dgeom.dxdr_method = INTERP_X_ON_LOBATTO;     */
  /* element_data_init(p4est, geometric_factors, d4est_ops, &dgeom, degree, input.gauss_quad_deg); */
  /* element_data_init_ext(p4est, */
                        /* input.deg, */
                        /* input.deg_gauss); */



    
  /* srand(102230213); */
  /* int num_randrefs = input.num_randrefs; */

  
  /* for (int randrefs = 0; randrefs < num_randrefs; ++randrefs){ */
  /*     p4est_refine_ext(p4est, */
  /*                      0, */
  /*                      -1, */
  /*                      random_h_refine, */
  /*                      NULL, */
  /*                      NULL */
  /*                     ); */

  /*     p4est_balance_ext */
  /*       ( */
  /*        p4est, */
  /*        P4EST_CONNECT_FULL, */
  /*        NULL, */
  /*        NULL */
  /*       ); */

  /*     p4est_ghost_destroy(ghost); */
  /*     P4EST_FREE(ghost_data); */

  /*     ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE); */
  /*     ghost_data = P4EST_ALLOC(element_data_t, ghost->ghosts.elem_count); */

  /*   element_data_init(p4est, geometric_factors, d4est_ops, &dgeom, degree, input.gauss_quad_deg); */

      
  /* } */

  /* element_data_init(p4est, */
  /*                          geometric_factors, */
  /*                          d4est_ops, */
  /*                          &dgeom, degree, */
  /*                          degree_gauss_diff[0], */
  /*                          GAUSS_INTEG); */



    
  local_nodes = element_data_get_local_nodes(p4est);

  Au = P4EST_REALLOC(Au, double, local_nodes);
  u = P4EST_REALLOC(u, double, local_nodes);
  rhs = P4EST_REALLOC(rhs, double, local_nodes);
  u_analytic = P4EST_REALLOC(u_analytic, double, local_nodes);

  prob_vecs.Au = Au;
  prob_vecs.u = u;
  prob_vecs.rhs = rhs;
  prob_vecs.local_nodes = local_nodes;


  prob_vecs.vector_flux_fcn_data = sipg_flux_vector_dirichlet_fetch_fcns
                                          (
                                           zero_fcn,
                                           &ip_flux_params
                                          );

  /* prob_vecs.curved_vector_flux_fcn_data = curved_gauss_central_flux_vector_dirichlet_fetch_fcns */
  /*                                         ( */
  /*                                          zero_fcn, */
  /*                                          &central_flux_params */
  /*                                         ); */


    
  prob_vecs.scalar_flux_fcn_data = sipg_flux_scalar_dirichlet_fetch_fcns
                                          (zero_fcn);

    
  /* d4est_linalg_fill_vec(u, 0., local_nodes); */
  /* element_data_init_node_vec(p4est,f,f_fcn,d4est_ops); */

  /* double total_volume = 0.; */
  /* for (p4est_topidx_t tt = p4est->first_local_tree; */
  /*      tt <= p4est->last_local_tree; */
  /*      ++tt) */
  /*   { */
  /*     p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
  /*     sc_array_t* tquadrants = &tree->quadrants; */
  /*     int Q = (p4est_locidx_t) tquadrants->elem_count; */
  /*     for (int q = 0; q < Q; ++q) { */
  /*       p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q); */
  /*       element_data_t* ed = quad->p.user_data; */
  /*       printf("Tree, Element, Volume = %d, %d, %.25f\n", ed->tree, ed->id, ed->volume); */
  /*       total_volume += ed->volume; */
  /*     } */
  /*   } */

  /* printf("Total volume = %.25f\n", total_volume); */
  /* printf("Theoretical volume = %.25f\n", (4./3.)*M_PI*pow((input.Rinf),3)); */
    
  /* serial_matrix_sym_tester */
  /*   ( */
  /*    p4est, */
  /*    &prob_vecs, /\* only needed for # of nodes *\/ */
  /*    (void*)&prob_fcns, */
  /*    .0000000001, */
  /*    d4est_ops, */
  /*    1, /\* is it curved *\/ */
  /*    3, /\* should we print *\/ */
  /*    &dgeom */
  /*   ); */

  /* p4est_vtk_write_all */
  /*   (p4est, p4est_geom,     /\* we do not need to transform from the vertex space into physical space, so we do not need a p4est_geometry_t * pointer *\/ */
  /*    0.99,    /\* draw each quadrant at almost full scale *\/ */
  /*    1,       /\* do not write the tree id's of each quadrant (there is only one tree in this example) *\/ */
  /*    1,       /\* do write the refinement level of each quadrant *\/ */
  /*    1,       /\* do write the mpi process id of each quadrant *\/ */
  /*    0,       /\* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) *\/ */
  /*    0,       /\* write one scalar field: the solution value *\/ */
  /*    0,       /\* write no vector fields *\/ */
  /*    "non-constant-jacobian" */
  /*   ); */


  /* element_data_init_node_vec */
  /*   ( */
  /*    p4est, */
  /*    u, */
  /*    analytic_fcn, */
  /*    d4est_ops */
  /*   ); */


    
  /* p4est_vtk_write_file */
  /*   (p4est, */
  /*    p4est_geom, */
  /*    "compact-sphere" */
  /*   ); */

  problem_build_rhs
    (
     p4est,
     &prob_vecs,
     &prob_fcns,
     ghost,
     ghost_data,
     d4est_ops,
     &ip_flux_params
    );


  
  
  for (int i = 0; i < local_nodes; i++){
    u[i] = d4est_util_uniform_rand(14234232, 0., 1.);
  }
  
  /* d4est_linalg_fill_vec(u, 1., local_nodes); */

  prob_fcns.apply_lhs(p4est, ghost, ghost_data, &prob_vecs, d4est_ops);

  /* DEBUG_PRINT_2ARR_DBL(prob_vecs.Au, prob_vecs.rhs, prob_vecs.local_nodes); */
  
  printf(
         "d4est_linalg_vec_sum Au, rhs = %.15f, %.15f\n",
         d4est_linalg_vec_sum(prob_vecs.Au, prob_vecs.local_nodes),
         d4est_linalg_vec_sum(prob_vecs.rhs, prob_vecs.local_nodes)
        );



 penalty_calc_t bi_u_penalty_fcn = bi_u_prefactor_conforming_maxp_minh;
  penalty_calc_t bi_u_dirichlet_penalty_fcn = bi_u_prefactor_conforming_maxp_minh;
  penalty_calc_t bi_gradu_penalty_fcn = bi_gradu_prefactor_maxp_minh;
  
  bi_estimator_compute
    (
     p4est,
     &prob_vecs,
     &prob_fcns,
     bi_u_penalty_fcn,
     bi_u_dirichlet_penalty_fcn,
     bi_gradu_penalty_fcn,
     boundary_fcn,
     ip_flux_params.sipg_penalty_prefactor,
     ghost,
     ghost_data,
     d4est_ops
    );


  estimator_stats_t stats;
  estimator_stats_compute(p4est, &stats,0);

  if(world_rank == 0)
    estimator_stats_print(&stats, 0); 
  
  
  /* DEBUG_PRINT_ARR_DBL(prob_vecs.rhs, prob_vecs.local_nodes); */
  /* DEBUG_PRINT_ARR_DBL_SUM(prob_vecs.rhs, prob_vecs.local_nodes); */
  /* krylov_petsc_params_t params = {input.krylov_type, 1, 0, NULL, NULL, NULL}; */
     
  /* krylov_petsc_info_t info = */
  /*   krylov_petsc_solve */
  /*   ( */
  /*    p4est, */
  /*    &prob_vecs, */
  /*    (void*)&prob_fcns, */
  /*    &ghost, */
  /*    (void**)&ghost_data, */
  /*    d4est_ops, */
  /*    NULL, */
  /*    NULL, */
  /*    NULL, */
  /*    NULL */
  /*    /\* &params *\/ */
  /*   ); */


    /* cg_solver_params_t cg_params; */
    /* cg_params.max_iter = 100;// d4est_util_int_pow_int(10, level+1); */
    /* cg_params.rtol = 1e-18; */
    /* cg_params.atol = 1e-18; */
    /* cg_params.mpi_rank = world_rank; */
    /* cg_params.monitor = 1; */
  

    /* cg_solve */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs, */
    /*    &prob_fcns, */
    /*    d4est_ops, */
    /*    ghost, */
    /*    ghost_data, */
    /*    &cg_params */
    /*   ); */

  
  /* double* vertex_u = P4EST_ALLOC(double, (P4EST_CHILDREN)*p4est->local_num_quadrants); */
  /* double* vertex_u_analytic = P4EST_ALLOC(double, (P4EST_CHILDREN)*p4est->local_num_quadrants); */
  /* element_data_init_node_vec(p4est, u_analytic, analytic_fcn, d4est_ops); */
  
  /* element_data_store_nodal_vec_in_vertex_array */
  /*   ( */
  /*    p4est, */
  /*    u, */
  /*    vertex_u */
  /*   ); */

  /* element_data_store_nodal_vec_in_vertex_array */
  /*   ( */
  /*    p4est, */
  /*    u_analytic, */
  /*    vertex_u_analytic */
  /*   ); */


  /* sc_array_t* vertex_u_sc = sc_array_new_data((void*)vertex_u, sizeof(double), (P4EST_CHILDREN)*p4est->local_num_quadrants); */
  /* sc_array_t* vertex_u_analytic_sc = sc_array_new_data((void*)vertex_u_analytic, sizeof(double), (P4EST_CHILDREN)*p4est->local_num_quadrants); */
      
  /* d4est_vtk_context_t* vtk_ctx = d4est_vtk_context_new(p4est, "compact-sphere"); */
  /* d4est_vtk_context_set_geom(vtk_ctx, p4est_geom); */
  /* d4est_vtk_context_set_scale(vtk_ctx, .95); */

  /* vtk_ctx = d4est_vtk_write_header(vtk_ctx); */

  /* vtk_ctx = d4est_vtk_write_point_dataf(vtk_ctx, 2, 0, "vertex_u",vertex_u_sc, "vertex_u_analytic",vertex_u_analytic_sc,vtk_ctx); */

  /* d4est_vtk_write_footer(vtk_ctx); */
  /* /\* d4est_vtk_context_destroy(vtk_ctx); *\/ */


  /* sc_array_destroy(vertex_u_sc); */
  /* sc_array_destroy(vertex_u_analytic_sc);      */

  /* P4EST_FREE(vertex_u); */
  /* P4EST_FREE(vertex_u_analytic); */


  int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
  int vtk_nodes = 0;
     
  int stride = 0;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        /* k++; */
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        element_data_t* ed = quad->p.user_data;
        deg_array[stride] = ed->deg;
        vtk_nodes = d4est_util_int_pow_int(deg_array[stride], (P4EST_DIM))*(P4EST_CHILDREN);
        stride++;
      }
    }


  /* for (int i = 0; i < p4est->local_num_quadrants; i++){ */
  /*   printf("d4est_util_int_pow_int(deg_array[i], (P4EST_DIM)) = %d\n",d4est_util_int_pow_int(deg_array[i], (P4EST_DIM))); */
  /*   printf("deg_array[i] = %d\n",deg_array[i]); */
  /*   printf("(P4EST_DIM) = %d\n",(P4EST_DIM)); */
  /* } */
     


  /* sc_array_t* u_sc = sc_array_new_data((void*)u, sizeof(double), (P4EST_CHILDREN)*p4est->local_num_quadrants); */
  /* sc_array_t* u_analytic_sc = sc_array_new_data((void*)u_analytic, sizeof(double), (P4EST_CHILDREN)*p4est->local_num_quadrants); */

  double* error = P4EST_ALLOC(double, local_nodes);

  element_data_init_node_vec(p4est, u_analytic, analytic_fcn, d4est_ops);

  for (int i = 0; i < local_nodes; i++){
    error[i] = fabs(error[i]);
  }
  
          
  d4est_linalg_vec_axpyeqz(-1., u, u_analytic, error, local_nodes);

  p4est_geometry_t* geom_vtk = p4est_geom;
     

  d4est_vtk_context_t* vtk_ctx = d4est_vtk_dg_context_new(p4est, d4est_ops, "compact-sphere");
  d4est_vtk_context_set_geom(vtk_ctx, geom_vtk);
  d4est_vtk_context_set_scale(vtk_ctx, .99);
  d4est_vtk_context_set_deg_array(vtk_ctx, deg_array);
  vtk_ctx = d4est_vtk_write_dg_header(vtk_ctx, d4est_ops);
  vtk_ctx = d4est_vtk_write_dg_point_dataf(vtk_ctx, 3, 0, "u",u, "u_analytic",u_analytic,"error", error, vtk_ctx);
  vtk_ctx = d4est_vtk_write_dg_cell_dataf
            (
             vtk_ctx,
             1,
             1,
             1,
             0,
             1,
             0,
             0,
             vtk_ctx
            );


  
  d4est_vtk_write_footer(vtk_ctx);
  P4EST_FREE(deg_array);
  P4EST_FREE(error);
  /* p4est_geometry_destroy(geom_vtk); */
  
  d4est_linalg_vec_axpy(-1., u, u_analytic, local_nodes);
  /* d4est_linalg_copy_1st_to_2nd(u, error, local_nodes); */
  /* d4est_linalg_vec_axpy(-1., u, u_analytic, local_nodes); */
  
  
  double local_l2_norm_sqr =  element_data_compute_l2_norm_sqr
                              (
                               p4est,
                               u_analytic,
                               /* local_nodes, */
                               d4est_ops
                              );

  double local_nodes_dbl = (double)local_nodes;
  double local_reduce [2];
  local_reduce[0] = local_nodes_dbl;
  local_reduce[1] = local_l2_norm_sqr;

  double global_reduce [2];

  sc_reduce
    (
     &local_reduce[0],
     &global_reduce[0],
     2,
     sc_MPI_DOUBLE,
     sc_MPI_SUM,
     0,
     sc_MPI_COMM_WORLD
    );

  double global_nodes_dbl = global_reduce[0];
  double global_l2_norm_sqr = global_reduce[1];
    
  if (world_rank == 0){
    /* end = clock(); */
    /* double time_spent = (double)(end - begin) / CLOCKS_PER_SEC; */
    printf
      (
       "\n\n[AMR]: %d %d %.25f %d %.25f %f \n\n",
       /* degree, */
       (int)p4est->global_num_quadrants,
       (int)global_nodes_dbl,
       sqrt(global_l2_norm_sqr),
       -1,//info.iterations,
       -1.,//info.residual_norm,
       0.
      );
  }

  /* geometric_factors_destroy(geometric_factors); */

  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }

  P4EST_FREE(Au);
  P4EST_FREE(rhs);
  P4EST_FREE(u_analytic);
  P4EST_FREE(u);
}
