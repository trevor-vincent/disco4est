#include <sc_reduce.h>
#include <pXest.h>
#include <util.h>
#include <d4est_linalg.h>
#include <problem.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_estimator_bi.h>
#include <d4est_hp_amr.h>
#include <d4est_hp_amr_smooth_pred.h>
#include <d4est_geometry.h>
#include <d4est_geometry_disk.h>
#include <d4est_vtk.h>
#include <d4est_mesh.h>
#include <ini.h>
#include <d4est_element_data.h>
#include <d4est_poisson.h>
#include <d4est_poisson_flux_sipg.h>
#include <d4est_cg.h>
#include <util.h>
#include <time.h>



#define NASTY_DEBUG

/* soon to be in the input files */
static const double pi = 3.1415926535897932384626433832795;
d4est_geometry_disk_attr_t global_disk_attrs;
d4est_poisson_flux_sipg_params_t global_ip_flux_params;


typedef struct {

  double* u;
  double* u_analytic;
  double* error;
  double* jacobian;

} vtk_nodal_vecs_t;


void
vtk_field_plotter
(
 d4est_vtk_context_t* vtk_ctx,
 void* user
)
{
 vtk_nodal_vecs_t* vecs = user;
  vtk_ctx = d4est_vtk_write_dg_point_dataf(vtk_ctx,
                                           4,
                                           0,
                                           "u",
                                           vecs->u,
                                           "u_analytic",
                                           vecs->u_analytic,
                                           "error",
                                           vecs->error,
                                           "jacobian",
                                           vecs->jacobian,
                                           vtk_ctx
                                          );


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
  
}


static int
uni_refine_function
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t *quadrant
)
{
  return 1;
}

typedef struct {

  int num_unifrefs;
  int deg;
  int deg_quad;
  int deg_stiffness;
  int rhs_use_lobatto;
  
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
  if (util_match_couple(section,"problem",name,"num_unifrefs")) {
    D4EST_ASSERT(pconfig->num_unifrefs == -1);
    pconfig->num_unifrefs = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg")) {
    D4EST_ASSERT(pconfig->deg == -1);
    pconfig->deg = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_quad")) {
    D4EST_ASSERT(pconfig->deg_quad == -1);
    pconfig->deg_quad = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_stiffness")) {
    D4EST_ASSERT(pconfig->deg_stiffness == -1);
    pconfig->deg_stiffness = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"rhs_use_lobatto")) {
    D4EST_ASSERT(pconfig->rhs_use_lobatto == -1);
    pconfig->rhs_use_lobatto = atoi(value);
  }  
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


static
problem_input_t
problem_input
(
 const char* input_file
)
{
  
  problem_input_t input;
  input.num_unifrefs = -1;
  input.deg = -1;
  input.deg_quad = -1; 
  input.deg_stiffness = -1; 
  input.rhs_use_lobatto = -1; 

  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("problem", input.num_unifrefs, -1);
  D4EST_CHECK_INPUT("problem", input.deg, -1);
  D4EST_CHECK_INPUT("problem", input.deg_quad, -1);

  return input;
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
     sizeof(d4est_element_data_t),
     NULL,
     NULL
    );
}



void
problem_set_degrees
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  // d4est_element_data_t* elem_data = elem_data_tmp;
  problem_input_t* input = user_ctx;
  /* outer shell */
  elem_data->deg = input->deg;
  elem_data->deg_quad = input->deg_quad;
}

static
double helmholtz_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 double u,
 void* user
)
{
  /* return 2.*x*y*z*x*y; */
  return 0.;
}

static
double analytic_solution_fcn
(
 double x,
 double y
 #if (P4EST_DIM)==3
 ,
 double z
#endif
)
{
  double r2 = x*x + y*y;
  return 1./sqrt(r2);
}

static
double boundary_fcn
(
 double x,
 double y
#if (P4EST_DIM)==3
 ,
 double z
#endif
)
{
  return analytic_solution_fcn(x,y
#if (P4EST_DIM)==3
                               ,z
#endif
                              );
}

static
double f_fcn
(
 double x,
 double y
#if (P4EST_DIM)==3
 ,
 double z
#endif
)
{
  double r2 = x*x + y*y;
  double r3 = pow(r2, 1.5);
  /* double u = analytic_solution_fcn(x,y); */
  return -(1./r3);// + helmholtz_fcn(x,y,u,NULL)*u;
}


static
double f_fcn_ext
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 double u0,
 void* user
)
{
  double r2 = x*x + y*y;
  double r3 = pow(r2, 1.5);
  double u = analytic_solution_fcn(x,y);
  return -(1./r3) + helmholtz_fcn(x,y,u,NULL)*u;
}

static
void problem_apply_lhs
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_mortar_fcn_ptrs_t* flux_fcn_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* user
)
{  
  d4est_poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, flux_fcn_data, d4est_ops, d4est_geom, d4est_quad);
}

static
void problem_build_rhs
(
 p4est_t* p4est,
 d4est_elliptic_data_t* prob_vecs,
 d4est_elliptic_eqns_t* prob_fcns,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_poisson_flux_sipg_params_t* sipg_params
)
{
  double* f = P4EST_ALLOC(double, prob_vecs->local_nodes);
  d4est_mesh_init_field
    (
     p4est,
     f,
     f_fcn,
     d4est_ops,
     d4est_geom
    );

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

        /* if (util_match(input->rhs_compute_method,"COMPUTE_RHS_ON_LOBATTO")){ */
        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq = ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = ed->q[2];
#endif
        
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
           ed->J_quad,
           ed->deg_quad,
           &prob_vecs->rhs[ed->nodal_stride]
          );
        
        printf("elem_id, rhs sum = %d %.25f\n", ed->id, d4est_linalg_vec_sum(&prob_vecs->rhs[ed->nodal_stride], d4est_lgl_get_nodes((P4EST_DIM), ed->deg)));
      }
    }
   
  int local_nodes = prob_vecs->local_nodes;
  double* u_eq_0 = P4EST_ALLOC_ZERO(double, local_nodes);
  double* tmp = prob_vecs->u;

  d4est_mortar_fcn_ptrs_t flux_fcn_data = d4est_poisson_flux_sipg_fetch_fcns(boundary_fcn, sipg_params);
  
  prob_vecs->u = u_eq_0;
  problem_apply_lhs(p4est, ghost, ghost_data, prob_vecs, &flux_fcn_data, d4est_ops, d4est_geom, d4est_quad, NULL);
  d4est_linalg_vec_axpy(-1., prob_vecs->Au, prob_vecs->rhs, local_nodes);

  prob_vecs->u = tmp;
  P4EST_FREE(u_eq_0);
  
  P4EST_FREE(f);
}

static
void
build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_mortar_fcn_ptrs_t* flux_fcn_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* user
)
{
  problem_apply_lhs(p4est, ghost, ghost_data, prob_vecs, flux_fcn_data, d4est_ops, d4est_geom, d4est_quad, user);
  d4est_linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);
}


void
problem_init
(
 const char* input_file,
 p4est_t* p4est,
 d4est_geometry_t* d4est_geom,
 d4est_operators_t* d4est_ops,
 int proc_size,
 sc_MPI_Comm mpicomm
)
{
  problem_input_t input = problem_input(input_file);
  global_disk_attrs = *((d4est_geometry_disk_attr_t*)d4est_geom->user);

  
  int level;
  
  D4EST_ASSERT((P4EST_DIM) == 2 || (P4EST_DIM) == 3);
  int world_rank, world_size;
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &world_rank);
  sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &world_size);

  double* Au = P4EST_ALLOC_ZERO(double, 1);
  double* rhs = P4EST_ALLOC_ZERO(double, 1);
  double* u = P4EST_ALLOC_ZERO(double, 1);
  double* u_analytic = P4EST_ALLOC_ZERO(double, 1);
  double* error = P4EST_ALLOC_ZERO(double, 1);
  int local_nodes = 1;

  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);

  p4est_partition(p4est, 0, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FACE, NULL);
  
  d4est_elliptic_data_t prob_vecs;
  prob_vecs.rhs = rhs;
  prob_vecs.Au = Au;
  prob_vecs.u = u;
  prob_vecs.local_nodes = local_nodes;

  d4est_poisson_flux_sipg_params_t* sipg_params = d4est_poisson_flux_sipg_params_new(p4est, "[SIPG_FLUX]", input_file);
  d4est_mortar_fcn_ptrs_t flux_fcn_data = d4est_poisson_flux_sipg_fetch_fcns(zero_fcn, sipg_params);


  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = build_residual;
  prob_fcns.apply_lhs = problem_apply_lhs;
  prob_fcns.flux_fcn_data = &flux_fcn_data;
  
  d4est_mesh_geometry_storage_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, input_file, "quadrature", "[QUADRATURE]");


  for (level = 0; level < input.num_unifrefs; ++level){

    if (world_rank == 0)
      printf("[D4EST_INFO]: UNIFORM REFINEMENT LEVEL %d\n", level);
    
      p4est_refine_ext(p4est,
                       0,
                       -1,
                       uni_refine_function,
                       NULL,
                       NULL
                      );

      p4est_partition(p4est, 0, NULL);
      p4est_balance_ext
        (
         p4est,
         P4EST_CONNECT_FACE,
         NULL,
         NULL
        );

      p4est_ghost_destroy(ghost);
      P4EST_FREE(ghost_data);

      ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
      ghost_data = P4EST_ALLOC(d4est_element_data_t, ghost->ghosts.elem_count);

     
  }

  d4est_mesh_update
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     geometric_factors,
     INITIALIZE_QUADRATURE_DATA,
     INITIALIZE_GEOMETRY_DATA,
     INITIALIZE_GEOMETRY_ALIASES,
     problem_set_degrees,
     (void*)&input
    );
    
    local_nodes = d4est_mesh_get_local_nodes(p4est);

    Au = P4EST_REALLOC(Au, double, local_nodes);
    u = P4EST_REALLOC(u, double, local_nodes);
    rhs = P4EST_REALLOC(rhs, double, local_nodes);
    u_analytic = P4EST_REALLOC(u_analytic, double, local_nodes);
    error = P4EST_REALLOC(error, double, local_nodes);

    prob_vecs.Au = Au;
    prob_vecs.u = u;
    prob_vecs.rhs = rhs;
    prob_vecs.local_nodes = local_nodes;

    d4est_mesh_init_field
      (
       p4est,
       u,
       analytic_solution_fcn,
       d4est_ops,
       d4est_geom
      );  
    
    d4est_mesh_init_field
      (
       p4est,
       u_analytic,
       analytic_solution_fcn,
       d4est_ops,
       d4est_geom
      );  
    d4est_linalg_vec_axpyeqz(-1., u, u_analytic, error, local_nodes);

    double initial_l2_norm_sqr_local = d4est_mesh_compute_l2_norm_sqr
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       error,
       local_nodes,
       DO_NOT_STORE_LOCALLY
      );
    printf("Initial local l2 norm = %.25f\n", sqrt(initial_l2_norm_sqr_local));


    problem_build_rhs
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       ghost,
       ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       sipg_params
      );
    /* DEBUG_PRINT_2ARR_DBL(prob_vecs.u, prob_vecs.rhs, local_nodes); */
    
    clock_t begin = 0;
    clock_t end = -1;

    if (world_rank == 0){
      begin = clock();
    }


   d4est_cg_params_t cg_params;
    d4est_cg_input(p4est, input_file, "d4est_cg", "[D4EST_CG]", &cg_params);


    DEBUG_PRINT_2ARR_DBL(prob_vecs.u,
                         prob_vecs.rhs,
                         prob_vecs.local_nodes);
    
    d4est_cg_solve
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       &ghost,
       (void**)&ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &cg_params
      );
    
    /* } */
    
    d4est_mesh_init_field
      (
       p4est,
       u_analytic,
       analytic_solution_fcn,
       d4est_ops,
       d4est_geom
      );
    d4est_linalg_vec_axpyeqz(-1., u, u_analytic, error, local_nodes);

    
    double local_l2_norm_sqr = d4est_mesh_compute_l2_norm_sqr
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       error,
       local_nodes,
       DO_NOT_STORE_LOCALLY
      );

    double* jacobian_lgl = P4EST_ALLOC(double, local_nodes);
    d4est_mesh_compute_jacobian_on_lgl_grid(p4est, d4est_ops, d4est_geom, jacobian_lgl);


    /* DEBUG_PRINT_ARR_DBL(jacobian_lgl, local_nodes); */
    
    d4est_linalg_vec_fabs(error, local_nodes); /* not needed, except for vtk output */
    vtk_nodal_vecs_t vtk_nodal_vecs;
    vtk_nodal_vecs.u = u;
    vtk_nodal_vecs.u_analytic = u_analytic;
    vtk_nodal_vecs.error = error;
    vtk_nodal_vecs.jacobian = jacobian_lgl;

    int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
    d4est_mesh_get_array_of_degrees(p4est, deg_array);

    char save_as [500];
    sprintf(save_as, "%s_level_%d", "test_uniform_cubedsphere", input.num_unifrefs);
    
    /* vtk output */
    d4est_vtk_save_geometry_and_dg_fields
      (
       save_as,
       p4est,
       d4est_ops,
       deg_array,
       input_file,
       "d4est_vtk_geometry",
       vtk_field_plotter,
       (void*)&vtk_nodal_vecs
      );

    P4EST_FREE(jacobian_lgl);
    P4EST_FREE(deg_array);
    double local_nodes_dbl = (double)local_nodes;
    double local_reduce [2];
    double global_reduce [2];

    local_reduce[0] = local_l2_norm_sqr;
    local_reduce[1] = local_nodes_dbl;
    
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

    double global_l2_norm_sqr = global_reduce[0];
    double global_nodes_dbl = global_reduce[1];

    if (world_rank == 0){
      end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf
      (
       "[HP_AMR]: %d, %d, %d, %.25f, %f\n",
       level,
       (int)p4est->global_num_quadrants,
       (int)global_nodes_dbl,
       sqrt(global_l2_norm_sqr),
       time_spent
      );
    }

  d4est_mesh_geometry_storage_destroy(geometric_factors);
  d4est_poisson_flux_sipg_params_destroy(sipg_params);
  
  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }

  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  
  P4EST_FREE(Au);
  P4EST_FREE(rhs);
  P4EST_FREE(u_analytic);
  P4EST_FREE(error);
  P4EST_FREE(u);
}

/*  */
