#include <sc_reduce.h>
#include <pXest.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <d4est_element_data.h>
#include <sipg_flux_vector_fcns.h>
#include <problem.h>
#include <problem_data.h>
#include <d4est_elliptic_eqns.h>
#include <central_flux_params.h>
#include <curved_bi_estimator.h>
#include <krylov_petsc.h>
#include <matrix_sym_tester.h>
#include <dg_norm.h>
#include <hp_amr.h>
#include <hp_amr_curved_smooth_pred.h>
#include <hp_amr_curved_uniform.h>
#include <d4est_geometry.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_geometry_disk.h>
#include <curved_poisson_debug_vecs.h>
#include <d4est_vtk.h>
#include <bi_estimator_flux_fcns.h>
#include <newton_petsc.h>
#include <ini.h>
#include <d4est_poisson.h>
#include <curved_Gauss_primal_sipg_kronbichler_flux_fcns.h>
#include <multigrid_matrix_operator.h>
#include <multigrid_smoother_cheby_d4est.h>
#include <multigrid_smoother_krylov_petsc.h>
#include <multigrid_bottom_solver_cg_d4est.h>
#include <multigrid_bottom_solver_cheby_d4est.h>
#include <multigrid_bottom_solver_krylov_petsc.h>
#include <multigrid_logger_residual.h>
#include <multigrid_element_data_updater_curved.h>
#include <krylov_pc_multigrid.h>
#include <ip_flux_params.h>
#include <ip_flux.h>
#include "time.h"
#include "d4est_util.h"

d4est_geometry_cubed_sphere_attr_t global_cubed_sphere_attrs;
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
  int deg_R0;
  int deg_quad_R0;
  int deg_stiffness_R0;
  int deg_R1;
  int deg_quad_R1;
  int deg_stiffness_R1;
  int deg_R2;
  int deg_quad_R2;
  int deg_stiffness_R2;
  char rhs_compute_method [50];
  
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
  if (d4est_util_match_couple(section,"problem",name,"num_unifrefs")) {
    D4EST_ASSERT(pconfig->num_unifrefs == -1);
    pconfig->num_unifrefs = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"rhs_compute_method")) {
    D4EST_ASSERT(pconfig->rhs_compute_method[0] == '*');
    snprintf (pconfig->rhs_compute_method, sizeof(pconfig->rhs_compute_method), "%s", value);
    D4EST_ASSERT(d4est_util_match(pconfig->rhs_compute_method, "COMPUTE_RHS_ON_GAUSS") ||
               d4est_util_match(pconfig->rhs_compute_method, "COMPUTE_RHS_ON_LOBATTO") );
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_R0")) {
    D4EST_ASSERT(pconfig->deg_R0 == -1);
    pconfig->deg_R0 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_quad_R0")) {
    D4EST_ASSERT(pconfig->deg_quad_R0 == -1);
    pconfig->deg_quad_R0 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_stiffness_R0")) {
    D4EST_ASSERT(pconfig->deg_stiffness_R0 == -1);
    pconfig->deg_stiffness_R0 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_R1")) {
    D4EST_ASSERT(pconfig->deg_R1 == -1);
    pconfig->deg_R1 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_quad_R1")) {
    D4EST_ASSERT(pconfig->deg_quad_R1 == -1);
    pconfig->deg_quad_R1 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_stiffness_R1")) {
    D4EST_ASSERT(pconfig->deg_stiffness_R1 == -1);
    pconfig->deg_stiffness_R1 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_R2")) {
    D4EST_ASSERT(pconfig->deg_R2 == -1);
    pconfig->deg_R2 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_quad_R2")) {
    D4EST_ASSERT(pconfig->deg_quad_R2 == -1);
    pconfig->deg_quad_R2 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_stiffness_R2")) {
    D4EST_ASSERT(pconfig->deg_stiffness_R2 == -1);
    pconfig->deg_stiffness_R2 = atoi(value);
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
  input.deg_R0 = -1;
  input.deg_quad_R0 = -1;
  input.deg_stiffness_R0 = -1;
  input.deg_R1 = -1;
  input.deg_quad_R1 = -1;
  input.deg_stiffness_R1 = -1;
  input.deg_R2 = -1;
  input.deg_quad_R2 = -1; 
  input.deg_stiffness_R2 = -1; 
  input.rhs_compute_method[0] = '*'; 

  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("problem", input.num_unifrefs, -1);
  D4EST_CHECK_INPUT("problem", input.deg_R0, -1);
  D4EST_CHECK_INPUT("problem", input.deg_quad_R0, -1);
  D4EST_CHECK_INPUT("problem", input.deg_stiffness_R0, -1);
  D4EST_CHECK_INPUT("problem", input.deg_R1, -1);
  D4EST_CHECK_INPUT("problem", input.deg_quad_R1, -1);
  D4EST_CHECK_INPUT("problem", input.deg_stiffness_R1, -1);
  D4EST_CHECK_INPUT("problem", input.deg_R2, -1);
  D4EST_CHECK_INPUT("problem", input.deg_quad_R2, -1);
  D4EST_CHECK_INPUT("problem", input.deg_stiffness_R2, -1);  
  D4EST_CHECK_INPUT("problem", input.rhs_compute_method[0], '*');  

  printf("[PROBLEM]: test_uniform_cubedsphere\n");
  printf("[PROBLEM]: rhs_compute_method = %s\n", input.rhs_compute_method);
  
  
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
 void* elem_data_tmp,
 void* user_ctx
)
{
  d4est_element_data_t* elem_data = elem_data_tmp;
  problem_input_t* input = user_ctx;
  /* outer shell */
  if (elem_data->tree < 6){
    elem_data->deg = input->deg_R2;
    elem_data->deg_quad = input->deg_quad_R2;
    elem_data->deg_stiffness = input->deg_stiffness_R2;
  }
  /* inner shell */
  else if(elem_data->tree < 12){
    elem_data->deg = input->deg_R1;
    elem_data->deg_quad = input->deg_quad_R1;
    elem_data->deg_stiffness = input->deg_stiffness_R1;
  }
  /* center cube */
  else {
    elem_data->deg = input->deg_R0;
    elem_data->deg_quad = input->deg_quad_R0;
    elem_data->deg_stiffness = input->deg_stiffness_R0;
  } 
}



static
double helmholtz_fcn
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  /* return 2.*x*y*z*x*y; */
  return 1.;
}

static
double analytic_solution_fcn
(
 double x,
 double y,
 double z
)
{
  double R2 = global_cubed_sphere_attrs.R2;
  double r2 = x*x + y*y + z*z;
  return -r2;
}

static
double boundary_fcn
(
 double x,
 double y,
 double z
)
{
  /* printf("x, y, z, boundary_fcn = %f, %f, %f, %f\n", x,y,z,analytic_solution_fcn(x,y,z)); */
  return analytic_solution_fcn(x,y,z);
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
  return 6 + helmholtz_fcn(x,y,z,0,NULL)*analytic_solution_fcn(x,y,z);
}

static
double f_fcn_ext
(
 double x,
 double y,
 double z,
 double u0,
 void* user
)
{
  return 6 + helmholtz_fcn(x,y,z,0,NULL)*analytic_solution_fcn(x,y,z);
}

static
void problem_apply_lhs
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
  /* prob_vecs->flux_fcn_data = curved_gauss_primal_sipg_kronbichler_flux_dirichlet_fetch_fcns */
                                             /* (zero_fcn, &global_ip_flux_params); */

  d4est_poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom);

 double* M_helmf_u = P4EST_ALLOC(double, prob_vecs->local_nodes);
 

 int nodal_stride = 0;
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
       int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);

       d4est_element_data_apply_fofufofvlilj_gaussnodes
         (
          d4est_ops,
          d4est_geom,
          &prob_vecs->u[nodal_stride],
          NULL,
          NULL,
          ed,
          ed->deg_stiffness, // + params->deg_offset_for_nonlinear_quad,
          (P4EST_DIM),
          &M_helmf_u[nodal_stride],           
          helmholtz_fcn,
          prob_vecs->user,
          NULL,
          NULL
         );

       nodal_stride += volume_nodes_lobatto;
     }
   }
  
  d4est_linalg_vec_axpy(1.0, M_helmf_u, prob_vecs->Au, prob_vecs->local_nodes);
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
 problem_input_t* input,
 void* user
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

  prob_vecs->flux_fcn_data = curved_gauss_primal_sipg_kronbichler_flux_dirichlet_fetch_fcns
                                             (boundary_fcn, &global_ip_flux_params);

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

        if (d4est_util_match(input->rhs_compute_method,"COMPUTE_RHS_ON_LOBATTO")){
        d4est_operators_apply_curvedgaussMass(d4est_ops,
                                     &f[ed->nodal_stride],
                                     ed->deg,
                                     ed->J_quad,
                                     ed->deg_quad,
                                     (P4EST_DIM),
                                     &prob_vecs->rhs[ed->nodal_stride]
                                    );
        }
        else if(d4est_util_match(input->rhs_compute_method,"COMPUTE_RHS_ON_GAUSS")){
          d4est_element_data_apply_fofufofvlj_gaussnodes
            (
             d4est_ops,
             d4est_geom,
             NULL,
             NULL,
             ed,
             ed->deg_stiffness,
             (P4EST_DIM),
             &prob_vecs->rhs[ed->nodal_stride],
             f_fcn_ext,
             NULL,
             NULL,
             NULL
            );
        }
        else {
          D4EST_ABORT("Should not happen\n");
        }


        printf("elem_id, rhs sum = %d %.25f\n", ed->id, d4est_linalg_vec_sum(&prob_vecs->rhs[ed->nodal_stride], d4est_lgl_get_nodes((P4EST_DIM), ed->deg)));
        
        
      }
    }    

  int local_nodes = prob_vecs->local_nodes;
  double* u_eq_0 = P4EST_ALLOC_ZERO(double, local_nodes);
  double* tmp = prob_vecs->u;
  
  prob_vecs->u = u_eq_0; 
  problem_apply_lhs(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom);
  /* DEBUG_PRINT_ARR_DBL(prob_vecs->Au, prob_vecs->local_nodes); */
  d4est_linalg_vec_axpy(-1., prob_vecs->Au, prob_vecs->rhs, local_nodes);
  
  prob_vecs->u = tmp;
  P4EST_FREE(u_eq_0);

  prob_vecs->flux_fcn_data = curved_gauss_primal_sipg_kronbichler_flux_dirichlet_fetch_fcns
                                             (zero_fcn, &global_ip_flux_params);
  
  P4EST_FREE(f);
}


static
void
problem_build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
  problem_apply_lhs(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom);
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
  global_cubed_sphere_attrs = *((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user);
  
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

  ip_flux_t* ip_flux = ip_flux_dirichlet_new(p4est, "[IP_FLUX]", input_file, zero_fcn);
  prob_vecs.flux_fcn_data = ip_flux->mortar_fcn_ptrs;
  global_ip_flux_params = *(ip_flux->ip_flux_params);
  prob_vecs.user = &input;

  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = problem_build_residual;
  prob_fcns.apply_lhs = problem_apply_lhs;
  
  d4est_mesh_geometry_storage_t* geometric_factors = geometric_factors_init(p4est);


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




      d4est_element_data_init_new(p4est,
                                   geometric_factors,
                                   d4est_ops,
                                   d4est_geom,
                                   problem_set_degrees,
                                   (void*)&input, 1, 1);
    
      local_nodes = d4est_element_data_get_local_nodes(p4est);

      Au = P4EST_REALLOC(Au, double, local_nodes);
      u = P4EST_REALLOC(u, double, local_nodes);
      rhs = P4EST_REALLOC(rhs, double, local_nodes);
      prob_vecs.Au = Au;
      prob_vecs.u = u;
      prob_vecs.rhs = rhs;
      prob_vecs.local_nodes = local_nodes;
      
  }


  d4est_element_data_init_new(p4est,
                               geometric_factors,
                               d4est_ops,
                               d4est_geom,
                               problem_set_degrees,
                               (void*)&input,
                               1,
                               1);


    
    local_nodes = d4est_element_data_get_local_nodes(p4est);

    Au = P4EST_REALLOC(Au, double, local_nodes);
    u = P4EST_REALLOC(u, double, local_nodes);
    rhs = P4EST_REALLOC(rhs, double, local_nodes);
    u_analytic = P4EST_REALLOC(u_analytic, double, local_nodes);
    error = P4EST_REALLOC(error, double, local_nodes);

    prob_vecs.Au = Au;
    prob_vecs.u = u;
    prob_vecs.rhs = rhs;
    prob_vecs.local_nodes = local_nodes;

    d4est_linalg_fill_vec(prob_vecs.u, 1., prob_vecs.local_nodes);

    d4est_mesh_init_field(
                                      p4est,
                                      u_analytic,
                                      analytic_solution_fcn,
                                      d4est_ops,
                                      d4est_geom
    );  
    d4est_linalg_vec_axpyeqz(-1., u, u_analytic, error, local_nodes);

    DEBUG_PRINT_ARR_DBL_SUM(error, local_nodes);
    DEBUG_PRINT_ARR_DBL_SUM(prob_vecs.u, local_nodes);
    DEBUG_PRINT_ARR_DBL_SUM(u_analytic, local_nodes);
    
    double initial_l2_norm_sqr_local = d4est_element_data_compute_l2_norm_sqr
      (
       p4est,
       error,
       local_nodes,
       d4est_ops,
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
       &input,
       NULL
      );   

    
    clock_t begin = 0;
    clock_t end = -1;

    if (world_rank == 0){
      begin = clock();
    }

    krylov_petsc_params_t petsc_params;
    krylov_petsc_input(p4est, input_file, "krylov_petsc", "[KRYLOV_PETSC]", &petsc_params);      
       
    krylov_info_t info =
      krylov_petsc_solve
      (
       p4est,
       &prob_vecs,
       (void*)&prob_fcns,
       &ghost,
       (void**)&ghost_data,
       d4est_ops,
       d4est_geom,
       &petsc_params,
       NULL
      );
    
    d4est_mesh_init_field(
                                      p4est,
                                      u_analytic,
                                      analytic_solution_fcn,
                                      d4est_ops,
                                      d4est_geom
                                     );  
    d4est_linalg_vec_axpyeqz(-1., u, u_analytic, error, local_nodes);

    
    double local_l2_norm_sqr = d4est_element_data_compute_l2_norm_sqr
                                (
                                 p4est,
                                 error,
                                 local_nodes,
                                 d4est_ops,
                                 DO_NOT_STORE_LOCALLY
                                );

    
    double* jacobian_lgl = P4EST_ALLOC(double, local_nodes);
    d4est_element_data_compute_jacobian_on_lgl_grid(p4est,d4est_geom, d4est_ops, jacobian_lgl);


    d4est_linalg_vec_fabs(error, local_nodes); /* not needed, except for vtk output */
    vtk_nodal_vecs_t vtk_nodal_vecs;
    vtk_nodal_vecs.u = u;
    vtk_nodal_vecs.u_analytic = u_analytic;
    vtk_nodal_vecs.error = error;
    vtk_nodal_vecs.jacobian = jacobian_lgl;

    int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
    d4est_element_data_get_array_of_degrees(p4est, deg_array);

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

  geometric_factors_destroy(geometric_factors);
  ip_flux_dirichlet_destroy(ip_flux);
  
  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }

  P4EST_FREE(Au);
  P4EST_FREE(rhs);
  P4EST_FREE(u_analytic);
  P4EST_FREE(error);
  P4EST_FREE(u);
}
