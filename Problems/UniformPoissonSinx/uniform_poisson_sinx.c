#include <sc_reduce.h>
#include <pXest.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <problem.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_estimator_bi.h>
#include <d4est_solver_cg.h>
#include <d4est_amr.h>
#include <d4est_amr_smooth_pred.h>
#include <d4est_geometry.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_vtk.h>
#include <d4est_output.h>
#include <d4est_mesh.h>
#include <ini.h>
#include <d4est_element_data.h>
#include <d4est_poisson.h>
#include <d4est_poisson_flux_sipg.h>
#include <newton_petsc.h>
#include <krylov_petsc.h>
#include <d4est_util.h>
#include <time.h>

typedef struct {

  int deg;
  int deg_quad;
  
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
  if (d4est_util_match_couple(section,"problem",name,"deg")) {
    D4EST_ASSERT(pconfig->deg == -1);
    pconfig->deg = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_quad")) {
    D4EST_ASSERT(pconfig->deg_quad == -1);
    pconfig->deg_quad = atoi(value);
  }
  else {
    return 0;
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
  input.deg = -1;
  input.deg_quad = -1;
  
  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("problem", input.deg, -1);
  D4EST_CHECK_INPUT("problem", input.deg_quad, -1);
  printf("[PROBLEM]: deg = %d\n",input.deg);
  printf("[PROBLEM]: deg_quad = %d\n",input.deg_quad);
  return input;
}


void
problem_set_degrees_init
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  problem_input_t* input = user_ctx;
  elem_data->deg = input->deg;
  elem_data->deg_quad = input->deg_quad;
}


void
problem_set_degrees_after_amr
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  elem_data->deg_quad = elem_data->deg;
}


double
problem_analytic_solution
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
#if (P4EST_DIM)==3
  return sin((M_PI)*x)*sin((M_PI)*y)*sin((M_PI)*z);
#else
  return sin((M_PI)*x)*sin((M_PI)*y);
#endif
}
double
problem_f_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return (M_PI)*(M_PI)*(P4EST_DIM)*problem_analytic_solution(x,y,z,user);
}

double
problem_initial_guess
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return 0.;
}

double
problem_boundary_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return problem_analytic_solution(x,y,z,user);
}

void
problem_apply_lhs
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* user
)
{
  d4est_poisson_flux_data_t* flux_fcn_data = user;
  d4est_poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, flux_fcn_data, d4est_ops, d4est_geom, d4est_quad);
}


void
problem_build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* user
)
{
  problem_apply_lhs(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom, d4est_quad, user);
  d4est_linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);
}

void
problem_build_rhs
(
 p4est_t* p4est,
 d4est_elliptic_data_t* prob_vecs,
 d4est_elliptic_eqns_t* prob_fcns,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_poisson_flux_data_t* flux_fcn_data_for_build_rhs
)
{
  double* f = P4EST_ALLOC(double, prob_vecs->local_nodes);
  d4est_mesh_init_field
    (
     p4est,
     f,
     problem_f_fcn,
     d4est_ops,
     d4est_geom,
     NULL
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
        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq = ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
        mesh_object.q[2] = ed->q[2];
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
      }
    }    

  int local_nodes = prob_vecs->local_nodes;
  double* u_eq_0 = P4EST_ALLOC_ZERO(double, local_nodes);
  double* tmp = prob_vecs->u; 
  prob_vecs->u = u_eq_0; 
  d4est_poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, flux_fcn_data_for_build_rhs, d4est_ops, d4est_geom, d4est_quad);
  d4est_linalg_vec_axpy(-1., prob_vecs->Au, prob_vecs->rhs, local_nodes);  
  prob_vecs->u = tmp;
  P4EST_FREE(u_eq_0);
  P4EST_FREE(f);
}

void
problem_init
(
 p4est_t* p4est,
 p4est_ghost_t** ghost,
 d4est_element_data_t** ghost_data,
 d4est_geometry_t* d4est_geom,
 d4est_operators_t* d4est_ops,
 const char* input_file,
 sc_MPI_Comm mpicomm
)
{
  D4EST_ASSERT(d4est_geom->geom_type == GEOM_BRICK);
  problem_input_t input = problem_input(input_file);

  d4est_mesh_geometry_storage_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, input_file, "quadrature", "[QUADRATURE]");
  d4est_poisson_flux_data_t* flux_data_for_apply_lhs = d4est_poisson_flux_new(p4est, input_file, zero_fcn, NULL);
  d4est_poisson_flux_data_t* flux_data_for_build_rhs = d4est_poisson_flux_new(p4est, input_file, problem_boundary_fcn, NULL);
  
  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = problem_build_residual;
  prob_fcns.apply_lhs = problem_apply_lhs;
  prob_fcns.user = flux_data_for_apply_lhs;
  
  double* error = NULL;
  double* u_analytic = NULL;
  int local_nodes = 0;
  
  d4est_elliptic_data_t prob_vecs;
  prob_vecs.Au = NULL;
  prob_vecs.u = NULL;
  prob_vecs.rhs = NULL;
  prob_vecs.local_nodes = 0;

  d4est_poisson_flux_sipg_params_t* sipg_params = flux_data_for_apply_lhs->user;
  d4est_estimator_bi_penalty_data_t penalty_data;
  penalty_data.u_penalty_fcn = bi_u_prefactor_conforming_maxp_minh;
  penalty_data.u_dirichlet_penalty_fcn = bi_u_prefactor_conforming_maxp_minh;
  penalty_data.gradu_penalty_fcn = bi_gradu_prefactor_maxp_minh;
  penalty_data.penalty_prefactor = sipg_params->sipg_penalty_prefactor;
  penalty_data.sipg_flux_h = sipg_params->sipg_flux_h;

  d4est_amr_t* d4est_amr =
    d4est_amr_init
    (
     p4est,
     input_file,
     "[D4EST_AMR]:",
     NULL
    );
  
  for (int level = 0; level < d4est_amr->num_of_amr_steps + 1; ++level){

    local_nodes = d4est_mesh_update
                  (
                   p4est,
                   *ghost,
                   *ghost_data,
                   d4est_ops,
                   d4est_geom,
                   d4est_quad,
                   geometric_factors,
                   INITIALIZE_QUADRATURE_DATA,
                   INITIALIZE_GEOMETRY_DATA,
                   INITIALIZE_GEOMETRY_ALIASES,
                   (level == 0) ? problem_set_degrees_init : problem_set_degrees_after_amr,
                   (void*)&input
                  );

    if (level == 0){
      prob_vecs.u = P4EST_REALLOC(prob_vecs.u, double, local_nodes);
      d4est_mesh_init_field
        (
         p4est,
         prob_vecs.u,
         problem_initial_guess,
         d4est_ops,
         d4est_geom,
         NULL
        );
    }
    
    prob_vecs.Au = P4EST_REALLOC(prob_vecs.Au, double, local_nodes);
    prob_vecs.rhs = P4EST_REALLOC(prob_vecs.rhs, double, local_nodes);
    prob_vecs.local_nodes = local_nodes;
    
    problem_build_rhs
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       *ghost,
       *ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       flux_data_for_build_rhs
      );
    
    d4est_estimator_bi_compute
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       penalty_data,
       zero_fcn,
       *ghost,
       *ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad
      );

    d4est_estimator_stats_t stats;
    d4est_estimator_stats_compute(p4est, &stats);
    d4est_estimator_stats_print(&stats);
    
    d4est_solver_cg_params_t params;
    d4est_solver_cg_input
      (
       p4est,
       input_file,
       "d4est_solver_cg",
       "[D4EST_SOLVER_CG]",
       &params
      );
    
    d4est_solver_cg_solve
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       ghost,
       ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &params
      );

    d4est_output_vtk_with_analytic_error
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &prob_vecs,
       input_file,
       "uniform_poisson_sinx",
       problem_analytic_solution,
       NULL,
       OUTPUT_ESTIMATOR,
       level
      );
    
    d4est_output_norms_using_analytic_solution
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &stats,
       &prob_vecs,
       problem_analytic_solution,
       NULL,
       level);

    if (level != d4est_amr->num_of_amr_steps){

      if (p4est->mpirank == 0)
        printf("[D4EST_INFO]: AMR REFINEMENT LEVEL %d\n", level+1);

      d4est_amr_step
        (
         p4est,
         ghost,
         ghost_data,
         d4est_ops,
         d4est_amr,
         &prob_vecs.u,
         NULL
        );
      
    }  
    
  }

  printf("[D4EST_INFO]: Starting garbage collection...\n");


  d4est_amr_destroy(d4est_amr);
  d4est_mesh_geometry_storage_destroy(geometric_factors);
  d4est_poisson_flux_destroy(flux_data_for_apply_lhs);  
  d4est_poisson_flux_destroy(flux_data_for_build_rhs);  
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);

  P4EST_FREE(error);
  P4EST_FREE(u_analytic);
}
