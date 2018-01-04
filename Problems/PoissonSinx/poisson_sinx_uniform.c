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
#include "poisson_sinx_fcns.h"

typedef struct {
  
  dirichlet_bndry_eval_method_t eval_method;
  
} poisson_sinx_init_params_t;

static
int skip_element_fcn
(
 d4est_element_data_t* ed
)
{
  if(ed->tree != 12)
    return 1;
  else
    return 0;
}

static
int poisson_sinx_init_params_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  poisson_sinx_init_params_t* pconfig = (poisson_sinx_init_params_t*)user;
  if (d4est_util_match_couple(section,"problem",name,"eval_method")) {
    if (d4est_util_match(value,"EVAL_BNDRY_FCN_ON_QUAD")){
      pconfig->eval_method = EVAL_BNDRY_FCN_ON_QUAD;
    }
    else if (d4est_util_match(value,"EVAL_BNDRY_FCN_ON_LOBATTO")){
      pconfig->eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;
    }
    else {
      D4EST_ABORT("Not a supported eval method");
    }
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


static
poisson_sinx_init_params_t
poisson_sinx_init_params_input
(
 const char* input_file
)
{
  poisson_sinx_init_params_t input;
  input.eval_method = EVAL_BNDRY_FCN_NOT_SET;

  if (ini_parse(input_file, poisson_sinx_init_params_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("problem", input.eval_method, EVAL_BNDRY_FCN_NOT_SET);  
  
  return input;
}

int
problem_set_mortar_degree
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  return elem_data->deg_vol_quad;
}


void
problem_init
(
 p4est_t* p4est,
 p4est_ghost_t** ghost,
 d4est_element_data_t** ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_geometry_storage_t* d4est_factors,
 d4est_mesh_initial_extents_t* initial_extents,
 const char* input_file,
 sc_MPI_Comm mpicomm
)
{

  int initial_nodes = initial_extents->initial_nodes;
  
  poisson_sinx_init_params_t init_params = poisson_sinx_init_params_input(input_file
                                                                           );
  dirichlet_bndry_eval_method_t eval_method = init_params.eval_method;
  
  d4est_poisson_dirichlet_bc_t bc_data_for_lhs;
  bc_data_for_lhs.dirichlet_fcn = zero_fcn;
  bc_data_for_lhs.eval_method = eval_method;
  
  d4est_poisson_dirichlet_bc_t bc_data_for_rhs;
  bc_data_for_rhs.dirichlet_fcn = poisson_sinx_boundary_fcn;
  bc_data_for_rhs.eval_method = eval_method;
  
  d4est_poisson_flux_data_t* flux_data_for_apply_lhs = d4est_poisson_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_for_lhs, problem_set_mortar_degree, NULL);
  
  d4est_poisson_flux_data_t* flux_data_for_build_rhs = d4est_poisson_flux_new(p4est, input_file,  BC_DIRICHLET, &bc_data_for_rhs, problem_set_mortar_degree, NULL);

  problem_ctx_t ctx;
  ctx.flux_data_for_apply_lhs = flux_data_for_apply_lhs;
  ctx.flux_data_for_build_rhs = flux_data_for_build_rhs;
                           
  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = poisson_sinx_build_residual;
  prob_fcns.apply_lhs = poisson_sinx_apply_lhs;
  prob_fcns.user = &ctx;
  
  double* error = NULL;
  double* u_analytic = NULL;
  
  d4est_elliptic_data_t prob_vecs;
  prob_vecs.Au = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.u = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.rhs = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.local_nodes = initial_nodes;

  d4est_poisson_flux_sipg_params_t* sipg_params = flux_data_for_apply_lhs->flux_data;
 
  d4est_amr_t* d4est_amr =
    d4est_amr_init
    (
     p4est,
     input_file,
     "[D4EST_AMR]:",
     NULL
    );

  D4EST_ASSERT(d4est_amr->scheme->amr_scheme_type == AMR_UNIFORM_H ||
               d4est_amr->scheme->amr_scheme_type == AMR_UNIFORM_P);

  d4est_mesh_init_field
    (
     p4est,
     prob_vecs.u,
     poisson_sinx_initial_guess,
     d4est_ops,
     d4est_geom,
     INIT_FIELD_ON_LOBATTO,
     NULL
    );
    
  d4est_poisson_build_rhs_with_strong_bc
    (
     p4est,
     *ghost,
     *ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &prob_vecs,
     flux_data_for_build_rhs,
     prob_vecs.rhs,
     poisson_sinx_rhs_fcn,
     (init_params.eval_method == EVAL_BNDRY_FCN_ON_QUAD) ? INIT_FIELD_ON_QUAD : INIT_FIELD_ON_LOBATTO,
     &ctx
    );

  
  for (int level = 0; level < d4est_amr->num_of_amr_steps + 1; ++level){


    krylov_petsc_params_t krylov_petsc_params;
    krylov_petsc_input(p4est, input_file, "krylov_petsc", "[KRYLOV_PETSC]", &krylov_petsc_params);
    
    krylov_petsc_solve
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       ghost,
       ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       &krylov_petsc_params,
       NULL
      );

    double* u_analytic = P4EST_ALLOC(double, prob_vecs.local_nodes);
    double* error = P4EST_ALLOC(double, prob_vecs.local_nodes);

    d4est_mesh_init_field
      (
       p4est,
       u_analytic,
       poisson_sinx_analytic_solution,
       d4est_ops,
       d4est_geom,
       INIT_FIELD_ON_LOBATTO,
       NULL
      );


    for (int i = 0; i < prob_vecs.local_nodes; i++){
      error[i] = fabs(prob_vecs.u[i] - u_analytic[i]);
    }
    
    d4est_vtk_save
      (
       p4est,
       d4est_ops,
       input_file,
       "d4est_vtk",
       "[D4EST_VTK]",
       (const char * []){"u","u_analytic","error", NULL},
       (double* []){prob_vecs.u, u_analytic, error},
       (const char * []){NULL},
       (double* []){NULL}
      );

    P4EST_FREE(u_analytic);
    P4EST_FREE(error);
    
    d4est_ip_energy_norm_data_t ip_norm_data;
    ip_norm_data.u_penalty_fcn = sipg_params->sipg_penalty_fcn;
    ip_norm_data.sipg_flux_h = sipg_params->sipg_flux_h;
    ip_norm_data.penalty_prefactor = sipg_params->sipg_penalty_prefactor;
    printf("ip_norm_data.penalty_prefactor = %f\n", ip_norm_data.penalty_prefactor);
    
    d4est_output_norms_using_analytic_solution
      (
      p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
      d4est_factors,
       *ghost,
       *ghost_data,
      -1.,
       &prob_vecs,
       &ip_norm_data,
       poisson_sinx_analytic_solution,
      &ctx,NULL,NULL);

    
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


    prob_vecs.local_nodes = d4est_mesh_update
                  (
                   p4est,
                   *ghost,
                   *ghost_data,
                   d4est_ops,
                   d4est_geom,
                   d4est_quad,
                   d4est_factors,
                   INITIALIZE_QUADRATURE_DATA,
                   INITIALIZE_GEOMETRY_DATA,
                   INITIALIZE_GEOMETRY_ALIASES,
                   d4est_mesh_set_quadratures_after_amr,
                   initial_extents
                  );

    
    prob_vecs.Au = P4EST_REALLOC(prob_vecs.Au, double, prob_vecs.local_nodes);
    prob_vecs.rhs = P4EST_REALLOC(prob_vecs.rhs, double, prob_vecs.local_nodes);
    
    
    d4est_poisson_build_rhs_with_strong_bc
      (
       p4est,
       *ghost,
       *ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       &prob_vecs,
       flux_data_for_build_rhs,
       prob_vecs.rhs,
       poisson_sinx_rhs_fcn,
       (init_params.eval_method == EVAL_BNDRY_FCN_ON_QUAD) ? INIT_FIELD_ON_QUAD : INIT_FIELD_ON_LOBATTO,
       &ctx
      );



  }

  printf("[D4EST_INFO]: Starting garbage collection...\n");
  d4est_amr_destroy(d4est_amr);
  d4est_poisson_flux_destroy(flux_data_for_apply_lhs);  
  d4est_poisson_flux_destroy(flux_data_for_build_rhs);  
  P4EST_FREE(error);
  P4EST_FREE(u_analytic);
  P4EST_FREE(prob_vecs.u);
  P4EST_FREE(prob_vecs.Au);
  P4EST_FREE(prob_vecs.rhs);
}
