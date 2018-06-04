#include <sc_reduce.h>
#include <pXest.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <problem.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_estimator_bi.h>
#include <d4est_solver_cg.h>
#include <d4est_solver_fcg.h>
#include <d4est_amr.h>
#include <d4est_amr_smooth_pred.h>
#include <d4est_geometry.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_vtk.h>
#include <d4est_norms.h>
#include <d4est_mesh.h>
#include <ini.h>
#include <d4est_element_data.h>
#include <d4est_poisson.h>
#include <d4est_poisson_flux_sipg.h>
#include <newton_petsc.h>
#include <krylov_petsc.h>
#include <krylov_pc_multigrid.h>
#include <d4est_power_method.h>
#include <d4est_solver_full_matrix.h>
#include <d4est_solver_matrix_symmetry.h>
#include <cg_eigs.h>
#include <multigrid_logger_residual.h>
#include <multigrid_element_data_updater.h>
#include <multigrid.h>
#include <d4est_util.h>
#include <time.h>
#include <zlog.h>
#include "poisson_lorentzian_fcns.h"


typedef struct {
  
  int use_dirichlet;
  int cg_eig_its;
  int pow_eig_its;
  
} poisson_lorentzian_init_params_t;



static
int poisson_lorentzian_init_params_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  poisson_lorentzian_init_params_t* pconfig = (poisson_lorentzian_init_params_t*)user;
  if (d4est_util_match_couple(section,"problem",name,"use_dirichlet")) {
    D4EST_ASSERT(pconfig->use_dirichlet == -1);
    pconfig->use_dirichlet = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"pow_eig_its")) {
    D4EST_ASSERT(pconfig->pow_eig_its == -1);
    pconfig->pow_eig_its = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"cg_eig_its")) {
    D4EST_ASSERT(pconfig->cg_eig_its == -1);
    pconfig->cg_eig_its = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

static
poisson_lorentzian_init_params_t
poisson_lorentzian_init_params_input
(
 const char* input_file
)
{
  poisson_lorentzian_init_params_t input;
  input.use_dirichlet = -1;
  input.cg_eig_its = -1;
  input.pow_eig_its = -1;

  if (ini_parse(input_file, poisson_lorentzian_init_params_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("problem", input.use_dirichlet, -1);
  D4EST_CHECK_INPUT("problem", input.cg_eig_its, -1);
  D4EST_CHECK_INPUT("problem", input.pow_eig_its, -1);
  
  return input;
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
 d4est_mesh_data_t* d4est_factors,
 d4est_mesh_initial_extents_t* initial_extents,
 const char* input_file,
 sc_MPI_Comm mpicomm
)
{
  zlog_category_t *c_default = zlog_get_category("problem");

  int initial_nodes = initial_extents->initial_nodes;



  poisson_lorentzian_init_params_t init_params = poisson_lorentzian_init_params_input
                                            (
                                             input_file
                                            ); 

  
  dirichlet_bndry_eval_method_t eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;
  
  // Setup boundary conditions

  d4est_poisson_robin_bc_t bc_data_robin_for_lhs;
  bc_data_robin_for_lhs.robin_coeff = poisson_lorentzian_robin_coeff_fcn;
  bc_data_robin_for_lhs.robin_rhs = poisson_lorentzian_robin_bc_rhs_fcn;
  
  d4est_poisson_dirichlet_bc_t bc_data_dirichlet_for_lhs;
  bc_data_dirichlet_for_lhs.dirichlet_fcn = zero_fcn;
  bc_data_dirichlet_for_lhs.eval_method = eval_method;
  
  d4est_poisson_dirichlet_bc_t bc_data_dirichlet_for_rhs;
  bc_data_dirichlet_for_rhs.dirichlet_fcn = poisson_lorentzian_boundary_fcn;
  bc_data_dirichlet_for_rhs.eval_method = eval_method;
  
  d4est_poisson_flux_data_t* flux_data_for_lhs = NULL;//d4est_poisson_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_dirichlet_for_lhs);
  
  d4est_poisson_flux_data_t* flux_data_for_rhs = NULL;//d4est_poisson_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_dirichlet_for_rhs);

  if(init_params.use_dirichlet){
    flux_data_for_lhs
      = d4est_poisson_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_dirichlet_for_lhs);
  
    flux_data_for_rhs
      = d4est_poisson_flux_new(p4est, input_file,  BC_DIRICHLET, &bc_data_dirichlet_for_rhs);
  }
  else {  
    flux_data_for_lhs = d4est_poisson_flux_new(p4est, input_file, BC_ROBIN, &bc_data_robin_for_lhs);
    flux_data_for_rhs = d4est_poisson_flux_new(p4est, input_file,  BC_ROBIN, &bc_data_robin_for_lhs);
  }
  

  problem_ctx_t ctx;
  ctx.flux_data_for_apply_lhs = flux_data_for_lhs;
  ctx.flux_data_for_build_rhs = flux_data_for_rhs;


  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = poisson_lorentzian_build_residual;
  prob_fcns.apply_lhs = poisson_lorentzian_apply_lhs;
  prob_fcns.user = &ctx;


  d4est_elliptic_data_t prob_vecs;
  prob_vecs.Au = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.u = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.rhs = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.local_nodes = initial_nodes;

  d4est_poisson_flux_sipg_params_t* sipg_params = flux_data_for_lhs->flux_data;
  
  
  // Setup norm function contexts
  
  d4est_norms_fcn_L2_ctx_t L2_norm_ctx;
  L2_norm_ctx.p4est = p4est;
  L2_norm_ctx.d4est_ops = d4est_ops;
  L2_norm_ctx.d4est_geom = d4est_geom;
  L2_norm_ctx.d4est_quad = d4est_quad;
  L2_norm_ctx.d4est_factors = d4est_factors;

  if (p4est->mpirank == 0)
    d4est_norms_write_headers(
      (const char * []){"u", NULL},
      (const char * []){"L_2", "L_infty", NULL}
    );


  // Setup AMR
  d4est_amr_t* d4est_amr = d4est_amr_init(
    p4est,
    input_file,
    NULL
  );

  D4EST_ASSERT(
    d4est_amr->scheme->amr_scheme_type == AMR_UNIFORM_H ||
    d4est_amr->scheme->amr_scheme_type == AMR_UNIFORM_P
  );

  d4est_mesh_init_field(
    p4est,
    prob_vecs.u,
    poisson_lorentzian_initial_guess,
    d4est_ops,
    d4est_geom,
    d4est_factors,
    INIT_FIELD_ON_LOBATTO,
    NULL
  );
    
  d4est_poisson_build_rhs_with_strong_bc(
    p4est,
    *ghost,
    *ghost_data,
    d4est_ops,
    d4est_geom,
    d4est_quad,
    d4est_factors,
    &prob_vecs,
    flux_data_for_rhs,
    prob_vecs.rhs,
    poisson_lorentzian_rhs_fcn,
    INIT_FIELD_ON_LOBATTO,
    &ctx
  );


  for (int level = 0; level < d4est_amr->num_of_amr_steps + 1; level++) {

    // Setup multigrid
    /* krylov_pc_t* pc = NULL; */

    /* multigrid_data_t* mg_data = multigrid_data_init( */
    /*   p4est, */
    /*   d4est_ops, */
    /*   d4est_geom, */
    /*   d4est_quad, */
    /*   ghost, */
    /*   ghost_data, */
    /*   d4est_factors, */
    /*   initial_extents, */
    /*   input_file */
    /* ); */

      
    /* pc = krylov_pc_multigrid_create(mg_data, NULL); */


    int cg_its = 30;
    int pow_its = 1000;
  
    double eig_power =
      d4est_power_method
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       *ghost,
       ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       init_params.pow_eig_its
      );

  
    double eig_cg;
    cg_eigs
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       *ghost,
       ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       init_params.cg_eig_its,
       &eig_cg
      );
      

      printf("eig_power = %.15f after %d its, eig_cg = %.15f after %d its\n", eig_power, init_params.pow_eig_its, eig_cg, init_params.cg_eig_its);
    
    // Krylov PETSc solve
    
    /* krylov_petsc_params_t krylov_petsc_params; */
    /* krylov_petsc_input(p4est, input_file, "krylov_petsc", &krylov_petsc_params); */
    
    /* krylov_petsc_solve( */
    /*   p4est, */
    /*   &prob_vecs, */
    /*   &prob_fcns, */
    /*   ghost, */
    /*   ghost_data, */
    /*   d4est_ops, */
    /*   d4est_geom, */
    /*   d4est_quad, */
    /*   d4est_factors, */
    /*   &krylov_petsc_params, */
    /*   (mg_data->num_of_levels == 1) ? NULL : pc */
    /* ); */

    // Compute and save mesh data to a VTK file
    
    // Compute analytical field values
    double* u_analytic = P4EST_ALLOC(double, prob_vecs.local_nodes);
    d4est_mesh_init_field(
      p4est,
      u_analytic,
      poisson_lorentzian_analytic_solution,
      d4est_ops, // unnecessary?
      d4est_geom, // unnecessary?
      d4est_factors,
      INIT_FIELD_ON_LOBATTO,
      NULL
    );

    // Compute errors between numerical and analytical field values
    double* error = P4EST_ALLOC(double, prob_vecs.local_nodes);
    d4est_linalg_vec_fabsdiff(prob_vecs.u, u_analytic, error, prob_vecs.local_nodes);

    double* error_l2 = P4EST_ALLOC(double, p4est->local_num_quadrants);
    d4est_mesh_compute_l2_norm_sqr
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       error,
       prob_vecs.local_nodes,
       NULL,
       error_l2
      );
    
    
    // Save to VTK file
    d4est_vtk_save(
      p4est,
      d4est_ops,
      input_file,
      "d4est_vtk",
      (const char * []){"u","u_analytic","error", NULL},
      (double* []){prob_vecs.u, u_analytic, error},
      (const char * []){"error_l2",NULL},
      (double* []){error_l2},
      level
    );
    
    // Compute and save norms
    d4est_norms_save(
      p4est,
      d4est_factors,
      (const char * []){ "u", NULL },
      (double * []){ prob_vecs.u },
      (double * []){ u_analytic }, // Using precomputed analytic field values
      (d4est_xyz_fcn_t[]){ NULL },
      (void * []) { NULL },
      (const char * []){"L_2", "L_infty", NULL},
      (d4est_norm_fcn_t[]){ &d4est_norms_fcn_L2, &d4est_norms_fcn_Linfty },
      (void * []){ &L2_norm_ctx, NULL },
      NULL
    );

    P4EST_FREE(error_l2);    
    P4EST_FREE(error);
    P4EST_FREE(u_analytic);


    // Perform the next AMR step
  
    if (level != d4est_amr->num_of_amr_steps){

      if (p4est->mpirank == 0)
        zlog_info(c_default, "Performing AMR refinement level %d of %d...", level + 1, d4est_amr->num_of_amr_steps);

      d4est_amr_step(
        p4est,
        ghost,
        ghost_data,
        d4est_ops,
        d4est_amr,
        &prob_vecs.u,
        NULL,
        NULL
      );
      
      if (p4est->mpirank == 0)
        zlog_info(c_default, "AMR refinement level %d of %d complete.", level + 1, d4est_amr->num_of_amr_steps);
    }


    prob_vecs.local_nodes = d4est_mesh_update(
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
    
    
    d4est_poisson_build_rhs_with_strong_bc(
      p4est,
      *ghost,
      *ghost_data,
      d4est_ops,
      d4est_geom,
      d4est_quad,
      d4est_factors,
      &prob_vecs,
      flux_data_for_rhs,
      prob_vecs.rhs,
      poisson_lorentzian_rhs_fcn,
      INIT_FIELD_ON_LOBATTO,
      &ctx
    );


      /* krylov_pc_multigrid_destroy(pc); */
      /* multigrid_data_destroy(mg_data); */

    
  }

  if (p4est->mpirank == 0)
    zlog_info(c_default, "Finishing up. Starting garbage collection...");
    
  d4est_amr_destroy(d4est_amr);
  d4est_poisson_flux_destroy(flux_data_for_lhs);
  d4est_poisson_flux_destroy(flux_data_for_rhs);
  P4EST_FREE(prob_vecs.u);
  P4EST_FREE(prob_vecs.Au);
  P4EST_FREE(prob_vecs.rhs);
}
