#define _GNU_SOURCE
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
#include <d4est_laplacian.h>
#include <d4est_laplacian_flux_sipg.h>
#include <d4est_solver_newton_petsc.h>
#include <d4est_solver_krylov_petsc.h>
#include <d4est_krylov_pc_multigrid.h>
#include <d4est_solver_multigrid_logger_residual.h>
#include <d4est_solver_multigrid_element_data_updater.h>
#include <d4est_solver_multigrid.h>
#include <d4est_solver_schwarz.h>
#include <d4est_solver_schwarz_laplacian_ext.h>
#include <d4est_util.h>
#include <time.h>
#include <zlog.h>
#include "poisson_sinx_fcns.h"



static void
poisson_sinx_schwarz_apply_lhs
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* d4est_ghost,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data,
 int subdomain,
 double* u_restricted_field_over_subdomain,
 double* Au_restricted_field_over_subdomain,
 void* ctx
){

  d4est_solver_schwarz_laplacian_ext_apply_over_subdomain
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     d4est_ghost,
     schwarz_data,
     schwarz_ops,
     schwarz_geometric_data,
     ctx,
     u_restricted_field_over_subdomain,
     Au_restricted_field_over_subdomain,
     subdomain
    );
}


void
problem_init
(
 p4est_t* p4est,
 d4est_ghost_t** d4est_ghost,
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
  
  dirichlet_bndry_eval_method_t eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;
  
  // Setup boundary conditions
  d4est_laplacian_dirichlet_bc_t bc_data_for_lhs;
  bc_data_for_lhs.dirichlet_fcn = zero_fcn;
  bc_data_for_lhs.eval_method = eval_method;
  
  d4est_laplacian_dirichlet_bc_t bc_data_for_rhs;
  bc_data_for_rhs.dirichlet_fcn = poisson_sinx_boundary_fcn;
  bc_data_for_rhs.eval_method = eval_method;
  
  d4est_laplacian_flux_data_t* flux_data_for_apply_lhs = d4est_laplacian_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_for_lhs);
  
  d4est_laplacian_flux_data_t* flux_data_for_build_rhs = d4est_laplacian_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_for_rhs);

  problem_ctx_t ctx;
  ctx.flux_data_for_apply_lhs = flux_data_for_apply_lhs;
  ctx.flux_data_for_build_rhs = flux_data_for_build_rhs;


  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = poisson_sinx_build_residual;
  prob_fcns.apply_lhs = poisson_sinx_apply_lhs;
  prob_fcns.user = &ctx;


  d4est_elliptic_data_t prob_vecs;
  prob_vecs.Au = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.u = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.rhs = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.local_nodes = initial_nodes;

  d4est_laplacian_flux_sipg_params_t* sipg_params = flux_data_for_apply_lhs->flux_data;
  
  
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
      (const char * []){"L_2", "L_infty", NULL},
      NULL
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
    poisson_sinx_initial_guess,
    d4est_ops,
    d4est_geom,
    d4est_factors,
    INIT_FIELD_ON_LOBATTO,
    NULL
  );

  d4est_field_type_t field_type = NODAL;

  
  /* d4est_laplacian_build_rhs_with_strong_bc( */
  /*   p4est, */
  /*   *d4est_ghost, */
  /*   d4est_ghost_data, */
  /*   d4est_ops, */
  /*   d4est_geom, */
  /*   d4est_quad, */
  /*   d4est_factors, */
  /*   &prob_vecs, */
  /*   flux_data_for_build_rhs, */
  /*   prob_vecs.rhs, */
  /*   poisson_sinx_rhs_fcn, */
  /*   INIT_FIELD_ON_LOBATTO, */
  /*   &ctx, */
  /*   0 */
  /* ); */

  for (int level = 0; level < d4est_amr->num_of_amr_steps + 1; level++) {


    d4est_ghost_data_t* d4est_ghost_data = d4est_ghost_data_init(p4est,
                                                                 *d4est_ghost,
                                                                 &field_type,
                                                                 1);


    
    d4est_laplacian_build_rhs_with_strong_bc(
      p4est,
      *d4est_ghost,
      d4est_ghost_data,
      d4est_ops,
      d4est_geom,
      d4est_quad,
      d4est_factors,
      &prob_vecs,
      flux_data_for_build_rhs,
      prob_vecs.rhs,
      poisson_sinx_rhs_fcn,
      INIT_FIELD_ON_LOBATTO,
      &ctx,
      0
    );

    
    // Setup d4est_solver_multigrid
    d4est_krylov_pc_t* pc = NULL;
    d4est_solver_multigrid_t* mg_data = d4est_solver_multigrid_data_init(
      p4est,
      d4est_ops,
      d4est_geom,
      d4est_quad,
      d4est_ghost,
      &d4est_ghost_data,
      d4est_factors,
      initial_extents,
      input_file
    );

    pc = d4est_krylov_pc_multigrid_create(mg_data, NULL);

    // Krylov PETSc solve
    
    d4est_solver_krylov_petsc_params_t d4est_solver_krylov_petsc_params;
    d4est_solver_krylov_petsc_input(p4est, input_file, "d4est_solver_krylov_petsc", &d4est_solver_krylov_petsc_params);

    prob_vecs.field_types = &field_type;
    prob_vecs.num_of_fields = 1;
      
    /* d4est_solver_krylov_petsc_solve( */
    /*   p4est, */
    /*   &prob_vecs, */
    /*   &prob_fcns, */
    /*   d4est_ghost, */
    /*   &d4est_ghost_data, */
    /*   d4est_ops, */
    /*   d4est_geom, */
    /*   d4est_quad, */
    /*   d4est_factors, */
    /*   &d4est_solver_krylov_petsc_params, */
    /*   pc,level */
    /* ); */


  d4est_solver_schwarz_apply_lhs_t*
    apply_lhs = d4est_solver_schwarz_apply_lhs_init
    (
     poisson_sinx_schwarz_apply_lhs,
     flux_data_for_apply_lhs
    );
    
    d4est_solver_schwarz_t* schwarz =
      d4est_solver_schwarz_init
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       *d4est_ghost,
       d4est_factors,
       NULL,
       apply_lhs,
       input_file,
       "d4est_solver_schwarz"
      );

    double* residual = P4EST_ALLOC(double, d4est_factors->local_sizes.local_nodes);
    int schwarz_iter = 40;
    for (int si = 0; si < schwarz_iter; si++){

    d4est_elliptic_eqns_apply_lhs
      (
       p4est,
       *d4est_ghost,
       d4est_ghost_data,
       &prob_fcns,
       &prob_vecs,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors
      );

    d4est_linalg_vec_axpyeqz(-1., prob_vecs.Au, prob_vecs.rhs, residual, prob_vecs.local_nodes);    

    printf("r2 norm = %.15f\n", d4est_linalg_vec_dot(residual, residual, prob_vecs.local_nodes));
    
    /* DEBUG_PRINT_3ARR_DBL(prob_vecs.rhs, prob_vecs.u, residual, prob_vecs.local_nodes); */
    
    d4est_solver_schwarz_iterate
      (
       p4est,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       *d4est_ghost,
       schwarz,
       prob_vecs.u,
       residual
      );



    /* char* schwarz_folder; */
    /* asprintf(&schwarz_folder,"Schwarz_%d_%d/", 0, 0); */
    /* char* full_dir = d4est_util_add_cwd(schwarz_folder); */
    /* d4est_util_make_directory(full_dir,0); */
    
    /* d4est_solver_schwarz_debug_vtk */
    /*   ( */
    /*    p4est, */
    /*    schwarz, */
    /*    (char*)input_file, */
    /*    "d4est_vtk", */
    /*    "schwarz_mglevel", */
    /*    schwarz_folder, */
    /*    0, */
    /*    NULL,//(const char * []){"residual",NULL}, */
    /*    NULL */
    /*   ); */

    /* free(schwarz_folder); */
    /* free(full_dir); */

    }


    d4est_solver_schwarz_destroy
      (
       schwarz
      );

    P4EST_FREE(residual);

    d4est_solver_schwarz_apply_lhs_destroy
      (
       apply_lhs
      );

    
    /* d4est_solver_cg_params_t fcg_params; */
    /* d4est_solver_cg_input(p4est, input_file, "d4est_solver_cg", "[D4EST_SOLVER_CG]", &fcg_params); */
    
    /* d4est_solver_cg_solve */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs, */
    /*    &prob_fcns, */
    /*    ghost, */
    /*    ghost_data, */
    /*    d4est_ops, */
    /*    d4est_geom, */
    /*    d4est_quad, */
    /*    d4est_factors, */
    /*    &fcg_params, */
    /*    pc */
    /*   ); */

    // Compute and save mesh data to a VTK file
    
    // Compute analytical field values
    double* u_analytic = P4EST_ALLOC(double, prob_vecs.local_nodes);
    d4est_mesh_init_field(
      p4est,
      u_analytic,
      poisson_sinx_analytic_solution,
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
      NULL,
      NULL,
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
      NULL,
      NULL,
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
        d4est_ops,
        d4est_amr,
        &prob_vecs.u,
        NULL,
        NULL,
        NULL
      );
      
      if (p4est->mpirank == 0)
        zlog_info(c_default, "AMR refinement level %d of %d complete.", level + 1, d4est_amr->num_of_amr_steps);
    }



    d4est_mesh_local_sizes_t local_sizes = d4est_mesh_update(
      p4est,
      d4est_ghost,
      d4est_ops,
      d4est_geom,
      d4est_quad,
      d4est_factors,
      initial_extents,
      INITIALIZE_GHOST,
      INITIALIZE_QUADRATURE_DATA,
      INITIALIZE_GEOMETRY_DATA,
      INITIALIZE_GEOMETRY_ALIASES,
      d4est_mesh_set_quadratures_after_amr,
      initial_extents
    );

    prob_vecs.local_nodes = local_sizes.local_nodes;
      
    prob_vecs.Au = P4EST_REALLOC(prob_vecs.Au, double, prob_vecs.local_nodes);
    prob_vecs.rhs = P4EST_REALLOC(prob_vecs.rhs, double, prob_vecs.local_nodes);
    
    d4est_krylov_pc_multigrid_destroy(pc);
    d4est_solver_multigrid_data_destroy(mg_data);

    if (d4est_ghost_data != NULL){
      d4est_ghost_data_destroy(d4est_ghost_data);
      d4est_ghost_data = NULL;
    } 


    
  }

  if (p4est->mpirank == 0)
    zlog_info(c_default, "Finishing up. Starting garbage collection...");
    
  d4est_amr_destroy(d4est_amr);
  d4est_laplacian_flux_destroy(flux_data_for_apply_lhs);
  d4est_laplacian_flux_destroy(flux_data_for_build_rhs);
  P4EST_FREE(prob_vecs.u);
  P4EST_FREE(prob_vecs.Au);
  P4EST_FREE(prob_vecs.rhs);
}
