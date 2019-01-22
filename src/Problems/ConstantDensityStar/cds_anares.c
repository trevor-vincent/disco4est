#include <sc_reduce.h>
#include <pXest.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <problem.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_estimator_bi_new.h>
#include <d4est_solver_cg.h>
#include <d4est_amr.h>
#include <d4est_amr_smooth_pred.h>
#include <d4est_geometry.h>
#include <d4est_geometry_brick.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_vtk.h>
#include <d4est_norms.h>
#include <d4est_mesh.h>
#include <ini.h>
#include <d4est_element_data.h>
#include <d4est_estimator_stats.h>
#include <d4est_laplacian_with_opt.h>
#include <d4est_laplacian_with_opt_flux_sipg.h>
#include <d4est_solver_newton.h>
#include <d4est_solver_multigrid.h>
#include <d4est_krylov_pc_multigrid.h>
#include <d4est_solver_multigrid_logger_residual.h>
#include <d4est_solver_multigrid_element_data_updater.h>
#include <d4est_solver_multigrid_matrix_operator.h>
#include <d4est_solver_krylov_petsc.h>
#include <d4est_solver_newton_petsc.h>
#include <d4est_util.h>
#include <d4est_h5.h>
#include <d4est_checkpoint.h>
#include <p4est_vtk_ext.h>
#include <time.h>
#include "constant_density_star_fcns.h"


typedef struct {
  
  int do_not_solve;
  int amr_level_for_uniform_p;
  
} constant_density_star_init_params_t;

static
int constant_density_star_init_params_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  constant_density_star_init_params_t* pconfig = (constant_density_star_init_params_t*)user;
  if (d4est_util_match_couple(section,"problem",name,"do_not_solve")) {
    D4EST_ASSERT(pconfig->do_not_solve == -1);
    pconfig->do_not_solve = atoi(value);
  }
  else if (d4est_util_match_couple(section,"amr",name,"amr_level_for_uniform_p")) {
    D4EST_ASSERT(pconfig->amr_level_for_uniform_p == -1);
    pconfig->amr_level_for_uniform_p = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


static
constant_density_star_init_params_t
constant_density_star_init_params_input
(
 const char* input_file
)
{
  constant_density_star_init_params_t input;
  input.do_not_solve = -1;
  input.amr_level_for_uniform_p = -1;

  if (ini_parse(input_file, constant_density_star_init_params_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("problem", input.do_not_solve, -1);
  D4EST_CHECK_INPUT("amr", input.amr_level_for_uniform_p, -1);
  
  return input;
}

static
int
amr_mark_element
(
 p4est_t* p4est,
 double eta2,
 d4est_estimator_stats_t* stats,
 d4est_element_data_t* elem_data,
 d4est_amr_smooth_pred_params_t* params,
 void* user
)
{
  problem_ctx_t* ctx = user;

  double eta2_percentile = stats->estimator_at_percentile;
    /* = d4est_estimator_stats_get_percentile(stats,params->percentile); */
  
  return ((eta2 >= eta2_percentile) || fabs(eta2 - eta2_percentile) < eta2*1e-4);
}

static
gamma_params_t
amr_set_element_gamma
(
 p4est_t* p4est,
 d4est_estimator_stats_t* stats,
 d4est_element_data_t* elem_data,
 d4est_amr_smooth_pred_params_t* params,
 void* user
)
{
  problem_ctx_t* ctx = user;
  gamma_params_t gamma_hpn;
  gamma_hpn.gamma_h = params->gamma_h;
  gamma_hpn.gamma_p = params->gamma_p;
  gamma_hpn.gamma_n = params->gamma_n;
  return gamma_hpn;
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
  int initial_nodes = initial_extents->initial_nodes;
  constant_density_star_init_params_t init_params = constant_density_star_init_params_input(input_file);

  constant_density_star_params_t constant_density_star_params = constant_density_star_input(input_file);
  
  d4est_laplacian_with_opt_dirichlet_bc_t bc_data_for_jac;
  bc_data_for_jac.dirichlet_fcn = zero_fcn;
  bc_data_for_jac.user = &constant_density_star_params;
  bc_data_for_jac.eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;
  
  d4est_laplacian_with_opt_dirichlet_bc_t bc_data_for_res;
  bc_data_for_res.dirichlet_fcn = constant_density_star_boundary_fcn;
  bc_data_for_res.user = &constant_density_star_params;
  bc_data_for_res.eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;
                                          
  d4est_laplacian_with_opt_flux_data_t* flux_data_for_jac =
    d4est_laplacian_with_opt_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_for_jac);
  d4est_laplacian_with_opt_flux_data_t* flux_data_for_res = d4est_laplacian_with_opt_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_for_res);

  problem_ctx_t ctx;
  ctx.constant_density_star_params = &constant_density_star_params;
  ctx.flux_data_for_jac = flux_data_for_jac;
  ctx.flux_data_for_res = flux_data_for_res;
                           
  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = constant_density_star_build_residual;
  prob_fcns.apply_lhs = constant_density_star_apply_jac;
  prob_fcns.user = &ctx;
  
 
  d4est_elliptic_data_t prob_vecs;
  prob_vecs.Au = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.u = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.local_nodes = initial_nodes;

  double* error = P4EST_ALLOC(double, prob_vecs.local_nodes);
  double* u_analytic = P4EST_ALLOC(double, prob_vecs.local_nodes);
  
  d4est_laplacian_with_opt_flux_sipg_params_t* sipg_params = flux_data_for_jac->flux_data;
  
  d4est_estimator_bi_new_penalty_data_t penalty_data;
  penalty_data.u_penalty_fcn = houston_u_prefactor_maxp_minh;
  penalty_data.size_params = NULL;
  penalty_data.u_dirichlet_penalty_fcn = houston_u_dirichlet_prefactor_maxp_minh;
  penalty_data.gradu_penalty_fcn = houston_gradu_prefactor_maxp_minh;
  penalty_data.penalty_prefactor = sipg_params->sipg_penalty_prefactor;
  penalty_data.user = &constant_density_star_params;
  
  d4est_amr_smooth_pred_marker_t amr_marker;
  amr_marker.user = (void*)&ctx;
  amr_marker.mark_element_fcn = amr_mark_element;
  amr_marker.set_element_gamma_fcn = amr_set_element_gamma;

  d4est_amr_t* d4est_amr =
    d4est_amr_init
    (
     p4est,
     input_file,
     &amr_marker
    );

  int initial_level = 0;
  if (initial_extents->load_from_checkpoint == 0 || initial_extents->checkpoint_prefix == NULL){
    d4est_mesh_init_field
      (
       p4est,
       prob_vecs.u,
       constant_density_star_initial_guess,
       d4est_ops,
       d4est_geom,
       d4est_factors,
       INIT_FIELD_ON_LOBATTO,
       NULL
      );
  }
  else {


    d4est_checkpoint_read_dataset
      (
       p4est,
       initial_extents->checkpoint_prefix,
       "u",
       H5T_NATIVE_DOUBLE,
       prob_vecs.u,
       initial_extents->checkpoint_number
      );


    double sum = d4est_util_sum_array_dbl(prob_vecs.u, prob_vecs.local_nodes);

    d4est_checkpoint_check_dataset(p4est,
                           initial_extents->checkpoint_prefix,
                           "u",
                           H5T_NATIVE_DOUBLE,
                           (void*)&sum,
                           initial_extents->checkpoint_number
                          );

    initial_level = initial_extents->checkpoint_number + 1;
/* d4est_h5_read_dataset(p4est->mpirank,initial_extents->checkpoint_prefix,"u",H5T_NATIVE_DOUBLE, prob_vecs.u); */
  }
  /* keep track of mg_levels on the previous amr step 
   * to calculate the mg_levels on the next amr step */
  int num_of_mg_levels_last_step = -1;
  
  // Norm function contexts  
  d4est_norms_fcn_L2_ctx_t L2_norm_ctx;
  L2_norm_ctx.p4est = p4est;
  L2_norm_ctx.d4est_ops = d4est_ops;
  L2_norm_ctx.d4est_geom = d4est_geom;
  L2_norm_ctx.d4est_quad = d4est_quad;
  L2_norm_ctx.d4est_factors = d4est_factors;
  
  d4est_norms_fcn_energy_ctx_t energy_norm_ctx;
  energy_norm_ctx.p4est = p4est;
  energy_norm_ctx.d4est_ops = d4est_ops;
  energy_norm_ctx.d4est_geom = d4est_geom;
  energy_norm_ctx.d4est_quad = d4est_quad;
  energy_norm_ctx.d4est_factors = d4est_factors;
  /* energy_norm_ctx.fit = NULL; */
  // These are updated later
  energy_norm_ctx.which_field = 0;
  energy_norm_ctx.energy_norm_data = NULL;
  energy_norm_ctx.energy_estimator_sq_local = -1.;

  /* if (p4est->mpirank == 0) */
  /*   d4est_norms_write_headers( */
  /*                             (const char * []){"u", NULL}, */
  /*                             (const char * []){"L_2", "L_infty", "energy_norm", "energy_estimator", NULL} */
  /*   ); */

  zlog_category_t* zlog_points = zlog_get_category("d4est_points");
  d4est_util_copy_1st_to_2nd(prob_vecs.u, u_analytic, prob_vecs.local_nodes);
  
  d4est_norms_linear_fit_t* l2_linear_fit = d4est_norms_linear_fit_init();
  
  for (int level = initial_level; level < d4est_amr->num_of_amr_steps + 1; ++level){

    d4est_field_type_t field_type = NODAL;
    d4est_ghost_data_t* d4est_ghost_data = d4est_ghost_data_init(p4est,
                                                                 *d4est_ghost,
                                                                 &field_type,
                                                                 1);

    /* d4est_mesh_size_parameters_t size_params = d4est_mesh_get_size_parameters(d4est_factors_compactified); */
    d4est_ip_energy_norm_data_t ip_norm_data;
    penalty_data.size_params = NULL;
    ip_norm_data.size_params = NULL;
    sipg_params->size_params = NULL;

    /* if (init_params.use_compactified_size_params){ */
    /*   penalty_data.size_params = &size_params; */
    /*   ip_norm_data.size_params = &size_params; */
    /*   sipg_params->size_params = &size_params; */
    /* } */

    double* residual_pointwise_quad
      = P4EST_ALLOC(double, d4est_factors->local_sizes.local_nodes_quad);
    
    constant_density_star_pointwise_residual
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       prob_vecs.u,
       residual_pointwise_quad,
       0,
       &ctx
      );
    
    double* estimator =
      d4est_estimator_bi_new_compute
      (
       p4est,
       &prob_vecs,
       residual_pointwise_quad,
       penalty_data,
       zero_fcn,
       NULL,
       *d4est_ghost,
       d4est_ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_factors,
       d4est_geom,
       d4est_factors,
       d4est_quad,
       0,
       NULL,
       NULL,
       1,
       NULL
      );

    P4EST_FREE(residual_pointwise_quad);


    d4est_mesh_init_field(
                          p4est,
                          u_analytic,
                          constant_density_star_analytic_solution,
                          d4est_ops,
                          d4est_geom,
                          d4est_factors,
                          INIT_FIELD_ON_LOBATTO,
                          &ctx
    );
    d4est_amr_smooth_pred_params_t* sp_params = d4est_amr_smooth_pred_params_input
                                                (
                                                 input_file
                                                );

    d4est_estimator_stats_t* stats = P4EST_ALLOC(d4est_estimator_stats_t,1);
    d4est_estimator_stats_compute(p4est, estimator, stats, sp_params->percentile, 1, 0);
    d4est_linalg_vec_fabsdiff(prob_vecs.u, u_analytic, error, prob_vecs.local_nodes);
    double* error_l2 = P4EST_ALLOC(double, p4est->local_num_quadrants);
    P4EST_FREE(sp_params);
    
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

    /* if(init_params.use_error_l2_as_estimator){ */
      /* d4est_util_copy_1st_to_2nd(error_l2, estimator, p4est->local_num_quadrants); */
    /* }     */
    
    d4est_amr_smooth_pred_data_t* smooth_pred_data = (d4est_amr_smooth_pred_data_t*) (d4est_amr->scheme->amr_scheme_data);

    d4est_vtk_save
      (
       p4est,
       d4est_ops,
       input_file,
       "d4est_vtk",
       (const char * []){"u","u_analytic","error", NULL},
       (double* []){prob_vecs.u, u_analytic, error},
       (const char * []){"estimator","error_l2",NULL},
       (double* []){estimator,error_l2},
       NULL,
       NULL,
       level
      );

    double* u_vertex = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN));

    d4est_element_data_store_nodal_vec_in_vertex_array
      (
       p4est,
       prob_vecs.u,
       u_vertex
      );

    char* u_corner_file = P4EST_ALLOC(char, 100);
    sprintf(u_corner_file, "%s_%d", "u_corner", level);
    p4est_vtk_ext_write_all
      (p4est,
       NULL,
       0.99,
       1,
       1,
       1,
       1,
       1,
       0,
       u_corner_file,
       "u",
       u_vertex
      );

    P4EST_FREE(u_corner_file);
    P4EST_FREE(u_vertex);

    

    ip_norm_data.u_penalty_fcn = sipg_params->sipg_penalty_fcn;
    ip_norm_data.penalty_prefactor = sipg_params->sipg_penalty_prefactor;


    double total_est = stats->estimator_total;    
    energy_norm_ctx.energy_norm_data = &ip_norm_data;
    energy_norm_ctx.energy_estimator_sq_local = total_est;
    energy_norm_ctx.ghost = *d4est_ghost;
    energy_norm_ctx.ghost_data = d4est_ghost_data;

    d4est_norms_save
      (
       p4est,
       d4est_factors,
       (const char * []){ "u", NULL },
       (double * []){ prob_vecs.u },
       (double * []){ u_analytic },
       (d4est_xyz_fcn_t []){ NULL },
       (void * []){ NULL },
       (const char * []){"L_2", "L_infty", "energy_norm", "energy_estimator", NULL},
       (d4est_norm_fcn_t[]){ &d4est_norms_fcn_L2, &d4est_norms_fcn_Linfty, &d4est_norms_fcn_energy, &d4est_norms_fcn_energy_estimator },
       (void * []){ &L2_norm_ctx, NULL, &energy_norm_ctx, &energy_norm_ctx },
       (d4est_norms_linear_fit_t * []){ l2_linear_fit, NULL, NULL, NULL },
       NULL,
       NULL
      );
    
    if (level != d4est_amr->num_of_amr_steps && level != 0){

      d4est_amr_step
        (
         p4est,
         d4est_ops,
         d4est_amr,
         &prob_vecs.u,
         estimator,
         stats,
         input_file
        );
      
    }

    P4EST_FREE(stats);
    


      d4est_mesh_local_sizes_t local_sizes= d4est_mesh_update
                            (
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

    if (d4est_ghost_data != NULL){
      d4est_ghost_data_destroy(d4est_ghost_data);
      d4est_ghost_data = NULL;
    } 
    

    d4est_ghost_data = d4est_ghost_data_init(p4est,
                                             *d4est_ghost,
                                             &field_type,
                                             1);

    
    
    prob_vecs.Au = P4EST_REALLOC(prob_vecs.Au, double, prob_vecs.local_nodes);
    u_analytic = P4EST_REALLOC(u_analytic, double, prob_vecs.local_nodes);
    error = P4EST_REALLOC(error, double, prob_vecs.local_nodes);
    d4est_util_copy_1st_to_2nd(prob_vecs.u, u_analytic, prob_vecs.local_nodes);
      
      d4est_solver_multigrid_t* mg_data = d4est_solver_multigrid_data_init(p4est,
                                                      d4est_ops,
                                                      d4est_geom,
                                                      d4est_quad,
                                                      d4est_ghost,
                                                      &d4est_ghost_data,
                                                      d4est_factors,
                                                      initial_extents,
                                                      input_file
                                                     );



      /* This new multigrid level guesser is more accurate than the last */
      /* if(init_params.use_new_multigrid_level_guesser == 1){ */
      /*   if (level == 0){ */
      /*     mg_data->num_of_levels = d4est_solver_multigrid_get_h_coarsen_levels_initial */
      /*                              ( */
      /*                               p4est, */
      /*                               initial_extents */
      /*                              ); */
      /*   } */
      /*   else { */
      /*     int has_there_been_h_refinements_local = (d4est_amr->has_there_been_h_refinements > 0);                     */
      /*     int has_there_been_h_refinements_global = -1; */
      /*     sc_allreduce */
      /*       ( */
      /*        &has_there_been_h_refinements_local, */
      /*        &has_there_been_h_refinements_global, */
      /*        1, */
      /*        sc_MPI_INT, */
      /*        sc_MPI_MAX, */
      /*        p4est->mpicomm */
      /*     ); */

      /*     if (initial_extents->load_from_checkpoint == 1 */
      /*         && */
      /*         level == initial_extents->checkpoint_number) */
      /*       { */
      /*         d4est_checkpoint_read_dataset */
      /*           ( */
      /*            p4est, */
      /*            initial_extents->checkpoint_prefix, */
      /*            "multigrid_h_levels", */
      /*            H5T_NATIVE_DOUBLE, */
      /*            &num_of_mg_levels_last_step, */
      /*            initial_extents->checkpoint_number */
      /*           );             */
      /*       } */

      /*     mg_data->num_of_levels = d4est_solver_multigrid_get_h_coarsen_levels_post_initial */
      /*                              ( */
      /*                               p4est, */
      /*                               num_of_mg_levels_last_step, */
      /*                               has_there_been_h_refinements_global */
      /*                              ); */

      /*   } */
      /*   num_of_mg_levels_last_step = mg_data->num_of_levels; */
      /* }       */
      
      d4est_solver_multigrid_user_callbacks_t* user_callbacks = d4est_solver_multigrid_matrix_operator_init(p4est, mg_data->num_of_levels);

      d4est_solver_multigrid_set_user_callbacks(
                            mg_data,
                            user_callbacks
                           );     

      d4est_krylov_pc_t* pc = d4est_krylov_pc_multigrid_create(mg_data, constant_density_star_pc_setup_fcn);
      ctx.use_matrix_operator = 1;
      ctx.mg_data = mg_data;
      zlog_category_t *c_cds = zlog_get_category("cds");
    
    if (!init_params.do_not_solve){

      d4est_solver_newton_petsc_params_t newton_params;
      d4est_solver_newton_petsc_input(p4est, input_file, &newton_params);

      d4est_solver_krylov_petsc_params_t krylov_params;

      if (mg_data->num_of_levels <= 1){
        d4est_solver_krylov_petsc_input(p4est, input_file, "d4est_solver_krylov_petsc_no_mg", &krylov_params);
      }
      else {
        d4est_solver_krylov_petsc_input(p4est, input_file, "d4est_solver_krylov_petsc", &krylov_params);
      }

      
      prob_vecs.field_types = &field_type;
      prob_vecs.num_of_fields = 1;

      if (initial_extents->load_newton_checkpoint == 1){
        if (level == initial_level){
          if (p4est->mpirank == 0)
            zlog_info(c_cds,"Loading from newton checkpoint if level %d == initial_level %d",
                 level, initial_level);
        }
        else{
          if (p4est->mpirank == 0)
            zlog_info(c_cds,"We will not load from newton checkpoint because level %d != initial_level %d"
                 ,level, initial_level);
        }
      }
      else{
        if (p4est->mpirank == 0)
          zlog_info(c_cds,"We will not load from newton checkpoint");
      }
      if (level == initial_level && initial_extents->load_newton_checkpoint == 1){

        d4est_checkpoint_read_dataset
          (
           p4est,
           initial_extents->newton_checkpoint_prefix,
           "u",
           H5T_NATIVE_DOUBLE,
           prob_vecs.u,
           level
          );


        double sum = d4est_util_sum_array_dbl(prob_vecs.u, prob_vecs.local_nodes);

        d4est_checkpoint_check_dataset(p4est,
                                       initial_extents->newton_checkpoint_prefix,
                                       "u",
                                       H5T_NATIVE_DOUBLE,
                                       (void*)&sum,
                                       level
                                      );



        zlog_info(c_cds, "Loading u from a newton checkpoint %s at level %d", initial_extents->newton_checkpoint_prefix, level);
      }
      
      d4est_solver_newton_petsc_solve
        (
         p4est,
         &prob_vecs,
         &prob_fcns,
         d4est_ghost,
         &d4est_ghost_data,
         d4est_ops,
         d4est_geom,
         d4est_quad,
         d4est_factors,
         &krylov_params,
         &newton_params,
         (mg_data->num_of_levels <= 1) ? NULL : pc,
         level
        );
    }
    
    if (level != d4est_amr->num_of_amr_steps && level != 0){
      d4est_checkpoint_save
        (
         level,
         "checkpoint",
         p4est,
         d4est_amr,
         d4est_factors,
         (const char * []){"u", "predictor", "multigrid_h_levels", NULL},
         (hid_t []){H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_INT},
         (int []){prob_vecs.local_nodes, p4est->local_num_quadrants, 1},
         (void* []){prob_vecs.u, smooth_pred_data->predictor, &mg_data->num_of_levels}
        );
    }

    d4est_krylov_pc_multigrid_destroy(pc);
    d4est_solver_multigrid_data_destroy(mg_data);
    d4est_solver_multigrid_matrix_operator_destroy(user_callbacks);
   
    P4EST_FREE(error_l2);
    P4EST_FREE(estimator);

    if (d4est_ghost_data != NULL){
      d4est_ghost_data_destroy(d4est_ghost_data);
      d4est_ghost_data = NULL;
    } 
  }
  
  printf("[D4EST_INFO]: Starting garbage collection...\n");
  /* d4est_mesh_data_destroy(d4est_factors_compactified); */
  /* d4est_geometry_destroy(d4est_geom_compactified); */
  d4est_amr_destroy(d4est_amr);
  d4est_norms_linear_fit_destroy(l2_linear_fit);
  d4est_laplacian_with_opt_flux_destroy(flux_data_for_jac);
  d4est_laplacian_with_opt_flux_destroy(flux_data_for_res);
  P4EST_FREE(error);
  
  P4EST_FREE(u_analytic);
  P4EST_FREE(prob_vecs.u);
  P4EST_FREE(prob_vecs.Au);
}
