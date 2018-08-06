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
#include <d4est_norms.h>
#include <d4est_mesh.h>
#include <ini.h>
#include <d4est_element_data.h>
#include <d4est_poisson.h>
#include <d4est_poisson_flux_sipg.h>
#include <d4est_solver_newton_petsc.h>
#include <d4est_solver_krylov_petsc.h>
#include <d4est_util.h>
#include <time.h>
#include <d4est_solver_multigrid.h>
#include <d4est_solver_multigrid_logger_residual.h>
#include <d4est_krylov_pc_multigrid.h>
#include <d4est_solver_multigrid_element_data_updater.h>
#include "poisson_lorentzian_fcns.h"


static
double solve_for_c
(
 double c,
 void* user
)
{
  double* Rs = user;
  double R1 = Rs[0];
  double R2 = Rs[1];
  double Rc = Rs[2];
  double m = (2 - 1)/((1/R2) - (1/R1));
  double n = (1*R1 - 2*R2)/(R1 - R2);
  double R = m/(c - n);
  double pp = 2 - c;
  double q = R/sqrt(1 + 2*pp);
  double x = q;
  return x - Rc;  
}



static
double solve_for_c_outer
(
 double c,
 void* user
)
{
  double* Rs = user;
  double R1 = Rs[0];
  double R2 = Rs[1];
  double Rc = Rs[2];
  double m = (2 - 1)/((1/R2) - (1/R1));
  double n = (1*R1 - 2*R2)/(R1 - R2);
  double R = m/(c - n);
  double pp = 2 - c;
  double q = R;
  double x = q;
  return x - Rc;  
}


double
get_inverted_outer_wedge_point(double R1, double R2, double Rc, int compactified){
  D4EST_ASSERT(Rc >= R1 && Rc <= R2);
  if (compactified){
    double c;
    if (Rc == R2){
      c = 2;
    }
    else {
      double Rs [] = {R1,R2,Rc};
      int success = d4est_util_bisection(solve_for_c_outer, 1, 2, DBL_EPSILON, 100000, &c, &Rs[0]);
      D4EST_ASSERT(!success);
    }
    return c - 1;
  }
  else{
    D4EST_ABORT("get_inverted_outer_wedge_point not accepted yet");
  }
}


double
get_inverted_inner_wedge_point(double R1, double R2, double Rc, int compactified){
  D4EST_ASSERT(Rc >= R1 && Rc <= R2);
  if (compactified){
    double c;
    if (Rc == R2){
      c = 2;
    }
    else {
      double Rs [] = {R1,R2,Rc};
      int success = d4est_util_bisection(solve_for_c, 1, 2, DBL_EPSILON, 100000, &c, &Rs[0]);
      D4EST_ASSERT(!success);
    }
    return c - 1;
  }
  else{
    return ((2*pow(R1,2) - 3*R1*R2 + pow(R2,2) - pow(Rc,2) + sqrt(pow(Rc,2)*(pow(R1,2) - 4*R1*R2 + 3*pow(R2,2) + pow(Rc,2))))/pow(R1 - R2,2)) - 1;
  }
}

double
get_inverted_box_point(double R0, double x){
  double a = R0/sqrt(3);
  D4EST_ASSERT(x <= a && x >= -a);
  return (x + a)/(2*a);
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




typedef struct {
  
  int use_dirichlet;
  
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

  if (ini_parse(input_file, poisson_lorentzian_init_params_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("problem", input.use_dirichlet, -1);
  
  return input;
}

int
problem_set_mortar_degree
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  return elem_data->deg;
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
  
  d4est_poisson_flux_data_t* flux_data_for_lhs = NULL; //d4est_poisson_flux_new(p4est, input_file, BC_DIRICHLET, 
  d4est_poisson_flux_data_t* flux_data_for_rhs = NULL; //d4est_poisson_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_dirichlet_for_rhs);

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
  
  double* error = NULL;
  double* u_analytic = NULL;
  
  d4est_elliptic_data_t prob_vecs;
  prob_vecs.Au = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.u = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.rhs = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.local_nodes = initial_nodes;

  d4est_poisson_flux_sipg_params_t* sipg_params = flux_data_for_lhs->flux_data;
  
  d4est_estimator_bi_penalty_data_t penalty_data;
  penalty_data.u_penalty_fcn = houston_u_prefactor_maxp_minh;
  penalty_data.size_params = NULL;
  penalty_data.u_dirichlet_penalty_fcn = houston_u_dirichlet_prefactor_maxp_minh;
  penalty_data.gradu_penalty_fcn = houston_gradu_prefactor_maxp_minh;
  penalty_data.penalty_prefactor = sipg_params->sipg_penalty_prefactor;
  
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

  d4est_amr_t* d4est_amr_uniform = d4est_amr_init_uniform_h(p4est, d4est_amr->num_of_amr_steps);

  d4est_mesh_init_field
    (
     p4est,
     prob_vecs.u,
     poisson_lorentzian_initial_guess,
     d4est_ops,
     d4est_geom,
     d4est_factors,
     INIT_FIELD_ON_LOBATTO,
     NULL
    );
    

  d4est_field_type_t field_type = VOLUME_NODAL;
  d4est_ghost_data_t* d4est_ghost_data = d4est_ghost_data_init(p4est,
                                                               *d4est_ghost,
                                                               &field_type,
                                                               1);
  
  d4est_poisson_build_rhs_with_strong_bc
    (
     p4est,
     *d4est_ghost,
     d4est_ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &prob_vecs,
     flux_data_for_rhs,
     prob_vecs.rhs,
     poisson_lorentzian_rhs_fcn,
     INIT_FIELD_ON_LOBATTO,
     &ctx,
     0
    );


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
    energy_norm_ctx.ghost = *d4est_ghost;
    energy_norm_ctx.ghost_data = d4est_ghost_data;
    energy_norm_ctx.energy_norm_data = NULL;
    energy_norm_ctx.energy_estimator_sq_local = -1.;
    energy_norm_ctx.which_field = 0;

    if (p4est->mpirank == 0)
      d4est_norms_write_headers(
        (const char * []){"u", NULL},
        (const char * []){"L_2", "L_infty", "energy_norm", "energy_estimator", NULL}
      );

    zlog_category_t *c_geom = zlog_get_category("d4est_geometry_compactified");
  d4est_geometry_t* d4est_geom_compactified = d4est_geometry_new(p4est->mpirank, input_file,"compactified_geometry", c_geom);
  d4est_mesh_data_t* d4est_factors_compactified = d4est_mesh_data_init(p4est);



  double point [4][30];
  double point_diff [4][30];
  double point_spec_diff [4][30];
  double point_err [4];
  double point_dof [30];
  
  point[0][0] = 0;
  point_diff[0][0] = 0;
  point[1][0] = 0;
  point_diff[1][0] = 0;
  point[2][0] = 0;
  point_diff[2][0] = 0;
  point[3][0] = 0;
  point_diff[3][0] = 0;
  point_dof[0] = 0;
  point_spec_diff[0][0] = 0;
  point_spec_diff[1][0] = 0;
  point_spec_diff[2][0] = 0;
  point_spec_diff[3][0] = 0;

  int iterations = 1;
    
  for (int level = 0; level < d4est_amr->num_of_amr_steps + 1; ++level){



    d4est_mesh_data_realloc
      (
       p4est,
       *d4est_ghost,
       d4est_factors_compactified,
       d4est_factors->local_sizes
      );
    
    d4est_mesh_data_compute
      (
       p4est,
       *d4est_ghost,
       /* *ghost_data, */
       d4est_ops,
       d4est_geom_compactified,
       d4est_quad,
       d4est_factors_compactified,
       initial_extents->face_h_type,
       initial_extents->volume_h_type
      );
    
    double* estimator = d4est_estimator_bi_compute
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       penalty_data,
       poisson_lorentzian_boundary_fcn,
       *d4est_ghost,
       d4est_ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_factors,
       d4est_geom_compactified,
       d4est_factors_compactified,
       d4est_quad,
       0
      );

   d4est_amr_smooth_pred_params_t* sp_params = d4est_amr_smooth_pred_params_input
                                                (
                                                 input_file
                                                );
    
    d4est_estimator_stats_t* stats = P4EST_ALLOC(d4est_estimator_stats_t,1);
    d4est_estimator_stats_compute(p4est, estimator, stats, sp_params->percentile, 1, 0);
    d4est_estimator_stats_print(stats);

   double* error_l2 = P4EST_ALLOC(double, p4est->local_num_quadrants);

    
    
    double* u_analytic = P4EST_ALLOC(double, prob_vecs.local_nodes);
    double* error = P4EST_ALLOC(double, prob_vecs.local_nodes);

    d4est_mesh_init_field
      (
       p4est,
       u_analytic,
       poisson_lorentzian_analytic_solution,
       d4est_ops,
       d4est_geom,
       d4est_factors,
       INIT_FIELD_ON_LOBATTO,
       &ctx
      );

    for (int i = 0; i < prob_vecs.local_nodes; i++){
      error[i] = fabs(prob_vecs.u[i] - u_analytic[i]);
    }

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

    
    d4est_vtk_save
      (
       p4est,
       d4est_ops,
       input_file,
       "d4est_vtk",
       (const char * []){"u","u_analytic","error",  NULL},
       (double* []){prob_vecs.u, u_analytic, error},
       (const char * []){"estimator", "error_l2", NULL},
       (double* []){estimator, error_l2, NULL},
       level
      );

    P4EST_FREE(u_analytic);
    P4EST_FREE(error);
    P4EST_FREE(error_l2);


    d4est_ip_energy_norm_data_t ip_norm_data;
    ip_norm_data.u_penalty_fcn = sipg_params->sipg_penalty_fcn;
    ip_norm_data.penalty_prefactor = sipg_params->sipg_penalty_prefactor;
    ip_norm_data.size_params = NULL;

    energy_norm_ctx.energy_norm_data = &ip_norm_data;
    energy_norm_ctx.energy_estimator_sq_local = stats->estimator_total;
    energy_norm_ctx.ghost = *d4est_ghost;
    energy_norm_ctx.ghost_data = d4est_ghost_data;

    d4est_norms_save(
      p4est,
      d4est_factors,
      (const char * []){ "u", NULL },
      (double * []){ prob_vecs.u },
      (double * []){ NULL },
      (d4est_xyz_fcn_t []){ poisson_lorentzian_analytic_solution },
      (void * []){ &ctx },
      (const char * []){"L_2", "L_infty", "energy_norm", "energy_estimator", NULL},
      (d4est_norm_fcn_t[]){ &d4est_norms_fcn_L2, &d4est_norms_fcn_Linfty, &d4est_norms_fcn_energy, &d4est_norms_fcn_energy_estimator },
      (void * []){ &L2_norm_ctx, NULL, &energy_norm_ctx, &energy_norm_ctx },
      NULL
    );
        
    P4EST_FREE(stats);
    
    if (level != d4est_amr->num_of_amr_steps){

      if (p4est->mpirank == 0)
        printf("[D4EST_INFO]: AMR REFINEMENT LEVEL %d\n", level+1);

      d4est_amr_step
        (
         p4est,
         d4est_ops,
         (level > 1) ? d4est_amr : d4est_amr_uniform,
         &prob_vecs.u,
         estimator,
         stats
        );
      
    }



      d4est_mesh_local_sizes_t local_sizes = d4est_mesh_update
                  (
                   p4est,
                   d4est_ghost,
                   /* *ghost_data, */
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
    prob_vecs.rhs = P4EST_REALLOC(prob_vecs.rhs, double, prob_vecs.local_nodes);
    
    
    d4est_poisson_build_rhs_with_strong_bc
      (
       p4est,
       *d4est_ghost,
       d4est_ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       &prob_vecs,
       flux_data_for_rhs,
       prob_vecs.rhs,
       poisson_lorentzian_rhs_fcn,
       INIT_FIELD_ON_LOBATTO,
       &ctx,
       0
      );

   /* int min_level, max_level; */

    /* d4est_solver_multigrid_get_level_range(p4est, &min_level, &max_level); */
    /* printf("[min_level, max_level] = [%d,%d]\n", min_level, max_level); */

    /* need to do a reduce on min,max_level before supporting multiple proc */
    /* mpi_assert(proc_size == 1); */
    /* int num_of_levels = max_level + 1; */
    
    

    d4est_solver_multigrid_data_t* mg_data = d4est_solver_multigrid_data_init(
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




 
    /* d4est_solver_multigrid_logger_t* logger = d4est_solver_multigrid_logger_residual_init */
                                 /* ( */
                                 /* ); */

    /* d4est_solver_multigrid_element_data_updater_t* updater = d4est_solver_multigrid_element_data_updater_init */
                                                /* ( */
                                                 /* mg_data->num_of_levels, */
                                                 /* ghost, */
                                                 /* ghost_data, */
                                                 /* d4est_factors, */
                                                 /* problem_set_degrees_after_amr, */
                                                 /* NULL */
                                                /* ); */

    

    /* d4est_solver_multigrid_set_callbacks( */
                            /* mg_data, */
                            /* logger, */
                            /* NULL, */
                            /* updater */
    /* ); */
    
    d4est_krylov_pc_t* pc = d4est_krylov_pc_multigrid_create(mg_data, NULL);
    
    krylov_petsc_params_t krylov_petsc_params;
    krylov_petsc_input(p4est, input_file, "krylov_petsc", &krylov_petsc_params);


    prob_vecs.field_types = &field_type;
    prob_vecs.num_of_fields = 1;
      
    
    krylov_petsc_solve
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
       &krylov_petsc_params,
       pc
      );


    d4est_mesh_interpolate_data_t data;

    double R0 = ((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user)->R0;
    double R1 = ((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user)->R1;
    double R2 = ((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user)->R2;
    int compactify_inner_shell = ((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user)->compactify_inner_shell;
    int compactify_outer_shell = ((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user)->compactify_outer_shell;
    d4est_geometry_type_t geom_type =  d4est_geom->geom_type;
    
    data = d4est_mesh_interpolate_at_tree_coord(p4est, d4est_ops, d4est_geom, (double []){get_inverted_box_point(R0,0),.5,.5}, 12, prob_vecs.u,  1);
    point[0][iterations] = (data.err == 0) ? data.f_at_xyz : 0;
    point_err[0] = data.err;
    printf("1st point is at xyz = %.15f,%.15f,%.15f\n",data.xyz[0],data.xyz[1],data.xyz[2]);
    
    data = d4est_mesh_interpolate_at_tree_coord(p4est, d4est_ops, d4est_geom, (double []){get_inverted_box_point(R0,3),.5,.5}, 12, prob_vecs.u, 1);
    point[1][iterations] = (data.err == 0) ? data.f_at_xyz : 0;
    point_err[1] = data.err;
    printf("2nd point is at xyz = %.15f,%.15f,%.15f\n",data.xyz[0],data.xyz[1],data.xyz[2]);
    
    data =  d4est_mesh_interpolate_at_tree_coord(p4est, d4est_ops, d4est_geom, (double []){.5,.5,get_inverted_inner_wedge_point(R0,R1,10,compactify_inner_shell)}, 9, prob_vecs.u, 1);
    point[2][iterations] = (data.err == 0) ? data.f_at_xyz : 0;
    point_err[2] = data.err;
    printf("3rd point is at xyz = %.15f,%.15f,%.15f\n",data.xyz[0],data.xyz[1],data.xyz[2]);

    /* if (geom_type == GEOM_CUBED_SPHERE_7TREE){ */
    /* data =  d4est_mesh_interpolate_at_tree_coord(p4est, */
    /*                                              d4est_ops, */
    /*                                              d4est_geom, */
    /*                                              (double []){.5,.5,get_inverted_inner_wedge_point(R0, */
    /*                                                                                               R1, */
    /*                                                                                               (100 > R1) ? R1 : 100,compactify_inner_shell)}, */
    /*                                              3, */
    /*                                              prob_vecs.u, 1); */
    /* point[3][iterations] = (data.err == 0) ? data.f_at_xyz : 0; */
    /* point_err[3] = data.err; */
    /* printf("4th point is at xyz = %.15f,%.15f,%.15f\n",data.xyz[0],data.xyz[1],data.xyz[2]); */
    /* } */
    /* else if (geom_type == GEOM_CUBED_SPHERE_13TREE){ */
    data =  d4est_mesh_interpolate_at_tree_coord(p4est, d4est_ops, d4est_geom, (double []){.5,.5,get_inverted_outer_wedge_point(R1,R2, (100 > R2) ? R2 : 100, compactify_outer_shell)}, 3, prob_vecs.u, 1);
    point[3][iterations] = (data.err == 0) ? data.f_at_xyz : 0;
    point_err[3] = data.err;
    printf("4th point is at xyz = %.15f,%.15f,%.15f\n",data.xyz[0],data.xyz[1],data.xyz[2]);
    /* } */
    /* else { */
    /*   D4EST_ABORT("not support geom type"); */
    /* } */
    double* point0 = &point[0][0];
    double* point3 = &point[1][0];
    double* point10 = &point[2][0];
    double* point100 = &point[3][0];
    double* point0_diff = &point_diff[0][0];
    double* point3_diff = &point_diff[1][0];
    double* point10_diff = &point_diff[2][0];
    double* point100_diff = &point_diff[3][0];
    double* point0_spec_diff = &point_spec_diff[0][0];
    double* point3_spec_diff = &point_spec_diff[1][0];
    double* point10_spec_diff = &point_spec_diff[2][0];
    double* point100_spec_diff = &point_spec_diff[3][0];
    
    int global_nodes;
    sc_reduce(
              &prob_vecs.local_nodes,
              &global_nodes,
              1,
              sc_MPI_INT,
              sc_MPI_SUM,
              0,
              sc_MPI_COMM_WORLD
    );
    point_dof[iterations] = global_nodes;
    double* dof = &point_dof[0];
    double points_global [4];
    double points_local [4];
    points_local[0] = point[0][iterations];
    points_local[1] = point[1][iterations];
    points_local[2] = point[2][iterations];
    points_local[3] = point[3][iterations];

    sc_reduce
      (
       &points_local,
       &points_global,
       4,
       sc_MPI_DOUBLE,
       sc_MPI_MAX,
       0,
       sc_MPI_COMM_WORLD
      );

     
    if (p4est->mpirank == 0){
      for (int p = 0; p < 4; p++){
        point[p][iterations] = points_global[p];
        point_diff[p][iterations] = fabs(point[p][iterations] - point[p][iterations-1]);
      }
      point_spec_diff[0][iterations] = fabs(point[0][iterations] - 1.);
      point_spec_diff[1][iterations] = fabs(point[1][iterations] - 0.31622776601);
      point_spec_diff[2][iterations] = fabs(point[2][iterations] - 0.09950371902);
      point_spec_diff[3][iterations] = fabs(point[3][iterations] - 0.00999950003);
      
      DEBUG_PRINT_4ARR_DBL(dof, point0, point0_diff, point0_spec_diff, iterations+1);
      DEBUG_PRINT_4ARR_DBL(dof, point3, point3_diff, point3_spec_diff,iterations+1);
      DEBUG_PRINT_4ARR_DBL(dof, point10, point10_diff, point10_spec_diff,iterations+1);
      DEBUG_PRINT_4ARR_DBL(dof, point100, point100_diff, point100_spec_diff,iterations+1);
    }
    iterations++;
    

    d4est_krylov_pc_multigrid_destroy(pc);
    
    /* d4est_solver_multigrid_logger_residual_destroy(logger); */
    /* d4est_solver_multigrid_element_data_updater_destroy(updater, mg_data->num_of_levels); */
    d4est_solver_multigrid_data_destroy(mg_data);


    P4EST_FREE(estimator);


    
  }

  if (d4est_ghost_data != NULL){
    d4est_ghost_data_destroy(d4est_ghost_data);
    d4est_ghost_data = NULL;
  }
  
  d4est_mesh_data_destroy(d4est_factors_compactified);
  d4est_geometry_destroy(d4est_geom_compactified);  
  printf("[D4EST_INFO]: Starting garbage collection...\n");
  d4est_amr_destroy(d4est_amr);
  d4est_amr_destroy(d4est_amr_uniform);
  d4est_poisson_flux_destroy(flux_data_for_lhs);
  d4est_poisson_flux_destroy(flux_data_for_rhs);
  P4EST_FREE(error);
  P4EST_FREE(u_analytic);
  P4EST_FREE(prob_vecs.u);
  P4EST_FREE(prob_vecs.Au);
  P4EST_FREE(prob_vecs.rhs);
}
