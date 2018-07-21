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
#include <d4est_solver_newton_petsc.h>
#include <d4est_solver_krylov_petsc.h>
#include <d4est_krylov_pc_multigrid.h>
#include <d4est_solver_multigrid_logger_residual.h>
#include <d4est_solver_multigrid_element_data_updater.h>
#include <d4est_solver_multigrid.h>
#include <d4est_util.h>
#include <time.h>
#include <zlog.h>
#include "poisson_lorentzian_fcns.h"


typedef struct {
  
  int use_dirichlet;
  
} poisson_lorentzian_init_params_t;


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

  d4est_field_type_t field_type = VOLUME_NODAL;

  
  /* d4est_poisson_build_rhs_with_strong_bc( */
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
  /*   poisson_lorentzian_rhs_fcn, */
  /*   INIT_FIELD_ON_LOBATTO, */
  /*   &ctx, */
  /*   0 */
  /* ); */



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

  
  for (int level = 0; level < d4est_amr->num_of_amr_steps + 1; level++) {


    d4est_ghost_data_t* d4est_ghost_data = d4est_ghost_data_init(p4est,
                                                                 *d4est_ghost,
                                                                 &field_type,
                                                                 1);


    
    d4est_poisson_build_rhs_with_strong_bc(
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


    
    // Setup d4est_solver_multigrid
    d4est_krylov_pc_t* pc = NULL;
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

    pc = d4est_krylov_pc_multigrid_create(mg_data, NULL);

    // Krylov PETSc solve
    
    krylov_petsc_params_t krylov_petsc_params;
    krylov_petsc_input(p4est, input_file, "krylov_petsc", &krylov_petsc_params);

    prob_vecs.field_types = &field_type;
    prob_vecs.num_of_fields = 1;
      
    krylov_petsc_solve(
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


    double* rhs_fcn = P4EST_ALLOC(double, prob_vecs.local_nodes);
    d4est_mesh_init_field(
      p4est,
      rhs_fcn,
      poisson_lorentzian_rhs_fcn,
      d4est_ops, // unnecessary?
      d4est_geom, // unnecessary?
      d4est_factors,
      INIT_FIELD_ON_LOBATTO,
      NULL
    );

    
    
    // Compute errors between numerical and analytical field values
    double* error = P4EST_ALLOC(double, prob_vecs.local_nodes);
    double* error_nofabs = P4EST_ALLOC(double, prob_vecs.local_nodes);

    for (int i = 0; i < prob_vecs.local_nodes; i++){
      error_nofabs[i] = prob_vecs.u[i] - u_analytic[i];
    }
    if (init_params.use_dirichlet){
    poisson_lorentzian_apply_lhs_with_bc
      (
       p4est,
       *d4est_ghost,
       d4est_ghost_data,
       &prob_vecs,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       &ctx
      );
    }
    else {
    poisson_lorentzian_apply_lhs
      (
       p4est,
       *d4est_ghost,
       d4est_ghost_data,
       &prob_vecs,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       &ctx
      );
    }

    double* laplace_u_analytic = P4EST_ALLOC(double, prob_vecs.local_nodes);
    double* error_lap_u = P4EST_ALLOC(double, prob_vecs.local_nodes);


    d4est_elliptic_data_t prob_vecs_analytic;
    d4est_elliptic_data_copy_ptrs(&prob_vecs, &prob_vecs_analytic);

    prob_vecs_analytic.Au = laplace_u_analytic;
    prob_vecs_analytic.u = u_analytic;

    if (init_params.use_dirichlet){
    poisson_lorentzian_apply_lhs_with_bc
      (
       p4est,
       *d4est_ghost,
       d4est_ghost_data,
       &prob_vecs_analytic,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       &ctx
      );
    }
    else {
    poisson_lorentzian_apply_lhs
      (
       p4est,
       *d4est_ghost,
       d4est_ghost_data,
       &prob_vecs_analytic,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       &ctx
      );
    }
    double* Minv_Au = P4EST_ALLOC(double, prob_vecs.local_nodes);
    double* Minv_Au_analytic = P4EST_ALLOC(double, prob_vecs.local_nodes);
      
    
    d4est_mesh_apply_invM_on_field
      (
       p4est,
       d4est_ops,
       d4est_factors,
       prob_vecs.Au,
       Minv_Au
      );

    d4est_mesh_apply_invM_on_field
      (
       p4est,
       d4est_ops,
       d4est_factors,
       prob_vecs_analytic.Au,
       Minv_Au_analytic
      );
    
    
    
    for (int n = 0; n < prob_vecs.local_nodes; n++){
      error_lap_u[n] = Minv_Au[n] - Minv_Au_analytic[n];  
    }
    
    /* d4est_mesh_debug_boundary_elements(p4est, */
    /*                                    d4est_ops, */
    /*                                    d4est_factors, */
    /*                                    (const char * []){"u","u_analytic","laplace_u","laplace_analytic", "error", NULL}, */
    /*                                    (double* []){prob_vecs.u, prob_vecs_analytic.u, prob_vecs.Au, prob_vecs_analytic.Au, error_nofabs, NULL}, */
    /*                                    prob_vecs.local_nodes */
    /*                                   ); */



    
    /* d4est_mesh_debug_boundary_elements(p4est, */
    /*                                    d4est_ops, */
    /*                                    d4est_factors, */
    /*                                    (const char * []){"u","u_analytic","lap_u","lap_u_anal", "error_u", "error_lap_u", NULL}, */
    /*                                    (double* []){prob_vecs.u, prob_vecs_analytic.u, prob_vecs.Au, prob_vecs_analytic.Au, error_nofabs, error_lap_u, NULL}, */
    /*                                    prob_vecs.local_nodes */
    /*                                   ); */




    
    /* d4est_mesh_debug_boundary_elements(p4est, */
    /*                                    d4est_ops, */
    /*                                    d4est_factors, */
    /*                                    (const char * []){"u","u_analytic","Minv_Au", "Minv_Au_analytic", "rhs_fcn",  NULL}, */
    /*                                    (double* []){prob_vecs.u, prob_vecs_analytic.u, Minv_Au, Minv_Au_analytic, rhs_fcn, NULL}, */
    /*                                    prob_vecs.local_nodes */
    /*                                   ); */

    


    P4EST_FREE(error_lap_u);
    P4EST_FREE(rhs_fcn);
    P4EST_FREE(Minv_Au);
    P4EST_FREE(Minv_Au_analytic);
    
    /* d4est_mesh_debug_boundary_elements(p4est, */
    /*                                    d4est_ops, */
    /*                                    d4est_factors, */
    /*                                    (const char * []){"error", NULL}, */
    /*                                    (double* []){error_nofabs, NULL}, */
    /*                                    prob_vecs.local_nodes */
    /*                                   ); */

    
    P4EST_FREE(error_nofabs);
    P4EST_FREE(laplace_u_analytic);
    
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
        NULL,
        NULL,
        d4est_ops,
        d4est_amr,
        &prob_vecs.u,
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
  d4est_poisson_flux_destroy(flux_data_for_lhs);
  d4est_poisson_flux_destroy(flux_data_for_rhs);
  P4EST_FREE(prob_vecs.u);
  P4EST_FREE(prob_vecs.Au);
  P4EST_FREE(prob_vecs.rhs);
}
