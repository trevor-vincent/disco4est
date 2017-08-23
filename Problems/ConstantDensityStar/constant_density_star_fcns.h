#ifndef CONSTANT_DENSITY_STAR_FCNS_H
#define CONSTANT_DENSITY_STAR_FCNS_H 

#include <d4est_amr_smooth_pred.h>

typedef struct {

  /* set externally */
  double R;
  double cx;
  double cy;
  double cz;
  double rho0_div_rhoc;

  /* set internally */
  double alpha;
  double beta;
  double rho0;
  double C0;

} constant_density_star_params_t;

typedef struct {
  
  constant_density_star_params_t* constant_density_star_params;
  d4est_poisson_flux_data_t* flux_data_for_jac;
  d4est_poisson_flux_data_t* flux_data_for_res;
  d4est_amr_smooth_pred_params_t* smooth_pred_params;
  
} problem_ctx_t;

static
double u_alpha
(
 double x,
 double y,
 double z,
 void* user
)
{
  constant_density_star_params_t* params = user;
  double cx = params->cx;
  double cy = params->cy;
  double cz = params->cz;
  double alpha = params->alpha;
  double R = params->R;
  
  double dx = x - cx;
  double dy = y - cy;
  double dz = z - cz;
  double r2 = dx*dx + dy*dy + dz*dz;
  return sqrt(alpha*R)/sqrt(r2 + alpha*R*alpha*R);
}


static
int problem_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  constant_density_star_params_t* pconfig = (constant_density_star_params_t*)user;
  if (d4est_util_match_couple(section,"problem",name,"R")) {
    D4EST_ASSERT(pconfig->R == -1);
    pconfig->R = atof(value);
  } 
  else if (d4est_util_match_couple(section,"problem",name,"cx")) {
    D4EST_ASSERT(pconfig->cx == -1);
    pconfig->cx = atof(value);
  } 
  else if (d4est_util_match_couple(section,"problem",name,"cy")) {
    D4EST_ASSERT(pconfig->cy == -1);
    pconfig->cy = atof(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"cz")) {
    D4EST_ASSERT(pconfig->cz == -1);
    pconfig->cz = atof(value);
  }   
  else if (d4est_util_match_couple(section,"problem",name,"rho0_div_rhoc")) {
    D4EST_ASSERT(pconfig->rho0_div_rhoc == -1);
    pconfig->rho0_div_rhoc = atof(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

static
double solve_for_alpha
(
 double a,
 void* user
)
{
  constant_density_star_params_t* params = user;
  double R = params->R;
  double rho0 = params->rho0;
  
  double a5 = a*a*a*a*a;
  double opa2 = 1. + a*a;
  double f_of_a = a5/(opa2*opa2*opa2);
  double f2 = f_of_a*f_of_a;
  return rho0*R*R - (3./(2.*M_PI))*f2;
}

constant_density_star_params_t
constant_density_star_input
(
 const char* input_file
)
{
  constant_density_star_params_t input;

  input.R = -1;
  input.cx = -1;
  input.cy = -1;
  input.cz = -1;
  input.rho0_div_rhoc = -1;
  
  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("problem",input.R,-1);
  D4EST_CHECK_INPUT("problem",input.cx,-1);
  D4EST_CHECK_INPUT("problem",input.cy,-1);
  D4EST_CHECK_INPUT("problem",input.cz,-1);
  D4EST_CHECK_INPUT("problem",input.rho0_div_rhoc,-1);
  
  double R = input.R;
  double cx = input.cx;
  double cy = input.cy;
  double cz = input.cz;
  
  double alpha_crit = sqrt(5);
  double rhoc = (3./(2.*M_PI))*(1.0/(R*R))*((double)(5.*5.*5.*5.*5.)/(double)(6.*6.*6.*6.*6.*6.));  
  double rho0 = input.rho0_div_rhoc*rhoc;

  double C0 = pow(1./(2.*M_PI*rho0/3.),.25);
  double alpha = 386.266;

  input.rho0 = rho0;
  input.C0 = C0;

  
  printf("[CONSTANT_DENSITY_STAR_PARAMS]: R = %.25f\n", R);
  /* printf("[CONSTANT_DENSITY_STAR_PARAMS]: alpha = %.25f\n", alpha); */
  /* printf("[CONSTANT_DENSITY_STAR_PARAMS]: beta = %.25f\n", beta); */
  printf("[CONSTANT_DENSITY_STAR_PARAMS]: rho0 = %.25f\n", rho0);
  printf("[CONSTANT_DENSITY_STAR_PARAMS]: C0 = %.25f\n", C0);
  printf("[CONSTANT_DENSITY_STAR_PARAMS]: cx = %.25f\n", cx);
  printf("[CONSTANT_DENSITY_STAR_PARAMS]: cy = %.25f\n", cy);
  printf("[CONSTANT_DENSITY_STAR_PARAMS]: cz = %.25f\n", cz);
  printf("[CONSTANT_DENSITY_STAR_PARAMS]: rho0/rhoc = %.25f\n", input.rho0_div_rhoc);
  
  D4EST_ASSERT( !d4est_util_bisection(solve_for_alpha, alpha_crit, 1000*alpha_crit, DBL_EPSILON, 100000, &alpha, &input));

  double u_alpha_at_R = sqrt(alpha*R)/sqrt(R*R + alpha*R*alpha*R);
  double beta = R*(C0*u_alpha_at_R - 1.);

  input.alpha = alpha;
  input.beta = beta;
  
  D4EST_ASSERT(
             (C0*u_alpha(R + cx,cy,cz, &input) == 1. + beta/R)
             &&
             (C0*u_alpha(cx,R + cy,cz, &input) == 1. + beta/R)
             &&
             (C0*u_alpha(cx,cy,R + cz, &input) == 1. + beta/R)
            );

  /* printf("[CONSTANT_DENSITY_STAR_PARAMS]: R = %.25f\n", R); */
  printf("[CONSTANT_DENSITY_STAR_PARAMS]: alpha = %.25f\n", alpha);
  printf("[CONSTANT_DENSITY_STAR_PARAMS]: beta = %.25f\n", beta);
  /* printf("[CONSTANT_DENSITY_STAR_PARAMS]: rho0 = %.25f\n", rho0); */
  /* printf("[CONSTANT_DENSITY_STAR_PARAMS]: C0 = %.25f\n", C0); */
  /* printf("[CONSTANT_DENSITY_STAR_PARAMS]: cx = %.25f\n", cx); */
  /* printf("[CONSTANT_DENSITY_STAR_PARAMS]: cy = %.25f\n", cy); */
  /* printf("[CONSTANT_DENSITY_STAR_PARAMS]: cz = %.25f\n", cz); */
  /* printf("[CONSTANT_DENSITY_STAR_PARAMS]: rho0/rhoc = %.25f\n", input.rho0_div_rhoc); */
  
  return input;
}

static
double psi_fcn
(
 double x,
 double y,
 double z,
 void* user
)
{
  constant_density_star_params_t* params = user;
  double cx = params->cx;
  double cy = params->cy;
  double cz = params->cz;
  double beta = params->beta;
  double C0 = params->C0;
  double R = params->R;
  
  double dx = x - cx;
  double dy = y - cy;
  double dz = z - cz;
  double r2 = dx*dx + dy*dy + dz*dz;
  if (r2 > R*R)
    return 1. + beta/sqrt(r2);
  else
    return C0*u_alpha(x,y,z, user);
}

static
double rho_fcn
(
 double x,
 double y,
 double z,
 void* user
)
{
  constant_density_star_params_t* params = user;
  double cx = params->cx;
  double cy = params->cy;
  double cz = params->cz;
  double R = params->R;
  
  double dx = x - cx;
  double dy = y - cy;
  double dz = z - cz;
  double r2 = dx*dx + dy*dy + dz*dz;
  if (r2 > R*R){
    return 0.;
  }
  else{
    return params->rho0;
  }
}

static
double constant_density_star_analytic_solution
(
 double x,
 double y,
 double z,
 void* user
)
{
  /* psi = u + 1 */
  return psi_fcn(x,y,z,user) - 1;
}

static
double constant_density_star_boundary_fcn
(
 double x,
 double y,
 double z,
 void* user
)
{
  return zero_fcn(x,y,z,user);
}


static
double neg_10pi_rho_up1_neg4
(
 double x,
 double y,
 double z,
 double u,
 void* ctx
)
{
  return (-10.*M_PI)*rho_fcn(x,y,z,ctx)*(u+1)*(u+1)*(u+1)*(u+1);
}

static
double neg_2pi_rho_up1_neg5
(
 double x,
 double y,
 double z,
 double u,
 void* ctx
)
{
  return (-2.*M_PI)*rho_fcn(x,y,z,ctx)*(u+1)*(u+1)*(u+1)*(u+1)*(u+1);
}


static void
constant_density_star_build_residual_add_nonlinear_term
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
  double* M_neg_2pi_rho_up1_neg5_vec= P4EST_ALLOC(double, prob_vecs->local_nodes);

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

        d4est_quadrature_apply_fofufofvlj
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           &mesh_object,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &prob_vecs->u[ed->nodal_stride],
           NULL,
           ed->deg,
           ed->J_quad,
           ed->xyz_quad,
           ed->deg_vol_quad,
           &M_neg_2pi_rho_up1_neg5_vec[ed->nodal_stride],
           neg_2pi_rho_up1_neg5,
           user,
           NULL,
           NULL
          );

        
      }
    }

  d4est_linalg_vec_axpy(1.0, M_neg_2pi_rho_up1_neg5_vec, prob_vecs->Au, prob_vecs->local_nodes);

  P4EST_FREE(M_neg_2pi_rho_up1_neg5_vec); 
}

static
void
constant_density_star_build_residual
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
  problem_ctx_t* ctx = user;
  d4est_poisson_flux_data_t* flux_data = ctx->flux_data_for_res;
  
  d4est_poisson_apply_aij(p4est,
                          ghost,
                          ghost_data,
                          prob_vecs,
                          flux_data,
                          d4est_ops,
                          d4est_geom,
                          d4est_quad
                         );


  constant_density_star_build_residual_add_nonlinear_term
    (
     p4est,
     ghost,
     ghost_data,
     prob_vecs,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     user
    );
}


static
void constant_density_star_apply_jac_add_nonlinear_term
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
  double* M_neg_10pi_rho_up1_neg4_of_u0_u_vec = P4EST_ALLOC(double, prob_vecs->local_nodes);


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

        d4est_quadrature_apply_fofufofvlilj
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           &mesh_object,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &prob_vecs->u[ed->nodal_stride],
           &prob_vecs->u0[ed->nodal_stride],
           NULL,
           ed->deg,
           ed->xyz_quad,
           ed->J_quad,
           ed->deg_vol_quad,
           &M_neg_10pi_rho_up1_neg4_of_u0_u_vec[ed->nodal_stride],
           neg_10pi_rho_up1_neg4,
           user,
           NULL,
           NULL
          );
        
      }
    }

  d4est_linalg_vec_axpy(1.0, M_neg_10pi_rho_up1_neg4_of_u0_u_vec, prob_vecs->Au, prob_vecs->local_nodes);
  P4EST_FREE(M_neg_10pi_rho_up1_neg4_of_u0_u_vec);
}

static
void constant_density_star_apply_jac
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
  problem_ctx_t* ctx = user;
  d4est_poisson_flux_data_t* flux_data = ctx->flux_data_for_jac;
  
  d4est_poisson_apply_aij(p4est,
                          ghost,
                          ghost_data,
                          prob_vecs,
                          flux_data,
                          d4est_ops,
                          d4est_geom,
                          d4est_quad
                         );
  
 constant_density_star_apply_jac_add_nonlinear_term
    (
     p4est,
     ghost,
     ghost_data,
     prob_vecs,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     user
    );
  
}


static double
constant_density_star_initial_guess
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return 1.;
}


#endif
