#ifndef CDS_FCNS_H
#define CDS_FCNS_H 


static const double pi = 3.1415926535897932384626433832795;



typedef struct {

  /* set internally */
  double R;
  double alpha;
  double beta;
  double rho0;
  double C0;
  double cx;
  double cy;
  double cz;

  /* set in input file */
  int num_unifrefs;
  int num_of_amr_levels;
  int deg_R0;
  int deg_quad_R0;
  int deg_R1;
  int deg_quad_R1;
  int deg_R2;
  int deg_quad_R2; 
  double rho0_div_rhoc;
  double ip_flux_penalty;

  ip_flux_params_t* ip_flux_params;

} constantdensitystar_params_t;


static
double u_alpha
(
 double x,
 double y,
 double z,
 void* user
)
{
  constantdensitystar_params_t* params = user;
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
  constantdensitystar_params_t* pconfig = (constantdensitystar_params_t*)user;
  if (util_match_couple(section,"amr",name,"num_of_amr_levels")) {
    mpi_assert(pconfig->num_of_amr_levels == -1);
    pconfig->num_of_amr_levels = atoi(value);
  }
  else if (util_match_couple(section,"amr",name,"num_unifrefs")) {
    mpi_assert(pconfig->num_unifrefs == -1);
    pconfig->num_unifrefs = atoi(value);
  }
  else if (util_match_couple(section,"flux",name,"ip_flux_penalty")) {
    mpi_assert(pconfig->ip_flux_penalty == -1);
    pconfig->ip_flux_penalty = atof(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_R0")) {
    mpi_assert(pconfig->deg_R0 == -1);
    pconfig->deg_R0 = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_quad_R0")) {
    mpi_assert(pconfig->deg_quad_R0 == -1);
    pconfig->deg_quad_R0 = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_R1")) {
    mpi_assert(pconfig->deg_R1 == -1);
    pconfig->deg_R1 = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_quad_R1")) {
    mpi_assert(pconfig->deg_quad_R1 == -1);
    pconfig->deg_quad_R1 = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_R2")) {
    mpi_assert(pconfig->deg_R2 == -1);
    pconfig->deg_R2 = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_quad_R2")) {
    mpi_assert(pconfig->deg_quad_R2 == -1);
    pconfig->deg_quad_R2 = atoi(value);
  }  
  else if (util_match_couple(section,"problem",name,"rho0_div_rhoc")) {
    mpi_assert(pconfig->rho0_div_rhoc == -1);
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
  constantdensitystar_params_t* params = user;
  double R = params->R;
  double rho0 = params->rho0;
  
  double a5 = a*a*a*a*a;
  double opa2 = 1. + a*a;
  double f_of_a = a5/(opa2*opa2*opa2);
  double f2 = f_of_a*f_of_a;
  return rho0*R*R - (3./(2.*pi))*f2;
}

constantdensitystar_params_t
constantdensitystar_input
(
 const char* input_file,
 d4est_geometry_t* d4est_geom
)
{
  
  constantdensitystar_params_t input;

  input.num_unifrefs = -1;
  input.num_of_amr_levels = -1;
  input.ip_flux_penalty = -1;
  input.deg_R0 = -1;
  input.deg_quad_R0 = -1;
  input.deg_R1 = -1;
  input.deg_quad_R1 = -1;
  input.deg_R2 = -1;
  input.deg_quad_R2 = -1; 
 
  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  D4EST_CHECK_INPUT("problem", input.num_unifrefs, -1);
  D4EST_CHECK_INPUT("problem", input.num_of_amr_levels, -1);
  D4EST_CHECK_INPUT("problem", input.ip_flux_penalty, -1);
  D4EST_CHECK_INPUT("problem", input.deg_R0, -1);
  D4EST_CHECK_INPUT("problem", input.deg_quad_R0, -1);
  D4EST_CHECK_INPUT("problem", input.deg_R1, -1);
  D4EST_CHECK_INPUT("problem", input.deg_quad_R1, -1);
  D4EST_CHECK_INPUT("problem", input.deg_R2, -1);
  D4EST_CHECK_INPUT("problem", input.deg_quad_R2, -1);

  ip_flux_params_t ip_flux_params;
  ip_flux_params.ip_flux_penalty_prefactor = input.ip_flux_penalty;
  ip_flux_params.ip_flux_penalty_calculate_fcn = sipg_flux_vector_calc_penalty_maxp2_over_minh;
  input.ip_flux_params = &ip_flux_params;
  
  double R1 = ((d4est_geometry_sphere_attr_t*)(d4est_geom->p4est_geom->user))->R1;
  double R = R1;
  double alpha_crit = sqrt(5);
  double rhoc = (3./(2.*pi))*(1.0/(R*R))*((double)(5.*5.*5.*5.*5.)/(double)(6.*6.*6.*6.*6.*6.));  
  double rho0 = input.rho0_div_rhoc*rhoc;
  double cx = 0;
  double cy = 0;
  double cz = 0;
  double C0 = pow(1./(2.*pi*rho0/3.),.25);
  double alpha = 386.266;

  mpi_assert( !util_bisection(solve_for_alpha, alpha_crit, 1000*alpha_crit, DBL_EPSILON, 100000, &alpha, &input));

  double u_alpha_at_R = sqrt(alpha*R)/sqrt(R*R + alpha*R*alpha*R);
  double beta = R*(C0*u_alpha_at_R - 1.);
  
  mpi_assert(
             (C0*u_alpha(R + 0.,0.,0., &input) == 1. + beta/R)
             &&
             (C0*u_alpha(0.,R + 0.,0., &input) == 1. + beta/R)
             &&
             (C0*u_alpha(0.,0.,R + 0., &input) == 1. + beta/R)
            );

  input.R = R;
  input.alpha = alpha;
  input.beta = beta;
  input.rho0 = rho0;
  input.C0 = C0;
  input.cx = cx;
  input.cy = cy;
  input.cz = cz;

  printf("[CDS_PARAMS]: R = %.25f\n", R);
  printf("[CDS_PARAMS]: alpha = %.25f\n", alpha);
  printf("[CDS_PARAMS]: beta = %.25f\n", beta);
  printf("[CDS_PARAMS]: rho0 = %.25f\n", rho0);
  printf("[CDS_PARAMS]: C0 = %.25f\n", C0);
  printf("[CDS_PARAMS]: cx = %.25f\n", cx);
  printf("[CDS_PARAMS]: cy = %.25f\n", cy);
  printf("[CDS_PARAMS]: cz = %.25f\n", cz);
  printf("[CDS_PARAMS]: rho0/rhoc = %.25f\n", input.rho0_div_rhoc);
  
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
  constantdensitystar_params_t* params = user;
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
  constantdensitystar_params_t* params = user;
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
double analytic_solution_fcn
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
double boundary_fcn
(
 double x,
 double y,
 double z
)
{
  /* return psi_fcn(x,y,z); */
  return zero_fcn(x,y,z);
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
  return (-10.*pi)*rho_fcn(x,y,z,ctx)*(u+1)*(u+1)*(u+1)*(u+1);
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
  return (-2.*pi)*rho_fcn(x,y,z,ctx)*(u+1)*(u+1)*(u+1)*(u+1)*(u+1);
}

static
void
cds_build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
  constantdensitystar_params_t* ctx = (constantdensitystar_params_t*)prob_vecs->user;
  
  prob_vecs->curved_scalar_flux_fcn_data = curved_Gauss_primal_sipg_flux_dirichlet_fetch_fcns
                                           (boundary_fcn, ctx->ip_flux_params);
  curved_poisson_operator_primal_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom);

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
        int deg_nonlinear = ed->deg;// + ctx->deg_offset_for_nonlinear_quad;

        d4est_operators_apply_fofufofvlj_Gaussnodes
          (
           d4est_ops,
           &prob_vecs->u[ed->nodal_stride],
           NULL,
           ed->deg,
           ed->J_quad,
           ed->xyz_quad,
           ed->deg_quad,
           (P4EST_DIM),
           &M_neg_2pi_rho_up1_neg5_vec[ed->nodal_stride],
           neg_2pi_rho_up1_neg5,
           ctx,
           NULL,
           NULL
          );
      }
    }

  linalg_vec_axpy(1.0, M_neg_2pi_rho_up1_neg5_vec, prob_vecs->Au, prob_vecs->local_nodes);

  P4EST_FREE(M_neg_2pi_rho_up1_neg5_vec); 
}

static
void cds_apply_jac
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
  constantdensitystar_params_t* ctx = (constantdensitystar_params_t*)prob_vecs->user;

  /* apply jac must always have zero boundary conditions for du in dirichlet problems */
  prob_vecs->curved_scalar_flux_fcn_data = curved_Gauss_primal_sipg_flux_dirichlet_fetch_fcns
                                           (zero_fcn, ctx->ip_flux_params);
  curved_poisson_operator_primal_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom);
  
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
        int deg_nonlinear = ed->deg;// + ctx->deg_offset_for_nonlinear_quad;
        d4est_operators_apply_fofufofvlilj_Gaussnodes
          (
           d4est_ops,
           &prob_vecs->u[ed->nodal_stride],
           &prob_vecs->u0[ed->nodal_stride],
           NULL,
           ed->deg,
           ed->J_quad,
           ed->xyz_quad,
           deg_nonlinear,
           (P4EST_DIM),
           &M_neg_10pi_rho_up1_neg4_of_u0_u_vec[ed->nodal_stride],
           neg_10pi_rho_up1_neg4,
           ctx,
           NULL,
           NULL
          );
      }
    }

  linalg_vec_axpy(1.0, M_neg_10pi_rho_up1_neg4_of_u0_u_vec, prob_vecs->Au, prob_vecs->local_nodes);
  P4EST_FREE(M_neg_10pi_rho_up1_neg4_of_u0_u_vec);
}


#endif
