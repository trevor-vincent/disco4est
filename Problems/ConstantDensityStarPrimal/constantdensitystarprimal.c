
/* #define TEST_SYM  /\* Test if matrix is symmetric *\/ */
/* #define NDEBUG /\* TURN DEBUG OFF/ON *\/ */

#define NDEBUG
#define _GNU_SOURCE
#include <sc_reduce.h>
#include <pXest.h>
#include <util.h>
#include <d4est_linalg.h>

#include <d4est_element_data.h>
#include <curved_gauss_primal_sipg_flux_fcns.h>
#include <sipg_flux_vector_fcns.h>
#include <curved_poisson_operator_primal.h>
#include <problem.h>
#include <problem_data.h>
#include <problem_weakeqn_ptrs.h>
#include <poisson_operator.h>
#include <hp_amr_curved_smooth_pred.h>
#include <curved_bi_estimator.h>
#include <estimator_stats.h>
#include <curved_bi_estimator_flux_fcns.h>
#include <bi_estimator_flux_fcns.h>
#include <hp_amr.h>
#include <newton_petsc.h>
#include <d4est_vtk.h>
#include <dg_norm.h>
#include <ini.h>
#include "time.h"

static const double pi = 3.1415926535897932384626433832795;
double R;
double alpha;
double beta;
double rho0;
double C0;
double cx;
double cy;
double cz;


typedef struct {

  int deg_off_set_for_nonlinear_quad;
  ip_flux_params_t* ip_flux_params;
  
} problem_ctx_t;

static
double u_alpha
(
 double x,
 double y,
 double z
)
{
  double dx = x - cx;
  double dy = y - cy;
  double dz = z - cz;
  double r2 = dx*dx + dy*dy + dz*dz;
  return sqrt(alpha*R)/sqrt(r2 + alpha*R*alpha*R);
}

static
double solve_for_alpha
(
 double a
)
{
  double a5 = a*a*a*a*a;
  double opa2 = 1. + a*a;
  double f_of_a = a5/(opa2*opa2*opa2);
  double f2 = f_of_a*f_of_a;
  return rho0*R*R - (3./(2.*pi))*f2;
}

static
double psi_fcn
(
 double x,
 double y,
 double z
)
{
  double dx = x - cx;
  double dy = y - cy;
  double dz = z - cz;
  double r2 = dx*dx + dy*dy + dz*dz;
  if (r2 > R*R)
    return 1. + beta/sqrt(r2);
  else
    return C0*u_alpha(x,y,z);
}

static
double rho_fcn
(
 double x,
 double y,
 double z
)
{
  double dx = x - cx;
  double dy = y - cy;
  double dz = z - cz;
  double r2 = dx*dx + dy*dy + dz*dz;
  if (r2 > R*R)
    return 0.;
  else
    return rho0;
}

static
double analytical_solution_fcn
(
 double x,
 double y,
 double z
)
{
  return psi_fcn(x,y,z);
}

static
double boundary_fcn
(
 double x,
 double y,
 double z
)
{
  return psi_fcn(x,y,z);
}


static
double neg_10pi_rho_up1_neg4
(
 double x,
 double y,
 double z,
 double psi,
 void* ctx
)
{
  return (-10.*pi)*rho_fcn(x,y,z)*(psi)*(psi)*(psi)*(psi);
}

static
double neg_2pi_rho_up1_neg5
(
 double x,
 double y,
 double z,
 double psi,
 void* ctx
)
{
  return (-2.*pi)*rho_fcn(x,y,z)*(psi)*(psi)*(psi)*(psi)*(psi);
}


static
void
build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
  problem_ctx_t* ctx = (problem_ctx_t*)prob_vecs->user;
  
  prob_vecs->curved_scalar_flux_fcn_data = curved_gauss_primal_sipg_flux_dirichlet_fetch_fcns
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
        int deg_nonlinear = ed->deg;

        d4est_operators_apply_fofufofvlj_gaussnodes
          (
           d4est_ops,
           &prob_vecs->u[ed->nodal_stride],
           NULL,
           ed->deg,
           ed->J_quad,
           ed->xyz_quad,
           deg_nonlinear,
           (P4EST_DIM),
           &M_neg_2pi_rho_up1_neg5_vec[ed->nodal_stride],
           neg_2pi_rho_up1_neg5,
           NULL,
           NULL,
           NULL
          );
      }
    }

  d4est_linalg_vec_axpy(1.0, M_neg_2pi_rho_up1_neg5_vec, prob_vecs->Au, prob_vecs->local_nodes);

  P4EST_FREE(M_neg_2pi_rho_up1_neg5_vec); 
}

static
void apply_jac
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
  problem_ctx_t* ctx = (problem_ctx_t*)prob_vecs->user;

  /* apply jac must always have zero boundary conditions for du in dirichlet problems */
  prob_vecs->curved_scalar_flux_fcn_data = curved_gauss_primal_sipg_flux_dirichlet_fetch_fcns
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
        d4est_operators_apply_fofufofvlilj_gaussnodes
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
           NULL,
           NULL,
           NULL
          );
      }
    }


  d4est_linalg_vec_axpy(1.0, M_neg_10pi_rho_up1_neg4_of_u0_u_vec, prob_vecs->Au, prob_vecs->local_nodes);

  P4EST_FREE(M_neg_10pi_rho_up1_neg4_of_u0_u_vec);

}


p4est_connectivity_t*
problem_build_conn()
{
  return p8est_connectivity_new_unitcube();
}

p4est_geometry_t*
problem_build_geom
(
 p4est_connectivity_t* conn
)
{
  return NULL;
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
     sizeof(element_data_t),
     NULL,
     NULL
    );
}

typedef struct {

  int endlevel;
  int degree;
  double gamma_h;
  double gamma_p;
  double domain_size;
  double rho0_div_rhoc;
  double ip_flux_penalty;
  int percentile;
  int use_gauss_quad;
  int degmax;
  int deg_offset_for_nonlinear_quad;
  KSPType krylov_type;
  int amr_inflation_size;
  
  int count;
  
} problem_input_t;

static
int
amr_mark_element(
                 p4est_t* p4est,
                 double eta2,
                 estimator_stats_t** stats,
                 d4est_element_data_t* elem_data,
                 void* user
                )
{
  problem_input_t* input = user;
  if (p4est->local_num_quadrants*p4est->mpisize < input->amr_inflation_size){
    double eta2_percentile
      = estimator_stats_get_percentile(*stats,25);
    return (eta2 >= eta2_percentile);
  }
  else{
    double eta2_percentile
      = estimator_stats_get_percentile(*stats,input->percentile);
    return (eta2 >= eta2_percentile);
  }
}

static
gamma_params_t
amr_set_element_gamma
(
 p4est_t* p4est,
 double eta2,
 estimator_stats_t** stats,
 d4est_element_data_t* elem_data,
 void* user
)
{
  problem_input_t* input = user;
  gamma_params_t gamma_hpn;
  gamma_hpn.gamma_h = input->gamma_h;
  gamma_hpn.gamma_p = input->gamma_p;
  gamma_hpn.gamma_n = 1.;

  return gamma_hpn;
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
  problem_input_t* pconfig = (problem_input_t*)user;
  if (util_match_couple(section,"amr",name,"amr_levels")) {
    mpi_assert(pconfig->endlevel == -1);
    pconfig->endlevel = atoi(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"amr",name, "initial_degree")) {
    mpi_assert(pconfig->degree == -1);
    pconfig->degree = atoi(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"amr",name,"percentile")) {
    mpi_assert(pconfig->percentile == -1);
    pconfig->percentile = atoi(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"amr",name,"degmax")) {
    mpi_assert(pconfig->degmax == -1);
    pconfig->degmax = atoi(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"amr",name,"gamma_h")) {
    mpi_assert(pconfig->gamma_h == -1);
    pconfig->gamma_h = atof(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"amr",name,"gamma_p")) {
    mpi_assert(pconfig->gamma_p == -1);
    pconfig->gamma_p = atof(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"flux",name,"ip_flux_penalty")) {
    mpi_assert(pconfig->ip_flux_penalty == -1);
    pconfig->ip_flux_penalty = atof(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"problem",name,"rho0_div_rhoc")) {
    mpi_assert(pconfig->rho0_div_rhoc == -1);
    pconfig->rho0_div_rhoc = atof(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"problem",name,"domain_size")) {
    mpi_assert(pconfig->domain_size == -1);
    pconfig->domain_size = atof(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"problem",name,"use_gauss_quad")) {
    mpi_assert(pconfig->use_gauss_quad == -1);
    pconfig->use_gauss_quad = atoi(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"problem",name,"deg_offset_for_nonlinear_quad")) {
    mpi_assert(pconfig->deg_offset_for_nonlinear_quad == -1);
    pconfig->deg_offset_for_nonlinear_quad = atoi(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"amr",name,"amr_inflation_size")) {
    mpi_assert(pconfig->amr_inflation_size == -1);
    pconfig->amr_inflation_size = atoi(value);
    pconfig->count += 1;
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
  int num_of_options = 12;
  
  problem_input_t input;
  input.degree = -1;
  input.domain_size = -1;
  input.endlevel = -1;
  input.gamma_h = -1;
  input.gamma_p = -1;
  input.ip_flux_penalty = -1;
  input.percentile = -1;
  input.rho0_div_rhoc = -1;
  input.use_gauss_quad = -1;
  input.deg_offset_for_nonlinear_quad = -1;
  input.amr_inflation_size = -1;
  input.degmax = -1;
  input.count = 0;
  
  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  mpi_assert(input.count == num_of_options);
  return input;
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
  elem_data->deg = input->degree;
  elem_data->deg_quad = input->degree;
}

void
problem_keep_degrees
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  problem_input_t* input = user_ctx;
  elem_data->deg_quad = elem_data->deg;
}

void
problem_save_to_vtk
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 double* u,
 double* u_analytic,
 double* u_error,
 int level,
 int with_eta
)
{
   int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
   double* eta_array = P4EST_ALLOC(double, p4est->local_num_quadrants);
   int vtk_nodes = 0;
     
    int stride = 0;
    for (p4est_topidx_t tt = p4est->first_local_tree;
         tt <= p4est->last_local_tree;
         ++tt)
      {
        p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
        sc_array_t* tquadrants = &tree->quadrants;
        int Q = (p4est_locidx_t) tquadrants->elem_count;
        for (int q = 0; q < Q; ++q) {
          /* k++; */
          p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
          d4est_element_data_t* ed = quad->p.user_data;
          deg_array[stride] = ed->deg;
          eta_array[stride] = ed->local_estimator;
          vtk_nodes = util_int_pow_int(deg_array[stride], (P4EST_DIM))*(P4EST_CHILDREN);
          stride++;
        }
      }

    p4est_connectivity_t* conn_vtk = p8est_connectivity_new_unitcube();
    p4est_geometry_t* geom_vtk = p8est_geometry_new_connectivity(conn_vtk);


    char sol_save_as [500];
    if (with_eta)
      sprintf(sol_save_as, "%s_hp_amr_level_%d_sols_witheta", "cds", level);
    else
      sprintf(sol_save_as, "%s_hp_amr_level_%d_sols_noeta", "cds", level);
    
    d4est_vtk_context_t* vtk_ctx = d4est_vtk_dg_context_new(p4est, d4est_ops, sol_save_as);
    d4est_vtk_context_set_geom(vtk_ctx, geom_vtk);
    d4est_vtk_context_set_scale(vtk_ctx, .99);
    d4est_vtk_context_set_deg_array(vtk_ctx, deg_array);
    vtk_ctx = d4est_vtk_write_dg_header(vtk_ctx, d4est_ops);    
    vtk_ctx = d4est_vtk_write_dg_point_dataf(vtk_ctx,
                                             3,
                                             0,
                                             "u",
                                             u,
                                             "u_analytic",
                                             u_analytic,
                                             "u_error",
                                             u_error,
                                             vtk_ctx
                                            );

    if (with_eta){
      vtk_ctx = d4est_vtk_write_dg_cell_dataf
                (
                 vtk_ctx,
                 1,
                 1,
                 1,
                 0,
                 1,
                 1,
                 0,
                 "eta",
                 eta_array,
                 vtk_ctx
                );
    }
    else {
      vtk_ctx = d4est_vtk_write_dg_cell_dataf
                (
                 vtk_ctx,
                 1,
                 1,
                 1,
                 0,
                 1,
                 0,
                 0,
                 vtk_ctx
                );
    }
    
    d4est_vtk_write_footer(vtk_ctx);
    P4EST_FREE(deg_array);
    P4EST_FREE(eta_array);
    p8est_connectivity_destroy(conn_vtk);
    p8est_geometry_destroy(geom_vtk);
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
  mpi_assert((P4EST_DIM) == 3);
  
  int world_rank,world_size;
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &world_rank);
  sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &world_size);

  problem_input_t input = problem_input(input_file);
  
  int level;
  int endlevel = input.endlevel;
  int degree = input.degree;         
  double gamma_h = input.gamma_h;
  double gamma_p = input.gamma_p;
  double domain_size = input.domain_size;
  int percentile = input.percentile;
  int inflation_percentile = 25;
  
  R = .5/domain_size;
  double alpha_crit = sqrt(5);
  /* double rhoc = .032/(R*R); */
  double rhoc = (3./(2.*pi))*(1.0/(R*R))*((double)(5*5*5*5*5)/(double)(6*6*6*6*6*6));  
  
  rho0 = input.rho0_div_rhoc*rhoc;
  cx = .5;
  cy = .5;
  cz = .5;
  C0 = pow(1./(2.*pi*rho0/3.),.25);
  alpha = 386.266;

  ip_flux_params_t ip_flux_params;
  ip_flux_params.ip_flux_penalty_prefactor = input.ip_flux_penalty;
  ip_flux_params.ip_flux_penalty_calculate_fcn = sipg_flux_vector_calc_penalty_maxp2_over_minh;

  penalty_calc_t bi_u_penalty_fcn = bi_u_prefactor_conforming_maxp_minh;
  penalty_calc_t bi_u_dirichlet_penalty_fcn = bi_u_prefactor_conforming_maxp_minh;
  penalty_calc_t bi_gradu_penalty_fcn = bi_gradu_prefactor_maxp_minh;
  
  mpi_assert( !util_bisection(solve_for_alpha, alpha_crit, 1000*alpha_crit, DBL_EPSILON, 100000, &alpha) );

  double u_alpha_at_R = sqrt(alpha*R)/sqrt(R*R + alpha*R*alpha*R);
  beta = R*(C0*u_alpha_at_R - 1.);
  
  mpi_assert(
             (C0*u_alpha(R + .5,.5,.5) == 1. + beta/R)
             &&
             (C0*u_alpha(.5,R + .5,.5) == 1. + beta/R)
             &&
             (C0*u_alpha(.5,.5,R + .5) == 1. + beta/R)
            );
  
  if (world_rank == 0){
    printf("\n");
    printf("Amr levels= %d\n", endlevel);
    printf("Initial Degree = %d\n", degree);
    printf("Ip Flux Penalty = %f\n", ip_flux_params.ip_flux_penalty_prefactor);
    printf("amr percentile = %d\n", percentile);
    printf("smooth_pred_gamma_h = %f\n", gamma_h);
    printf("smooth_pred_gamma_p = %f\n", gamma_p);
    printf("Domain size = %f\n", domain_size);
    printf("rho0/rhoc = %.25f\n", rho0/rhoc);
    printf("alpha = %.25f\n", alpha);
    printf("beta = %.25f\n", beta);
    printf("C0 = %.25f\n", C0);
    printf("R = %.25f\n", R);
  }
  
  /* if (!load_from_checkpoint){ */
  p4est_partition_ext(p4est, 0, NULL);
  p4est_balance_ext(p4est, P4EST_CONNECT_FACE, NULL, NULL);
  /* } */

  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  /* create space for storing the ghost data */
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                            ghost->ghosts.elem_count);

  /* if (!load_from_checkpoint){ */
    p4est_reset_data(p4est, sizeof(d4est_element_data_t), NULL, NULL);


  d4est_mesh_geometry_storage_t* geometric_factors = geometric_factors_init(p4est);
    
    d4est_element_data_init_new(p4est,
                                 geometric_factors,
                                 d4est_ops,
                                 d4est_geom,
                                 problem_set_degrees,
                                 (void*)&input);

  /* } */
  
    int local_nodes = d4est_element_data_get_local_nodes(p4est);
  double* Au = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u_analytic = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u_error = P4EST_ALLOC_ZERO(double, local_nodes);

  double local_eta2 = -1.;

  problem_data_t prob_vecs;
  prob_vecs.Au = Au;
  prob_vecs.u = u;
  prob_vecs.u0 = u;
  prob_vecs.local_nodes = local_nodes;
  prob_vecs.curved_scalar_flux_fcn_data = curved_gauss_primal_sipg_flux_dirichlet_fetch_fcns
                                           (boundary_fcn,
                                            &ip_flux_params);

  d4est_linalg_fill_vec(u, 1., prob_vecs.local_nodes);

  curved_weakeqn_ptrs_t prob_fcns;
  problem_ctx_t prob_ctx;
  mpi_assert(input.deg_offset_for_nonlinear_quad > -1);
  prob_ctx.deg_offset_for_nonlinear_quad = input.deg_offset_for_nonlinear_quad;
  prob_ctx.ip_flux_params = &ip_flux_params;

  prob_fcns.apply_lhs = apply_jac;
  prob_fcns.build_residual = build_residual;
  prob_vecs.user = (void*)&prob_ctx;

  curved_smooth_pred_marker_t amr_marker;
  amr_marker.user = (void*)&input;
  amr_marker.mark_element_fcn = amr_mark_element;
  amr_marker.set_element_gamma_fcn = amr_set_element_gamma;
  amr_marker.name = "cds_marker";
  amr_marker.user = (void*)&input;

  hp_amr_scheme_t* scheme =
    hp_amr_curved_smooth_pred_init
    (
     p4est,
     input.degmax,
     amr_marker
    );
    
  d4est_mesh_init_field(p4est, u_analytic, analytical_solution_fcn, d4est_ops, d4est_geom->p4est_geom);    
  d4est_linalg_vec_axpyeqz(-1., u, u_analytic, u_error, local_nodes);

    
    /* dg norm should always have the boundary fcn set to zero */
    double local_dg_norm_sqr = d4est_element_data_compute_dg_norm_sqr
                               (
                                p4est,
                                u_error,
                                local_nodes,                                
                                &ip_flux_params,
                                d4est_geom,
                                ghost,
                                ghost_data,
                                d4est_ops
                               );
    
    double local_l2_norm_sqr =  d4est_element_data_compute_l2_norm_sqr
                                (
                                 p4est,
                                 u_error,
                                 local_nodes,
                                 d4est_ops,
                                 DO_NOT_STORE_LOCALLY
                                );

    double local_nodes_dbl = (double)local_nodes;
    double local_reduce [3];
    local_reduce[0] = local_nodes_dbl;
    local_reduce[1] = local_l2_norm_sqr;
    local_reduce[2] = local_dg_norm_sqr;

    double global_reduce [3];

    sc_reduce
      (
       &local_reduce[0],
       &global_reduce[0],
       3,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       0,
       sc_MPI_COMM_WORLD
      );

    double global_nodes_dbl = global_reduce[0];
    double global_l2_norm_sqr = global_reduce[1];
    double global_dg_norm_sqr = global_reduce[2];
    
    if (world_rank == 0){
      printf
        (
         "[HP_AMR]: Level Elements Nodes Eta2 L2-Norm DG-Norm Time\n"
        );
      printf
        (
         "\n\n[HP_AMR]: %d %d %d %.25f %.25f %.25f %f \n\n",
         0,
         (int)p4est->global_num_quadrants,
         (int)global_nodes_dbl,
         -1.,
         sqrt(global_l2_norm_sqr),
         sqrt(global_dg_norm_sqr),
         0.
        );
    }

    double* dof_data_for_fit = P4EST_ALLOC(double, endlevel);
    double* dgerr_data_for_fit = P4EST_ALLOC(double, endlevel);
    
  for (level = 0; level < endlevel; ++level){

    
    curved_bi_estimator_compute
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       bi_u_penalty_fcn,
       bi_u_dirichlet_penalty_fcn,
       bi_gradu_penalty_fcn,
       boundary_fcn,
       ip_flux_params.ip_flux_penalty_prefactor,
       ghost,
       ghost_data,
       d4est_ops,
       d4est_geom
      );
    
    estimator_stats_t* stats = P4EST_ALLOC(estimator_stats_t, 1);
    estimator_stats_compute(p4est, stats,1);

    if(world_rank == 0)
      estimator_stats_print(stats);

    local_eta2 = stats->total;
        
    d4est_mesh_init_field(p4est, u_analytic, analytical_solution_fcn, d4est_ops, d4est_geom->p4est_geom);    
    d4est_linalg_vec_axpyeqz(-1., u, u_analytic, u_error, local_nodes);

    int with_eta = 1;
    problem_save_to_vtk(p4est, d4est_ops,u, u_analytic, u_error, level, with_eta);
    

    printf("p4est->local_num_quadrants*p4est->mpisize < input.amr_inflation_size -> %d < %d = %d\n", p4est->local_num_quadrants*p4est->mpisize,input.amr_inflation_size, p4est->local_num_quadrants*p4est->mpisize < input.amr_inflation_size);

    int curved = 1;
    hp_amr(p4est,
           d4est_ops,
           &u,
           &stats,
           scheme,
           curved
          );        

    P4EST_FREE(stats);
    p4est_ghost_destroy(ghost);
    P4EST_FREE(ghost_data);

    ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
    ghost_data = P4EST_ALLOC(d4est_element_data_t, ghost->ghosts.elem_count);
    
    d4est_element_data_init_new(p4est,
                                 geometric_factors,
                                 d4est_ops,
                                 d4est_geom,
                                 problem_keep_degrees,
                                 (void*)&input);
    
    local_nodes = d4est_element_data_get_local_nodes(p4est);

    Au = P4EST_REALLOC(Au, double, local_nodes);
    u_analytic = P4EST_REALLOC(u_analytic, double, local_nodes);
    u_error = P4EST_REALLOC(u_error, double, local_nodes);
  
    prob_vecs.Au = Au;
    prob_vecs.u = u;
    prob_vecs.u0 = u;
    prob_vecs.local_nodes = local_nodes;
    
    clock_t begin = 0;
    clock_t end = -1;

    if (world_rank == 0){
      begin = clock();
    }  

    /* krylov_petsc_params_t params; */
    /* params.user_defined_pc = 0; */
    /* params.ksp_monitor = 0; */
    /* params.ksp_type = input.ksp_type; */
    
  newton_petsc_solve
    (
     p4est,
     &prob_vecs,
     (void*)&prob_fcns,
     &ghost,
     (void**)&ghost_data,
     d4est_ops,
     d4est_geom,
     input_file,
     NULL
    );
    
  d4est_mesh_init_field(p4est, u_analytic, analytical_solution_fcn, d4est_ops, d4est_geom->p4est_geom);    
    d4est_linalg_vec_axpyeqz(-1., u, u_analytic, u_error, local_nodes);

    
    /* dg norm should always have the boundary fcn set to zero */
    double local_dg_norm_sqr = d4est_element_data_compute_dg_norm_sqr
                               (
                                p4est,
                                u_error,
                                local_nodes,
                                &ip_flux_params,
                                d4est_geom,
                                ghost,
                                ghost_data,
                                d4est_ops
                               );
    
    double local_l2_norm_sqr =  d4est_element_data_compute_l2_norm_sqr
                                (
                                 p4est,
                                 u_error,
                                 local_nodes,
                                 d4est_ops,
                                 DO_NOT_STORE_LOCALLY
                                );

    double local_nodes_dbl = (double)local_nodes;
    /* int stride = 0; */
    
    double local_reduce [4];
    local_reduce[0] = local_nodes_dbl;
    local_reduce[1] = local_l2_norm_sqr;
    local_reduce[2] = local_dg_norm_sqr;
    local_reduce[3] = local_eta2;

    double global_reduce [4];

    sc_reduce
      (
       &local_reduce[0],
       &global_reduce[0],
       4,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       0,
       sc_MPI_COMM_WORLD
      );

    double global_nodes_dbl = global_reduce[0];
    double global_l2_norm_sqr = global_reduce[1];
    double global_dg_norm_sqr = global_reduce[2];
    double global_eta2 = global_reduce[3];
    
    if (world_rank == 0){
      end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
      printf
        (
         "\n\n[HP_AMR]: %d %d %d %.25f %.25f %.25f %f \n\n",
         level,
         (int)p4est->global_num_quadrants,
         (int)global_nodes_dbl,
         sqrt(global_eta2),
         sqrt(global_l2_norm_sqr),
         sqrt(global_dg_norm_sqr),
         time_spent
        );

      dgerr_data_for_fit[level] = log(sqrt(global_dg_norm_sqr));
      dof_data_for_fit[level] = pow(global_nodes_dbl, 1./(2.*(P4EST_DIM)-1.));

      if (level > 0){
        double slope;
        double intercept;
        int num_of_hpamr_levels = level + 1;
        util_linear_regression
          (
           dgerr_data_for_fit,
           dof_data_for_fit,
           &slope,
           &intercept,
           num_of_hpamr_levels
          );
        printf("[HP_AMR_FIT](1): ||err||DG = C1*exp(-C2*DOF^(1/%d))\n",2*(P4EST_DIM)-1);
        printf("[HP_AMR_FIT](2): LEV SLOPE DG_ERR\n");
        printf("[HP_AMR_FIT](3): %d %.25f %.25f, \n\n", level, slope, sqrt(global_dg_norm_sqr));
      }      
    }


    /* int checkpoint_period = 3; */
    
    /* if (level % checkpoint_period == 0){ */
    /*   char* p4est_filename = NULL; */
    /*   D4EST_ASPRINTF(p4est_filename, "%s_%s_%d", "checkpoint","level", level); */
    /*   element_data_copy_from_vec_to_storage(p4est, u); */
    /*   p4est_save_ext(p4est_filename, p4est, 1, 1); */
    /*   free(p4est_filename); */
    /* } */

   

  }



  d4est_mesh_init_field(p4est, u_analytic, analytical_solution_fcn, d4est_ops, d4est_geom->p4est_geom);    
  d4est_linalg_vec_axpyeqz(-1., u, u_analytic, u_error, local_nodes);

  int with_eta = 0;
  problem_save_to_vtk(p4est, d4est_ops, u, u_analytic, u_error, level, with_eta);
  
  hp_amr_curved_smooth_pred_destroy(scheme);
  geometric_factors_destroy(geometric_factors);
  
  /* if (ghost) { */
  p4est_ghost_destroy (ghost);
  P4EST_FREE (ghost_data);
  ghost = NULL;
  ghost_data = NULL;
  /* } */

  P4EST_FREE(dof_data_for_fit);
  P4EST_FREE(dgerr_data_for_fit);  
  P4EST_FREE(Au);
  P4EST_FREE(u_analytic);
  P4EST_FREE(u);
  P4EST_FREE(u_error);
}
