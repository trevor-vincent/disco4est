#define NDEBUG
#define _GNU_SOURCE

#include <sc_reduce.h>
#include <pXest.h>
#include <util.h>
#include <linalg.h>
#include <element_data.h>
#include <sipg_flux_vector_fcns.h>
#include <sipg_flux_scalar_fcns.h>
#include <problem.h>
#include <problem_data.h>
#include <problem_weakeqn_ptrs.h>
#include <poisson_operator.h>
#include <hp_amr_smooth_pred.h>
#include <hp_amr_uniform.h>
#include <bi_estimator.h>
#include <estimator_stats.h>
#include <bi_estimator_flux_fcns.h>
#include <hp_amr.h>
#include <multigrid.h>
#include <hacked_p4est_vtk.h>
#include <dg_norm.h>
#include <ini.h>
#include <multigrid_smoother_cheby.h>
#include <multigrid_bottom_solver_cg.h>
#include <multigrid_logger_residual.h>
#include <multigrid_element_data_updater_nocurved.h>
#include "time.h"

static const double pi = 3.1415926535897932384626433832795;
static double c2x;// = .5;
static double c2y;// = .5;
static double c2z;// = .5;

static
double analytic_solution_fcn
(
 double x,
 double y
#if (P4EST_DIM)==3
 ,
 double z
#endif
)
{
  double xp = x - c2x;
  double yp = y - c2y;
#if (P4EST_DIM)==3
  double zp = z - c2z;
#endif
  double rp = sqrt(xp*xp
                   + yp*yp
#if (P4EST_DIM)==3
                   + zp*zp
#endif   
                  );

#if (P4EST_DIM)==2  
  return x*(1.-x)*y*(1.-y)*util_dbl_pow_int(rp, 3);
#else
  return x*(1.-x)*y*(1.-y)*z*(1.-z)*util_dbl_pow_int(rp, 3);
#endif
}

static
double boundary_fcn
(
 double x,
 double y
#if (P4EST_DIM)==3
 ,
 double z
#endif
)

{
  return 0.;
}


/**
 * returns - \nabla \cdot (1 + .5\sin(|\nabla|^2)) \nabla u
 *
 * @param x
 * @param (P4EST_DIM)
 * @param endif
 *
 * @return
 */
static
double f_fcn
(
 double x,
 double y
#if (P4EST_DIM)==3
 ,
 double z
#endif
)
{
  #if (P4EST_DIM)==2
  if (x == c2x && y == c2y){
    return 0.;
  }

  else
    return ((util_dbl_pow_int(c2x,2) + util_dbl_pow_int(c2y,2) - 2*c2x*x + util_dbl_pow_int(x,2) - 2*c2y*y + util_dbl_pow_int(y,2))*(-2*util_dbl_pow_int(c2x,2)*x - 6*c2y*x - 2*util_dbl_pow_int(c2y,2)*x + 4*c2x*util_dbl_pow_int(x,2) + 2*util_dbl_pow_int(c2x,2)*util_dbl_pow_int(x,2) + 6*c2y*util_dbl_pow_int(x,2) +
       2*util_dbl_pow_int(c2y,2)*util_dbl_pow_int(x,2) - 2*util_dbl_pow_int(x,3) - 4*c2x*util_dbl_pow_int(x,3) + 2*util_dbl_pow_int(x,4) - 6*c2x*y - 2*util_dbl_pow_int(c2x,2)*y - 2*util_dbl_pow_int(c2y,2)*y + 21*x*y + 16*c2x*x*y + 16*c2y*x*y - 29*util_dbl_pow_int(x,2)*y -
       16*c2y*util_dbl_pow_int(x,2)*y + 6*c2x*util_dbl_pow_int(y,2) + 2*util_dbl_pow_int(c2x,2)*util_dbl_pow_int(y,2) + 4*c2y*util_dbl_pow_int(y,2) + 2*util_dbl_pow_int(c2y,2)*util_dbl_pow_int(y,2) - 29*x*util_dbl_pow_int(y,2) - 16*c2x*x*util_dbl_pow_int(y,2) + 37*util_dbl_pow_int(x,2)*util_dbl_pow_int(y,2) -
                                                                                                                          2*util_dbl_pow_int(y,3) - 4*c2y*util_dbl_pow_int(y,3) + 2*util_dbl_pow_int(y,4)))/sqrt(util_dbl_pow_int(-c2x + x,2) + util_dbl_pow_int(-c2y + y,2));

#else
  if (x == c2x && y == c2y && z == c2z){
    return 0.;
  }
  else{
    return (-2*(1 - x)*x*(1 - y)*y*util_dbl_pow_int(util_dbl_pow_int(c2x - x,2) + util_dbl_pow_int(c2y - y,2) + util_dbl_pow_int(c2z - z,2),2) + 
     9*(1 - x)*x*(1 - y)*y*(util_dbl_pow_int(c2x - x,2) + util_dbl_pow_int(c2y - y,2) + util_dbl_pow_int(c2z - z,2))*(1 - z)*
     z + 6*(1 - x)*(-c2x + x)*(1 - y)*y*
     (util_dbl_pow_int(c2x - x,2) + util_dbl_pow_int(c2y - y,2) + util_dbl_pow_int(c2z - z,2))*(1 - z)*z - 
     6*x*(-c2x + x)*(1 - y)*y*(util_dbl_pow_int(c2x - x,2) + util_dbl_pow_int(c2y - y,2) + util_dbl_pow_int(c2z - z,2))*
     (1 - z)*z + 6*(1 - x)*x*(1 - y)*(-c2y + y)*
     (util_dbl_pow_int(c2x - x,2) + util_dbl_pow_int(c2y - y,2) + util_dbl_pow_int(c2z - z,2))*(1 - z)*z - 
     6*(1 - x)*x*y*(-c2y + y)*(util_dbl_pow_int(c2x - x,2) + util_dbl_pow_int(c2y - y,2) + util_dbl_pow_int(c2z - z,2))*
     (1 - z)*z - 2*(1 - x)*x*util_dbl_pow_int(util_dbl_pow_int(c2x - x,2) + util_dbl_pow_int(c2y - y,2) + util_dbl_pow_int(c2z - z,2),
                                              2)*(1 - z)*z - 2*(1 - y)*y*util_dbl_pow_int(util_dbl_pow_int(c2x - x,2) + util_dbl_pow_int(c2y - y,2) + 
                                                                                          util_dbl_pow_int(c2z - z,2),2)*(1 - z)*z - 
     3*util_dbl_pow_int(c2x - x,2)*(-1 + x)*x*(-1 + y)*y*(-1 + z)*z - 
     3*(-1 + x)*x*util_dbl_pow_int(c2y - y,2)*(-1 + y)*y*(-1 + z)*z - 
     3*(-1 + x)*x*(-1 + y)*y*util_dbl_pow_int(c2z - z,2)*(-1 + z)*z + 
     6*(1 - x)*x*(1 - y)*y*(util_dbl_pow_int(c2x - x,2) + util_dbl_pow_int(c2y - y,2) + util_dbl_pow_int(c2z - z,2))*(1 - z)*
     (-c2z + z) - 6*(1 - x)*x*(1 - y)*y*
     (util_dbl_pow_int(c2x - x,2) + util_dbl_pow_int(c2y - y,2) + util_dbl_pow_int(c2z - z,2))*z*(-c2z + z))/
      sqrt(util_dbl_pow_int(c2x - x,2) + util_dbl_pow_int(c2y - y,2) + util_dbl_pow_int(c2z - z,2));
  }
#endif
}

static
void problem_build_rhs
(
 p4est_t* p4est,
 problem_data_t* prob_vecs,
 weakeqn_ptrs_t* prob_fcns,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_operators_t* d4est_ops
)
{
  prob_vecs->scalar_flux_fcn_data.bndry_fcn = boundary_fcn;
  prob_vecs->vector_flux_fcn_data.bndry_fcn = boundary_fcn;

  double* f = P4EST_ALLOC(double, prob_vecs->local_nodes);
  element_data_init_node_vec
    (
     p4est,
     f,
     f_fcn,
     d4est_ops
    );
  
  element_data_apply_Mij_on_vec
    (
     p4est,
     f,
     prob_vecs->rhs,
     d4est_ops
    );
  
  linalg_vec_scale(-1., prob_vecs->rhs, prob_vecs->local_nodes);

  /* DEBUG_PRINT_ARR_DBL(prob_vecs->rhs, prob_vecs->local_nodes); */
  
  int local_nodes = prob_vecs->local_nodes;
  double* u_eq_0 = P4EST_ALLOC_ZERO(double, local_nodes);
  double* tmp = prob_vecs->u;
  
  prob_vecs->u = u_eq_0;
  poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops, NULL);
  linalg_vec_axpy(-1., prob_vecs->Au, prob_vecs->rhs, local_nodes);

  /* DEBUG_PRINT_ARR_DBL(prob_vecs->Au, prob_vecs->local_nodes); */
  
  prob_vecs->u = tmp;
  P4EST_FREE(u_eq_0);
  P4EST_FREE(f);

  prob_vecs->scalar_flux_fcn_data.bndry_fcn = zero_fcn;
  prob_vecs->vector_flux_fcn_data.bndry_fcn = zero_fcn;
}

static
void
build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
  poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops, NULL);
  linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);
}
 
p4est_connectivity_t*
problem_build_conn()
{
#if (P4EST_DIM)==3
  return p8est_connectivity_new_unitcube();
#elif (P4EST_DIM)==2
  return p4est_connectivity_new_unitsquare();
#else
  mpi_abort("[D4EST_ERROR]: Dim not supported");
#endif

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

p4est_t*
problem_load_p4est_from_checkpoint
(
 const char* filename,
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t** conn
){
  int autopartition = 1;
  int load_data = 1;
  int broadcasthead = 0;
  
  return p4est_load_ext (filename,
                mpicomm,
                sizeof(element_data_t),
                load_data,
                autopartition,
                broadcasthead,
                NULL,
                conn);
}

typedef struct {

  int endlevel;
  int degree;
  double gamma_h;
  double gamma_p;
  double sigma;
  double c2x;
  double c2y;
  double c2z;
  double ip_flux_penalty;
  
  int count;
  
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
  if (util_match_couple(section,"amr",name,"amr_levels")) {
    mpi_assert(pconfig->endlevel == -1);
    pconfig->endlevel = atoi(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"amr",name, "initial_degree")) {
    mpi_assert(pconfig->degree == -1);
    pconfig->degree = atoi(value);
    pconfig->count += 1;
  } else if (util_match_couple(section,"amr",name,"sigma")) {
    mpi_assert(pconfig->sigma == -1);
    pconfig->sigma = atof(value);
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
  } 
  else if (util_match_couple(section,"problem",name,"c2x")) {
    mpi_assert(pconfig->c2x == -1);
    pconfig->c2x = atof(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"c2y")) {
    mpi_assert(pconfig->c2y == -1);
    pconfig->c2y = atof(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"c2z")) {
    mpi_assert(pconfig->c2z == -1);
    pconfig->c2z = atof(value);
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
  int num_of_options = 9;
  
  problem_input_t input;
  input.degree = -1;
  input.endlevel = -1;
  input.gamma_h = -1;
  input.gamma_p = -1;
  input.ip_flux_penalty = -1;
  input.sigma = -1;
  input.c2x = -1;
  input.c2y = -1;
  input.c2z = -1;

  input.count = 0;
  
  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  mpi_assert(input.count == num_of_options);
  return input;
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
  mpi_assert((P4EST_DIM) == 3 || (P4EST_DIM) == 2);
  
  int world_rank,world_size;
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &world_rank);
  sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &world_size);

  problem_input_t input = problem_input(input_file);
  
  int level;
  int endlevel = input.endlevel;
  int degree = input.degree;         
  double gamma_h = input.gamma_h;
  double gamma_p = input.gamma_p;
  double sigma = input.sigma;
  c2x = input.c2x;
  c2y = input.c2y;
  c2z = input.c2z;

  
  ip_flux_params_t ip_flux_params;
  ip_flux_params.ip_flux_penalty_prefactor = input.ip_flux_penalty;
  ip_flux_params.ip_flux_penalty_calculate_fcn = sipg_flux_vector_calc_penalty_maxp2_over_minh;

  /* penalty_calc_t bi_u_penalty_fcn = bi_u_prefactor_conforming_maxp_minh; */
  /* penalty_calc_t bi_u_dirichlet_penalty_fcn = bi_u_prefactor_conforming_maxp_minh; */
  /* penalty_calc_t bi_gradu_penalty_fcn = bi_gradu_prefactor_maxp_minh; */

  penalty_calc_t bi_u_penalty_fcn = houston_u_prefactor_maxp_minh;
  penalty_calc_t bi_u_dirichlet_penalty_fcn = houston_u_dirichlet_prefactor_maxp_minh;
  penalty_calc_t bi_gradu_penalty_fcn = houston_gradu_prefactor_maxp_minh;
  
    p4est_partition_ext(p4est, 0, NULL);
    p4est_balance_ext(p4est, P4EST_CONNECT_FACE, NULL, NULL);

  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  /* create space for storing the ghost data */
  element_data_t* ghost_data = P4EST_ALLOC (element_data_t,
                                            ghost->ghosts.elem_count);

    p4est_reset_data(p4est, sizeof(element_data_t), NULL, NULL);
    element_data_init(p4est, degree);
  
  int local_nodes = element_data_get_local_nodes(p4est);
  double* Au = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u_analytic = P4EST_ALLOC_ZERO(double, local_nodes);
  double* rhs = P4EST_ALLOC_ZERO(double, local_nodes);

  double local_eta2 = -1.;

  problem_data_t prob_vecs;
  prob_vecs.Au = Au;
  prob_vecs.u = u;
  prob_vecs.rhs = rhs;
  prob_vecs.local_nodes = local_nodes;

  prob_vecs.vector_flux_fcn_data = sipg_flux_vector_dirichlet_fetch_fcns
                                          (
                                           boundary_fcn,
                                           &ip_flux_params
                                          );
  
  prob_vecs.scalar_flux_fcn_data = sipg_flux_scalar_dirichlet_fetch_fcns(boundary_fcn);
  linalg_fill_vec(u, 100., prob_vecs.local_nodes);
  
  weakeqn_ptrs_t prob_fcns;

  prob_fcns.apply_lhs = poisson_apply_aij;
  prob_fcns.build_residual = build_residual;


  problem_build_rhs
    (
     p4est,
     &prob_vecs,
     &prob_fcns,
     ghost,
     ghost_data,
     d4est_ops
    );
  
 
  hp_amr_scheme_t* scheme =
    hp_amr_smooth_pred_init
    (
     p4est,
     gamma_h,
     gamma_p,
     1.,
     7,
     hp_amr_smooth_pred_get_sigaverage_marker(&sigma)
    );


  hp_amr_scheme_t* scheme_uni = hp_amr_uniform();

  
  element_data_init_node_vec(p4est, u_analytic, analytic_solution_fcn, d4est_ops);    
  linalg_vec_axpy(-1., u, u_analytic, local_nodes);

    
    /* dg norm should always have the boundary fcn set to zero */
    double local_dg_norm_sqr = element_data_compute_DG_norm_sqr
                               (
                                p4est,
                                u_analytic,
                                zero_fcn,
                                &ip_flux_params,
                                d4est_ops,
                                ghost,
                                ghost_data
                               );
    
    double local_l2_norm_sqr =  element_data_compute_l2_norm_sqr_no_local
                                (
                                 p4est,
                                 u_analytic,
                                 d4est_ops
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

    bi_estimator_compute
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
       d4est_ops
      );
    
    estimator_stats_t* stats = P4EST_ALLOC(estimator_stats_t, 1);
    estimator_stats_compute(p4est, stats,0);

    /* if(world_rank == 0) */
      /* estimator_stats_print(&stats, 0); */

    local_eta2 = stats->total;
        
    element_data_init_node_vec(p4est, u_analytic, analytic_solution_fcn, d4est_ops);    

    double* u_analytic_vertex = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN));
    double* u_error_vertex = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN));
    double* u_vertex = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN));
    double* eta2_vertex = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN));
    element_data_store_local_estimator_in_corner_array
      (
       p4est,
       eta2_vertex
      );
    
    element_data_store_nodal_vec_in_vertex_array
      (
       p4est,
       u_analytic,
       u_analytic_vertex
      );

    element_data_store_nodal_vec_in_vertex_array
      (
       p4est,
       u,
       u_vertex
      );
    
    linalg_vec_axpy(-1., u, u_analytic, local_nodes);

    element_data_store_nodal_vec_in_vertex_array
      (
       p4est,
       u_analytic,
       u_error_vertex
      );

    linalg_vec_fabs(u_error_vertex, (P4EST_CHILDREN)*p4est->local_num_quadrants);
    
    char sol_save_as [500];
    sprintf(sol_save_as, "%s_hp_amr_level_%d_sols", P4EST_STRING, level);

    hacked_p4est_vtk_write_all
      (p4est,
       NULL,
       0.99,
       0,   
       1,   
       1,   
       0,
       4,
       0,  
       sol_save_as,
       "eta2",
       eta2_vertex,
       "psi",
       u_vertex,
       "psi_analytic",
       u_analytic_vertex,
       "psi_error",
       u_error_vertex
      );

    P4EST_FREE(u_analytic_vertex);
    P4EST_FREE(u_error_vertex);
    P4EST_FREE(u_vertex);   
    P4EST_FREE(eta2_vertex);   

    if(world_rank == 0)
      estimator_stats_print(stats);
    
    int curved = 0;
    hp_amr(p4est,
           d4est_ops,
           &u,
           &stats,
           (level > 1) ? scheme : scheme_uni,
           curved
          );
    
    P4EST_FREE(stats);
    
    p4est_ghost_destroy(ghost);
    P4EST_FREE(ghost_data);

    ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
    ghost_data = P4EST_ALLOC(element_data_t, ghost->ghosts.elem_count);
    
    element_data_init(p4est, -1);
    local_nodes = element_data_get_local_nodes(p4est);

    Au = P4EST_REALLOC(Au, double, local_nodes);
    u_analytic = P4EST_REALLOC(u_analytic, double, local_nodes);
    rhs = P4EST_REALLOC(rhs, double, local_nodes);
  
    prob_vecs.Au = Au;
    prob_vecs.u = u;
    prob_vecs.rhs = rhs;
    prob_vecs.local_nodes = local_nodes;

    problem_build_rhs
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       ghost,
       ghost_data,
       d4est_ops
      );    
    
    clock_t begin = 0;
    clock_t end = -1;

    if (world_rank == 0){
      begin = clock();
    }  

    int min_level, max_level;

    multigrid_get_level_range(p4est, &min_level, &max_level);
    printf("[min_level, max_level] = [%d,%d]\n", min_level, max_level);

    /* need to do a reduce on min,max_level before supporting multiple proc */
    /* mpi_assert(proc_size == 1); */
    int num_of_levels = max_level + 1;
    
    /* int vcycle_iter = 3; */
    /*   double vcycle_rtol = 1e-9; */
    /*   double vcycle_atol = 0.; */
    /*   int smooth_iter = 15; */
    /*   int cg_eigs_iter = 10; */
    /*   /\* double max_eig_factor = 1.1; *\/ */
    /*   double max_eig_factor = 1.0; */
    /*   int max_eig_reuse = 1; */
    /*   double lmax_lmin_rat = 30.; */
    /*   int coarse_iter = 10000; */
    /*   double coarse_rtol = 1e-10; */
    /*   int save_vtk_snapshot = 0; */
    /*   int perform_checksum = 0; */
    /*   int cg_eigs_use_zero_vec_as_initial = 0; */

    multigrid_smoother_t* smoother = multigrid_smoother_cheby_d4est_init
                                   (
                                    p4est,
                                    num_of_levels,
                                    input_file
                                   );


    multigrid_bottom_solver_t* bottom_solver = multigrid_bottom_solver_cg_init
                                               (
                                                p4est,
                                                input_file
                                               );

    multigrid_logger_t* logger = multigrid_logger_residual_init
                                 (
                                 );

    multigrid_element_data_updater_t* updater = multigrid_element_data_updater_nocurved_init(
                                                                                             &ghost,
                                                                                             &ghost_data
    );

    
    multigrid_data_t* mg_data = multigrid_data_init(p4est,
                                                    d4est_ops,
                                                    num_of_levels,
                                                    smoother,
                                                    bottom_solver,
                                                    logger,
                                                    NULL,
                                                    updater,
                                                    input_file
                                                   );


      element_data_init_node_vec(p4est, u_analytic, analytic_solution_fcn, d4est_ops);

      multigrid_solve
        (
         p4est,
         &prob_vecs,
         &prob_fcns,
         mg_data
        );


      multigrid_smoother_cheby_d4est_destroy(smoother);
      multigrid_bottom_solver_cg_destroy(bottom_solver);
      multigrid_logger_residual_destroy(logger);
      multigrid_element_data_updater_nocurved_destroy(updater);
      
  
      multigrid_data_destroy(mg_data);

      

    linalg_vec_axpy(-1., u, u_analytic, local_nodes);
    
    /* dg norm should always have the boundary fcn set to zero */
    double local_dg_norm_sqr = element_data_compute_DG_norm_sqr
                               (
                                p4est,
                                u_analytic,
                                zero_fcn,
                                &ip_flux_params,
                                d4est_ops,
                                ghost,
                                ghost_data
                               );
    
    double local_l2_norm_sqr =  element_data_compute_l2_norm_sqr_no_local
                                (
                                 p4est,
                                 u_analytic,
                                 d4est_ops
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
  }

  double* u_analytic_vertex = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN));
  double* u_error_vertex = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN));
  double* u_vertex = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN));
    
  element_data_store_nodal_vec_in_vertex_array
    (
     p4est,
     u_analytic,
     u_analytic_vertex
    );

  element_data_store_nodal_vec_in_vertex_array
    (
     p4est,
     u,
     u_vertex
    );
    
  linalg_vec_axpy(-1., u, u_analytic, local_nodes);

  element_data_store_nodal_vec_in_vertex_array
    (
     p4est,
     u_analytic,
     u_error_vertex
    );

  linalg_vec_fabs(u_error_vertex, (P4EST_CHILDREN)*p4est->local_num_quadrants);
  
  char sol_save_as [500];
  sprintf(sol_save_as, "%s_hp_amr_level_%d_sols_noeta2", P4EST_STRING, level);

  hacked_p4est_vtk_write_all
    (p4est,
     NULL,
     0.99,
     0,   
     1,   
     1,   
     0,
     3,
     0,  
     sol_save_as,
     "u",
     u_vertex,
     "u_analytic",
     u_analytic_vertex,
     "psi_error",
     u_error_vertex
    );


  
  P4EST_FREE(u_analytic_vertex);
  P4EST_FREE(u_error_vertex);
  P4EST_FREE(u_vertex);   

  hp_amr_smooth_pred_destroy(scheme);
  hp_amr_uniform_destroy(scheme_uni);
  
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
  P4EST_FREE(rhs);
}
