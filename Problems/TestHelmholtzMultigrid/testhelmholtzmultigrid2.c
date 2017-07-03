#include <sc_reduce.h>
#include <pXest.h>
#include <util.h>
#include <d4est_linalg.h>
#include <element_data.h>
#include <sipg_flux_vector_fcns.h>
#include <sipg_flux_scalar_fcns.h>
#include <problem.h>
#include <problem_data.h>
#include <d4est_elliptic_eqns.h>
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
#include <krylov_petsc.h>
#include <krylov_pc_multigrid.h>
#include <multigrid_matrix_operator.h>
#include "time.h"

double pi = 3.1415926535897932384626433832795;

static int
random_h_refine(p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * quadrant)
{
  /* d4est_element_data_t* data = (d4est_element_data_t*) quadrant->p.user_data; */
  return rand()%2;
}


static
double helmholtz_fcn
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  return 2*x*y*z;
}


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
#if (P4EST_DIM)==3
 return sin(pi*x)*sin(pi*y)*sin(pi*z);
#else
 return sin(pi*x)*sin(pi*y);
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
  return zero_fcn(x,y
#if (P4EST_DIM)==3
    ,z
#endif
    );
}

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
#if (P4EST_DIM)==3
  double u = sin(pi*x)*sin(pi*y)*sin(pi*z);
  return 3*pi*pi*u + helmholtz_fcn(x,y,z,0,NULL)*u;
#else
  double u = sin(pi*x)*sin(pi*y);
  return 2*pi*pi*u + helmholtz_fcn(x,y,z,0,NULL)*u;
#endif
}
 
static int
refine_function
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t *quadrant
)
{
  return 1;
}

typedef struct {

  double ip_flux_penalty;
  int deg;
  int endlevel;
  int degree;
  double gamma_h;
  double gamma_p;
  double sigma;
  
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
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}



static
void apply_helmholtz_matrix_2
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_elliptic_problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops
)
{  
  poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops);

  /* printf("Quadrants = %d\n", p4est->local_num_quadrants); */
  
  double* M_helmf_u = P4EST_ALLOC(double, prob_vecs->local_nodes);
  int matrix_stride = 0;

  /* int local_matrix_nodes = element_data_get_local_matrix_nodes(p4est); */

  /* DEBUG_PRINT_ARR_DBL(((multigrid_matrix_op_t*)prob_vecs->user)->matrix_at0, local_matrix_nodes); */
  /* double debug_sum = 0.; */
  
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        element_data_t* ed = quad->p.user_data;

        int deg_gauss = ed->deg;
        int volume_nodes_gauss = d4est_lgl_get_nodes((P4EST_DIM), deg_gauss);
        double* r_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 0);
        double* s_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 1);
        double* t_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 2);
        double *jac_gauss = P4EST_ALLOC(double, volume_nodes_gauss);

        double* x_GL = P4EST_ALLOC(double, volume_nodes_gauss);
        double* y_GL = P4EST_ALLOC(double, volume_nodes_gauss);
        double* z_GL = P4EST_ALLOC(double, volume_nodes_gauss);

        d4est_linalg_fill_vec(jac_gauss, ed->jacobian, volume_nodes_gauss);  
        d4est_reference_rtox_array(r_GL, ed->xyz_corner[0], ed->h, x_GL, volume_nodes_gauss);
        d4est_reference_rtox_array(s_GL, ed->xyz_corner[1], ed->h, y_GL, volume_nodes_gauss);
        d4est_reference_rtox_array(t_GL, ed->xyz_corner[2], ed->h, z_GL, volume_nodes_gauss);
        double* xyz_gauss [3] = {x_GL, y_GL, z_GL};

        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        /* double* debug_ones = P4EST_ALLOC(double, volume_nodes); */
        /* d4est_linalg_fill_vec(debug_ones, 1., volume_nodes); */
        
        d4est_linalg_matvec_plus_vec(1.,&((multigrid_matrix_op_t*)prob_vecs->user)->matrix[matrix_stride], &prob_vecs->u[ed->stride], 0., &M_helmf_u[ed->stride], volume_nodes, volume_nodes);


        /* double* matrix = &((multigrid_matrix_op_t*)prob_vecs->user)->matrix[matrix_stride]; */
        /* DEBUG_PRINT_ARR_DBL(matrix, volume_nodes*volume_nodes); */
        /* debug_sum += d4est_linalg_vec1_trans_mat_vec2(debug_ones, matrix, debug_ones, volume_nodes); */
        
        
        P4EST_FREE(z_GL);
        P4EST_FREE(y_GL);
        P4EST_FREE(x_GL);
        P4EST_FREE(jac_gauss);
        /* P4EST_FREE(debug_ones); */

        matrix_stride += volume_nodes*volume_nodes;
      }
    }

  /* printf("debug_sum = %.25f\n", debug_sum); */
  d4est_linalg_vec_axpy(1.0, M_helmf_u, prob_vecs->Au, prob_vecs->local_nodes);
  P4EST_FREE(M_helmf_u);
}



static
void apply_helmholtz_matrix
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_elliptic_problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops
)
{    
  poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops);

  printf("Quadrants = %d\n", p4est->local_num_quadrants);
  
  double* M_helmf_u = P4EST_ALLOC(double, prob_vecs->local_nodes);
 
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        element_data_t* ed = quad->p.user_data;

        int deg_gauss = ed->deg;
        int volume_nodes_gauss = d4est_lgl_get_nodes((P4EST_DIM), deg_gauss);
        double* r_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 0);
        double* s_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 1);
        double* t_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 2);
        double *jac_gauss = P4EST_ALLOC(double, volume_nodes_gauss);

        double* x_GL = P4EST_ALLOC(double, volume_nodes_gauss);
        double* y_GL = P4EST_ALLOC(double, volume_nodes_gauss);
        double* z_GL = P4EST_ALLOC(double, volume_nodes_gauss);

        d4est_linalg_fill_vec(jac_gauss, ed->jacobian, volume_nodes_gauss);  
        d4est_reference_rtox_array(r_GL, ed->xyz_corner[0], ed->h, x_GL, volume_nodes_gauss);
        d4est_reference_rtox_array(s_GL, ed->xyz_corner[1], ed->h, y_GL, volume_nodes_gauss);
        d4est_reference_rtox_array(t_GL, ed->xyz_corner[2], ed->h, z_GL, volume_nodes_gauss);
        double* xyz_gauss [3] = {x_GL, y_GL, z_GL};

        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        double* matrix = P4EST_ALLOC(double, volume_nodes*volume_nodes);

        d4est_operators_form_fofufofvlilj_matrix_gaussnodes
          (
           d4est_ops,
           NULL,
           NULL,
           ed->deg,
           xyz_gauss,
           jac_gauss,
           ed->deg,
           (P4EST_DIM),
           matrix,
           helmholtz_fcn,
           NULL,
           NULL,
           NULL
          );

        DEBUG_PRINT_ARR_DBL(matrix, volume_nodes*volume_nodes);
        
        d4est_linalg_matvec_plus_vec(1., matrix, &prob_vecs->u[ed->stride], 0., &M_helmf_u[ed->stride], volume_nodes, volume_nodes);

        
        P4EST_FREE(z_GL);
        P4EST_FREE(y_GL);
        P4EST_FREE(x_GL);
        P4EST_FREE(jac_gauss);        
      }
    }
  
  d4est_linalg_vec_axpy(1.0, M_helmf_u, prob_vecs->Au, prob_vecs->local_nodes);
  P4EST_FREE(M_helmf_u);
}


static
void apply_helmholtz
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_elliptic_problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops
)
{  
  poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops);
  
  double* M_helmf_u = P4EST_ALLOC(double, prob_vecs->local_nodes);
 
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        element_data_t* ed = quad->p.user_data;

        int deg_gauss = ed->deg;
        int volume_nodes_gauss = d4est_lgl_get_nodes((P4EST_DIM), deg_gauss);
        double* r_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 0);
        double* s_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 1);
        double* t_GL = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg_gauss, 2);
        double *jac_gauss = P4EST_ALLOC(double, volume_nodes_gauss);

        double* x_GL = P4EST_ALLOC(double, volume_nodes_gauss);
        double* y_GL = P4EST_ALLOC(double, volume_nodes_gauss);
        double* z_GL = P4EST_ALLOC(double, volume_nodes_gauss);

        d4est_linalg_fill_vec(jac_gauss, ed->jacobian, volume_nodes_gauss);  
        d4est_reference_rtox_array(r_GL, ed->xyz_corner[0], ed->h, x_GL, volume_nodes_gauss);
        d4est_reference_rtox_array(s_GL, ed->xyz_corner[1], ed->h, y_GL, volume_nodes_gauss);
        d4est_reference_rtox_array(t_GL, ed->xyz_corner[2], ed->h, z_GL, volume_nodes_gauss);
        double* xyz_gauss [3] = {x_GL, y_GL, z_GL};
        
        d4est_operators_apply_fofufofvlilj_gaussnodes
          (
           d4est_ops,
           &prob_vecs->u[ed->stride],
           NULL,
           NULL,
           ed->deg,
           jac_gauss,
           xyz_gauss,
           deg_gauss,
           (P4EST_DIM),
           &M_helmf_u[ed->stride],
           helmholtz_fcn,
           NULL,
           NULL,
           NULL
          );

        P4EST_FREE(z_GL);
        P4EST_FREE(y_GL);
        P4EST_FREE(x_GL);
        P4EST_FREE(jac_gauss);        
      }
    }
  
  d4est_linalg_vec_axpy(1.0, M_helmf_u, prob_vecs->Au, prob_vecs->local_nodes);
  P4EST_FREE(M_helmf_u);
}


static
void
build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_elliptic_problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops
)
{
  apply_helmholtz(p4est, ghost, ghost_data, prob_vecs, d4est_ops);
  d4est_linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);
}


/* static */
/* void apply_helmholtz */
/* ( */
/*  p4est_t* p4est, */
/*  p4est_ghost_t* ghost, */
/*  d4est_element_data_t* ghost_data, */
/*  d4est_elliptic_problem_data_t* prob_vecs, */
/*  d4est_operators_t* d4est_ops, */
/*  d4est_geometry_t* d4est_geom */
/* ) */
/* {   */
/*   d4est_poisson_apply_aij(p4est, */
/*                                            ghost, */
/*                                            ghost_data, */
/*                                            prob_vecs, */
/*                                            d4est_ops, */
/*                                            d4est_geom); */
  
/*   double* M_helmf_u = P4EST_ALLOC(double, prob_vecs->local_nodes); */

/*   for (p4est_topidx_t tt = p4est->first_local_tree; */
/*        tt <= p4est->last_local_tree; */
/*        ++tt) */
/*     { */
/*       p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
/*       sc_array_t* tquadrants = &tree->quadrants; */
/*       int Q = (p4est_locidx_t) tquadrants->elem_count; */
/*       for (int q = 0; q < Q; ++q) { */
/*         p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q); */
/*         d4est_element_data_t* ed = quad->p.user_data;         */
/*         d4est_operators_apply_fofufofvlilj_gaussnodes */
/*           ( */
/*            d4est_ops, */
/*            &prob_vecs->u[ed->nodal_stride], */
/*            NULL, */
/*            NULL, */
/*            ed->deg, */
/*            ed->J_quad, */
/*            ed->xyz_quad, */
/*            ed->deg_quad, */
/*            (P4EST_DIM), */
/*            &M_helmf_u[ed->nodal_stride], */
/*            helmholtz_fcn, */
/*            NULL, */
/*            NULL, */
/*            NULL */
/*           ); */
/*       } */
/*     } */
  
/*   d4est_linalg_vec_axpy(1.0, M_helmf_u, prob_vecs->Au, prob_vecs->local_nodes); */
/*   P4EST_FREE(M_helmf_u); */
/* } */




static
void problem_build_rhs
(
 p4est_t* p4est,
 d4est_elliptic_problem_data_t* prob_vecs,
 d4est_elliptic_eqns_t* prob_fcns,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 d4est_operators_t* d4est_ops
)
{
  double* f = P4EST_ALLOC(double, prob_vecs->local_nodes);
  element_data_init_node_vec
    (
     p4est,
     f,
     f_fcn,
     d4est_ops
    );
  
  element_data_apply_mij_on_vec
    (
     p4est,
     f,
     prob_vecs->rhs,
     d4est_ops
    );
  
  int local_nodes = prob_vecs->local_nodes;
  double* u_eq_0 = P4EST_ALLOC_ZERO(double, local_nodes);
  double* tmp = prob_vecs->u;
  
  prob_vecs->u = u_eq_0;
  apply_helmholtz(p4est, ghost, ghost_data, prob_vecs, d4est_ops);
  d4est_linalg_vec_axpy(-1., prob_vecs->Au, prob_vecs->rhs, local_nodes);

  prob_vecs->u = tmp;
  P4EST_FREE(u_eq_0);
  P4EST_FREE(f);
}


static
problem_input_t
problem_input
(
 const char* input_file
)
{
  int num_of_options = 6;
  
  problem_input_t input;
  input.ip_flux_penalty = -1;
  input.degree = -1;
  input.endlevel = -1; 
  input.gamma_h = -1;
  input.gamma_p = -1;
  input.sigma = -1;
  
  input.count = 0;
  
  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  mpi_assert(input.count == num_of_options);
  return input;
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

  problem_input_t input = problem_input("options.input");
  
  int level;
  int endlevel = input.endlevel;
  int degree = input.degree;         
  double gamma_h = input.gamma_h;
  double gamma_p = input.gamma_p;
  double sigma = input.sigma;
  
  d4est_poisson_flux_sipg_params_t ip_flux_params;
  ip_flux_params.sipg_penalty_prefactor = input.ip_flux_penalty;
  ip_flux_params.sipg_penalty_fcn = sipg_flux_vector_calc_penalty_maxp2_over_minh;

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

  d4est_elliptic_problem_data_t prob_vecs;
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
  /* d4est_linalg_fill_vec(u, 0., prob_vecs.local_nodes); */
  element_data_init_node_vec(p4est, u, analytic_solution_fcn, d4est_ops);
  d4est_elliptic_eqns_t prob_fcns;

  prob_fcns.apply_lhs = apply_helmholtz_matrix_2;
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
     6,
     hp_amr_smooth_pred_get_sigaverage_marker(&sigma)
    );


  hp_amr_scheme_t* scheme_uni = hp_amr_uniform();

  
  element_data_init_node_vec(p4est, u_analytic, analytic_solution_fcn, d4est_ops);    
  d4est_linalg_vec_axpy(-1., u, u_analytic, local_nodes);

    
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
       ip_flux_params.sipg_penalty_prefactor,
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
    
    d4est_linalg_vec_axpy(-1., u, u_analytic, local_nodes);

    element_data_store_nodal_vec_in_vertex_array
      (
       p4est,
       u_analytic,
       u_error_vertex
      );

    d4est_linalg_vec_fabs(u_error_vertex, (P4EST_CHILDREN)*p4est->local_num_quadrants);
    
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

    estimator_stats_print(stats);
    
    hp_amr(p4est,
           d4est_ops,
           &u,
           &stats,
           (level > 1) ? scheme : scheme_uni,
           0
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

    element_data_get_level_range(p4est, &min_level, &max_level);
    printf("number of elements = %d, [min_top_level, max_top_level] = [%d,%d]\n", p4est->local_num_quadrants, min_level, max_level);
    
    
    /* need to do a reduce on min,max_level before supporting multiple proc */
    /* mpi_assert(proc_size == 1); */
    int num_of_levels = max_level + 1;
    
    int vcycle_iter = 3;
    /* int vcycle_iter = 1; */
      double vcycle_rtol = 1e-9;
      double vcycle_atol = 0.;
      int smooth_iter = 15;
      /* int smooth_iter = 1; */
      int cg_eigs_iter = 15;
     /* int cg_eigs_iter = 1; */
      /* double max_eig_factor = 1.1; */
      double max_eig_factor = 1.0;
      int max_eig_reuse = 1;
      double lmax_lmin_rat = 30.;
      int coarse_iter = 10000;
      /* int coarse_iter = 1; */
      double coarse_rtol = 1e-10;
      int save_vtk_snapshot = 0;
      int perform_checksum = 0;
      int cg_eigs_use_zero_vec_as_initial = 0;
      
      multigrid_data_t* mg_data
        = multigrid_data_init
        (
         world_rank,
         num_of_levels,
         vcycle_iter,
         vcycle_rtol,
         vcycle_atol,
         smooth_iter,
         cg_eigs_iter,
         max_eig_factor,
         max_eig_reuse,
         lmax_lmin_rat,
         CG,
         coarse_iter,
         coarse_rtol,
         save_vtk_snapshot,
         perform_checksum,
         RESIDUAL_INFO,
         d4est_ops,
         cg_eigs_use_zero_vec_as_initial
        );

      
      multigrid_matrix_op_t* matrix_op = multigrid_matrix_operator_init(mg_data, p4est, d4est_ops, NULL);
      
      multigrid_matrix_fofu_fofv_mass_operator_setup_deg_quad_eq_deg(p4est, d4est_ops, NULL, NULL, helmholtz_fcn, NULL, NULL, NULL, matrix_op);


      element_data_init_node_vec(p4est, u_analytic, analytic_solution_fcn, d4est_ops);
      /* util_print_matrix(u_analytic, prob_vecs.local_nodes, 1, " u_analytic = ", 0); */

      /* element_data_print(p4est); */

      multigrid_data_set_analytical_solution
        (
         mg_data,
         analytic_solution_fcn
        );

      /* int local_matrix_nodes = element_data_get_local_matrix_nodes(p4est); */

      /* DEBUG_PRINT_ARR_DBL(matrix_op->matrix_at0, local_matrix_nodes); */
      prob_vecs.user = matrix_op;

      
      /* printf("STOP HERE DUCKERS\n"); */
      
      /* multigrid_solve */
      /*   ( */
      /*    p4est, */
      /*    &prob_vecs, */
      /*    &prob_fcns, */
      /*    mg_data, */
      /*    &ghost, */
      /*    &ghost_data */
      /*   ); */

      krylov_pc_t* pc = krylov_pc_multigrid_create(mg_data);

                
      krylov_petsc_info_t info =
        krylov_petsc_solve
        (
         p4est,
         &prob_vecs,
         (void*)&prob_fcns,
         &ghost,
         (void**)&ghost_data,
         d4est_ops,
         NULL,
         input_file,
         pc
        );
  
      multigrid_data_destroy(mg_data);
      multigrid_matrix_operator_destroy(matrix_op);
      krylov_pc_multigrid_destroy(pc);      

    d4est_linalg_vec_axpy(-1., u, u_analytic, local_nodes);
    
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
    
  d4est_linalg_vec_axpy(-1., u, u_analytic, local_nodes);

  element_data_store_nodal_vec_in_vertex_array
    (
     p4est,
     u_analytic,
     u_error_vertex
    );

  d4est_linalg_vec_fabs(u_error_vertex, (P4EST_CHILDREN)*p4est->local_num_quadrants);
  
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
