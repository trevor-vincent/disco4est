#ifndef HO_H
#define HO_H 

#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"
#include "../dGMath/dgmath.h"

typedef struct {

  int hrefine; /* holds a flag: 0 = no refinement, 1 = h refinement, 2 = coarsened then refined after 2:1 balance */
  int degh [(P4EST_CHILDREN)]; /* holds the degrees of the children elements */
  int degH;  /* holds the degree of the coarsened element */

} multigrid_refine_data_t;

typedef struct {

  int iter;
  double rtol;
  int mpi_rank;

} multigrid_cg_params_t;

typedef struct {

  int iter;
  double lmin;
  double lmax;
  int mpi_rank;

} multigrid_cheby_params_t;

/* determine which type of verbosity one desires */
typedef enum 
  {
    NONE,
    DEBUG,
    ERR_LOG, /* print error if analytical soln is available */
    RES_LOG, /* print residual */
    ERR_AND_EIG_LOG, /* print eigenvalues and error */
    RES_AND_EIG_LOG, /* print eigenvalues and residual */
    PRINT_NOT_SET
  } multigrid_log_option_t;

typedef enum
  {
    PREV_PRE_SMOOTH,
    PREV_POST_SMOOTH,
    DOWNV_POST_COARSEN,
    DOWNV_POST_BALANCE,
    DOWNV_PRE_SMOOTH,
    DOWNV_POST_SMOOTH,
    COARSE_PRE_SOLVE,
    COARSE_POST_SOLVE,
    UPV_PRE_REFINE,
    UPV_PRE_SMOOTH,
    POSTV_PRE_SMOOTH,
    POSTV_POST_SMOOTH
  } multigrid_state_t;


typedef enum 
  {
    CHEBY, 
    CG,
    SOLVER_NOT_SET
  } multigrid_coarse_solver_t;

typedef
void
(*multigrid_update_user_callback_fcn_t)
(
 p4est_t*,
 int, /* current level */
 problem_data_t*
);

typedef
void
(*multigrid_restrict_user_callback_fcn_t)
(
 p4est_iter_volume_info_t*,
 void*
);

typedef
void
(*multigrid_prolong_user_callback_fcn_t)
(
 p4est_t*,
 p4est_topidx_t,
 p4est_quadrant_t*
);

typedef struct {
 
  /* ******* REQUIRED EXTERNAL PARAMETERS ******* */
  /* ******* user should set all of these ******** */

  int smooth_iter; /* smoothing iterations */

  int coarse_iter; /* coarse iterations */
  double coarse_rtol;

  int num_of_levels;
  
  int mpi_rank;
  int save_vtk_snapshot; /* save each grid level to vtk */
  int perform_checksum; /* perform a checksum before&after vcycle */

  int vcycle_iter; /* max number of vcycles */
  double vcycle_rtol; /* residual tolerance for termination */
  double vcycle_atol;
  
  double lmax_lmin_rat; /* lmin = lmax/lmax_lmin_rat */
  int cg_eigs_iter;
  double max_eig_factor;
  int max_eig_reuse; // reuse same eigenvalues going up V (yes if you want MG-operator to be symmetric)

  int cg_eigs_use_zero_vec_as_initial;
  
  multigrid_log_option_t log_option;
  multigrid_coarse_solver_t coarse_solver_type;
  dgmath_jit_dbase_t* dgmath_jit_dbase;

  /* NOT REQUIRED EXTERNAL PARAMETERS */
  /* Set these if you want to prolong/restrict custom (user-defined) fields */
  int user_defined_fields;
  multigrid_prolong_user_callback_fcn_t mg_prolong_user_callback;
  multigrid_restrict_user_callback_fcn_t mg_restrict_user_callback;
  multigrid_update_user_callback_fcn_t mg_update_user_callback;
  void* user_ctx;

  /* If we are running a test problem with an analytical solution */
  grid_fcn_t analytical_solution;


  
  /* ******* INTERNAL PARAMETERS ******* */
  /* ******* no need to set these ******** */
  int fine_nodes;
  int coarse_nodes;

  /* eigenvalues on each level */
  double* max_eigs;
  int solve_for_eigs;
  
  /* multigrid vectors */
  double* Ae;
  double* err;
  double* res;
  double* rres;

  double* intergrid_ptr;
  
  /* To store mesh information for each level */
  multigrid_refine_data_t* coarse_grid_refinement;

  /* Helper strides */
  int stride;
  int temp_stride;
  int fine_stride;
  int coarse_stride;
  
  double vcycle_r2local;
  int final_vcycles; /* final number of vcycles */

  double lmax; /* max eigenvalue to sweep */
  double lmin; /* min eigenvalue to sweep */

  multigrid_state_t mg_state;
  
} multigrid_data_t;


multigrid_data_t*
multigrid_data_init
(
 int mpi_rank, 
 int num_of_levels,
 int vcycle_iter,
 /* Stopping condition := r2 <= vcycle_rtol*vcycle_rtol*r2initial + vcycle_atol*vcycle_atol */
 double vcycle_rtol,
 double vcycle_atol,
 int smooth_iter,
 /* Iterations of cg to find spectral radius */
int cg_eigs_iter,
 /* Multiplier for maximum eigenvalue estimate, = 1.0 for no multiplier */
 double max_eig_factor,
 /* Reuse eigenvalue estimates from the last vcycle, =1 is most sensible */
 int max_eig_reuse,
 double lmax_lmin_rat,
 multigrid_coarse_solver_t coarse_solver_type,
 int coarse_iter,
 double coarse_rtol,
 int save_vtk_snapshot,
 int perform_checksum,
 multigrid_log_option_t log_option,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 int cg_eigs_use_zero_vec_as_initial
);

void
multigrid_data_set_user_defined_fields
(
 multigrid_data_t* mg_data,
 multigrid_prolong_user_callback_fcn_t mg_prolong_user_callback,
 multigrid_restrict_user_callback_fcn_t mg_restrict_user_callback,
 multigrid_update_user_callback_fcn_t mg_update_user_callback,
 void* user_ctx
);

void
multigrid_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 multigrid_data_t* mg_data,
 p4est_ghost_t** ghost,
 element_data_t** ghost_data
);

void multigrid_data_destroy
(multigrid_data_t*);

void
multigrid_data_set_analytical_solution
(
 multigrid_data_t* mg_data,
 grid_fcn_t analytical_solution
);

#endif
