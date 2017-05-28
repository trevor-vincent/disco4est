#ifndef HO_H
#define HO_H 

#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"
#include "../dGMath/d4est_operators.h"

typedef struct multigrid_data_t multigrid_data_t;



typedef struct {

  int hrefine; /* holds a flag: 0 = no refinement, 1 = h refinement, 2 = coarsened then refined after 2:1 balance */
  int degh [(P4EST_CHILDREN)]; /* holds the degrees of the children elements */
  int degH;  /* holds the degree of the coarsened element */

} multigrid_refine_data_t;

typedef enum
  {
    START,
    PRE_V,
    DOWNV_PRE_SMOOTH,
    DOWNV_POST_SMOOTH,
    DOWNV_PRE_COARSEN,
    DOWNV_POST_COARSEN,
    DOWNV_PRE_BALANCE,
    DOWNV_POST_BALANCE,
    DOWNV_PRE_RESTRICTION,
    DOWNV_POST_RESTRICTION,
    COARSE_PRE_SOLVE,
    COARSE_POST_SOLVE,
    UPV_PRE_REFINE,
    UPV_POST_REFINE,
    UPV_PRE_SMOOTH,
    UPV_POST_SMOOTH,
    POST_V,
    POST_RESIDUAL_UPDATE,
    END
  } multigrid_state_t;

/* typedef */
/* void */
/* (*multigrid_user_setup_callback_fcn_t) */
/* ( */
/*  p4est_t*, */
/*  multigrid_data_t*, /\* current level *\/ */
/*  void* */
/* ); */

typedef
void
(*multigrid_update_callback_fcn_t)
(
 p4est_t*,
 int, /* current level */
 problem_data_t*
);

typedef
void
(*multigrid_restrict_callback_fcn_t)
(
 p4est_iter_volume_info_t*,
 void*,
 int*,
 int,
 int
);

typedef
void
(*multigrid_prolong_callback_fcn_t)
(
 p4est_t*,
 p4est_topidx_t,
 p4est_quadrant_t*
);




typedef
void
(*multigrid_prolong_callback_fcn_t)
(
 p4est_t*,
 p4est_topidx_t,
 p4est_quadrant_t*
);


typedef
void
(*multigrid_smoother_fcn_t)
(
 p4est_t*,
 problem_data_t*,
 weakeqn_ptrs_t*,
 double*,
 int
);

typedef
void
(*multigrid_bottom_solver_fcn_t)
(
 p4est_t*,
 problem_data_t*,
 weakeqn_ptrs_t*,
 double*
);

typedef struct {

  multigrid_update_callback_fcn_t update;
  void* user;

} multigrid_logger_t;

typedef struct {

  multigrid_update_callback_fcn_t update;
  multigrid_smoother_fcn_t smooth;
  void* user;

} multigrid_smoother_t;

typedef struct {

  multigrid_update_callback_fcn_t update;
  multigrid_bottom_solver_fcn_t solve;
  void* user;

} multigrid_bottom_solver_t;

typedef struct {

  multigrid_prolong_callback_fcn_t mg_prolong_user_callback;
  multigrid_restrict_callback_fcn_t mg_restrict_user_callback;
  multigrid_update_callback_fcn_t update;
  /* multigrid_user_setup_callback_fcn_t setup; */
  void* user;
  
} multigrid_user_callbacks_t;

typedef struct {


  int(*get_local_nodes)(p4est_t*);
  multigrid_update_callback_fcn_t update;
  void* user;
  p4est_ghost_t** ghost;
  void** ghost_data;
  
  /* extra stuff for use with curved infrastructure */
  d4est_geometry_t* d4est_geom;
  void(*element_data_init_user_fcn)(void*,void*);
  d4est_geometry_storage_t** geometric_factors;
  
} multigrid_element_data_updater_t;


struct multigrid_data_t {
 
  /* ******* REQUIRED EXTERNAL PARAMETERS ******* */
  /* ******* SET BY FUNCTION CALL ******** */
  int num_of_levels;

  /* ******* SET BY INPUT FILE ******** */
  int vcycle_imax; /* max number of vcycles */
  char smoother_name [50];
  char bottom_solver_name [50];
  double vcycle_rtol; /* residual tolerance for termination */
  double vcycle_atol;

  /* ******* INTERNAL PARAMETERS ******* */
  
  multigrid_state_t mg_state;
  multigrid_refine_data_t* coarse_grid_refinement;
  d4est_operators_t* d4est_ops;
  
  /* Helper strides */
  int stride;
  int temp_stride;
  int fine_stride;
  int coarse_stride;

  /* Helper sizes */
  int fine_nodes;
  int coarse_nodes;

  /* Helper alias */
  double* intergrid_ptr;

  /* Residuals and Vcycle Info */
  double vcycle_r2_local_current;
  double vcycle_r2_global_current;
  double vcycle_r2_global_last;
  double vcycle_r2_global_stoptol;
  int vcycle_num_finished;
  
  /* INTERNAL COMPONENTS */
  multigrid_smoother_t* smoother;
  multigrid_logger_t* logger;
  multigrid_bottom_solver_t* bottom_solver;
  multigrid_user_callbacks_t* user_callbacks;
  multigrid_element_data_updater_t* elem_data_updater;

  /* INTERNAL PARAMETERS FOR LOGGING*/
  double* Ae_at0;
  double* err_at0;
  double* rres_at0;
  double* res_at0;
  int* elements_on_level_of_multigrid;
  int* elements_on_level_of_surrogate_multigrid;
  int* nodes_on_level_of_multigrid;
  int* nodes_on_level_of_surrogate_multigrid;
  
};


void
multigrid_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 multigrid_data_t* mg_data
);

void multigrid_data_destroy
(multigrid_data_t*);


multigrid_data_t*
multigrid_data_init
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 int num_of_levels,
 multigrid_logger_t* logger,
 multigrid_user_callbacks_t* user_callbacks,
 multigrid_element_data_updater_t* updater,
 const char* input_file
);

void
multigrid_get_level_range
(
 p4est_t* p4est,
 int* min_level,
 int* max_level
);

#endif
