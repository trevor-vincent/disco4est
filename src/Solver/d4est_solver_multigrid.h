#ifndef D4EST_SOLVER_MULTIGRID_H
#define D4EST_SOLVER_MULTIGRID_H 

#include <pXest.h>
#include <d4est_mesh.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_operators.h>

typedef struct  d4est_solver_multigrid_data  d4est_solver_multigrid_data_t;

typedef struct {

  int hrefine; /* holds a flag: 0 = no refinement, 1 = h refinement, 2 = coarsened then refined after 2:1 balance */
  int degh [(P4EST_CHILDREN)]; /* holds the degrees of the children elements */
  int degH;  /* holds the degree of the coarsened element */

} d4est_solver_multigrid_refine_data_t;

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
  } d4est_solver_multigrid_state_t;

typedef
void
(*d4est_solver_multigrid_update_callback_fcn_t)
(
 p4est_t*,
 int, /* current level */
 d4est_elliptic_data_t*
);

typedef
void
(*d4est_solver_multigrid_restrict_callback_fcn_t)
(
 p4est_iter_volume_info_t*,
 void*,
 int*,
 int,
 int
);

typedef
void
(*d4est_solver_multigrid_prolong_callback_fcn_t)
(
 p4est_t*,
 p4est_topidx_t,
 p4est_quadrant_t*
);




typedef
void
(*d4est_solver_multigrid_prolong_callback_fcn_t)
(
 p4est_t*,
 p4est_topidx_t,
 p4est_quadrant_t*
);


typedef
void
(*d4est_solver_multigrid_smoother_fcn_t)
(
 p4est_t*,
 d4est_elliptic_data_t*,
 d4est_elliptic_eqns_t*,
 double*,
 int
);

typedef
void
(*d4est_solver_multigrid_bottom_solver_fcn_t)
(
 p4est_t*,
 d4est_elliptic_data_t*,
 d4est_elliptic_eqns_t*,
 double*
);

typedef struct {

  d4est_solver_multigrid_update_callback_fcn_t update;
  void* user;

} d4est_solver_multigrid_logger_t;

typedef struct {

  d4est_solver_multigrid_update_callback_fcn_t update;
  void* user;

} d4est_solver_multigrid_profiler_t;

typedef struct {

  d4est_solver_multigrid_update_callback_fcn_t update;
  d4est_solver_multigrid_smoother_fcn_t smooth;
  void* user;

} d4est_solver_multigrid_smoother_t;

typedef struct {

  d4est_solver_multigrid_update_callback_fcn_t update;
  d4est_solver_multigrid_bottom_solver_fcn_t solve;
  void* user;

} d4est_solver_multigrid_bottom_solver_t;

typedef struct {

  d4est_solver_multigrid_prolong_callback_fcn_t mg_prolong_user_callback;
  d4est_solver_multigrid_restrict_callback_fcn_t mg_restrict_user_callback;
  d4est_solver_multigrid_update_callback_fcn_t update;
  void* user;
  
} d4est_solver_multigrid_user_callbacks_t;


typedef struct {

  d4est_solver_multigrid_update_callback_fcn_t update;  
  int stride;
  int levels;
  unsigned int p4est_checksums [20];
  unsigned int deg_checksums [20];
  /* char* input_file; */

} d4est_solver_multigrid_mesh_analyzer_t;

typedef struct {
  
  d4est_solver_multigrid_update_callback_fcn_t update;
  d4est_mesh_data_t** d4est_factors_on_level;
  d4est_ghost_t** d4est_ghost_on_level;
  d4est_ghost_data_t** d4est_ghost_data_on_level;
  d4est_mesh_initial_extents_t* initial_extents;
  
  /* alias to currently used geometric factors, points to an element of geometric_factors above */
  d4est_mesh_data_t* current_d4est_factors;
  d4est_ghost_t* current_d4est_ghost;
  d4est_ghost_data_t* current_d4est_ghost_data;
  
  void(*element_data_init_user_fcn)(d4est_element_data_t*,void*);
  void* user;
  
} d4est_solver_multigrid_element_data_updater_t;

struct  d4est_solver_multigrid_data {

  char* input_file;
  
  /* ******* REQUIRED EXTERNAL PARAMETERS ******* */
  /* ******* SET BY INPUT FILE ******** */
  int vcycle_imax; /* max number of vcycles */
  char smoother_name [50];
  char bottom_solver_name [50];
  double vcycle_rtol; /* residual tolerance for termination */
  double vcycle_atol;
  int use_profiler;
  int use_analyzer;
  int print_state_info;
  int print_level_info;
  int use_power_method_debug;
  double power_atol;
  double power_rtol;
  double power_imax;
  double power_imin;
  int use_p_coarsen; /* p-d4est_solver_multigrid is untested */
    
  /* ******* INTERNAL PARAMETERS ******* */
  int num_of_levels;
  int num_of_p_coarsen_levels; /* p-d4est_solver_multigrid is untested */
  d4est_solver_multigrid_state_t mg_state;
  d4est_solver_multigrid_refine_data_t* coarse_grid_refinement;
  d4est_operators_t* d4est_ops;
  d4est_geometry_t* d4est_geom;
  d4est_quadrature_t* d4est_quad;
  
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
  d4est_solver_multigrid_smoother_t* smoother;
  d4est_solver_multigrid_profiler_t* profiler;
  d4est_solver_multigrid_mesh_analyzer_t* analyzer;
  d4est_solver_multigrid_logger_t* logger;
  d4est_solver_multigrid_bottom_solver_t* bottom_solver;
  d4est_solver_multigrid_user_callbacks_t* user_callbacks;
  d4est_solver_multigrid_element_data_updater_t* elem_data_updater;

  /* INTERNAL PARAMETERS FOR LOGGING*/
  double* Ae_at0;
  double* err_at0;
  double* rres_at0;
  double* res_at0;
  int* elements_on_level_of_multigrid;
  int* elements_on_level_of_surrogate_multigrid;
  int* nodes_on_level_of_multigrid;
  int* nodes_on_level_of_surrogate_multigrid;

  /* UPDATED BY EXTERNAL FUNCTIONS */
  /* if d4est_solver_multigrid is being used as a preconditioner, keep track of linear operator updates, this is useful so that we don't need to recompute operator related
 * quantities every time we run a v-cycle, only when the linear operator is updated
 * (e.g. a new newton iteration) 
 * This is used in d4est_solver_multigrid_matrix_operator for example. */
  int linear_operator_updates;
  
};

/* This file was automatically generated.  Do not edit! */
void d4est_solver_multigrid_solve(p4est_t *p4est,d4est_elliptic_data_t *vecs,d4est_elliptic_eqns_t *fcns,d4est_solver_multigrid_data_t *mg_data);
void d4est_solver_multigrid_data_destroy(d4est_solver_multigrid_data_t *mg_data);
void d4est_solver_multigrid_set_user_callbacks(d4est_solver_multigrid_data_t *mg_data,d4est_solver_multigrid_user_callbacks_t *user_callbacks);
d4est_solver_multigrid_data_t *d4est_solver_multigrid_data_init(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_ghost_t **d4est_ghost,d4est_ghost_data_t **d4est_ghost_data,d4est_mesh_data_t *d4est_factors,d4est_mesh_initial_extents_t *initial_extents,const char *input_file);
void d4est_solver_multigrid_destroy_bottom_solver(d4est_solver_multigrid_data_t *mg_data);
void d4est_solver_multigrid_set_bottom_solver(p4est_t *p4est,const char *input_file,d4est_solver_multigrid_data_t *mg_data);
void d4est_solver_multigrid_destroy_smoother(d4est_solver_multigrid_data_t *mg_data);
void d4est_solver_multigrid_set_smoother(p4est_t *p4est,const char *input_file,d4est_solver_multigrid_data_t *mg_data);
int d4est_solver_multigrid_get_h_coarsen_levels_old(p4est_t *p4est);
int d4est_solver_multigrid_get_h_coarsen_levels_initial(p4est_t *p4est,d4est_mesh_initial_extents_t *initial_extents);
int d4est_solver_multigrid_get_h_coarsen_levels_post_initial(p4est_t *p4est,int num_of_levels_last,int has_there_been_h_refinements);
int d4est_solver_multigrid_get_p_coarsen_levels(p4est_t *p4est);

#endif
