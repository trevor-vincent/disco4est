#ifndef HP_AMR_H
#define HP_AMR_H 

#include "../pXest/pXest.h"
#include "../dGMath/dgmath.h"
#include "../Estimators/estimator_stats.h"

typedef struct {

  void (*post_balance_callback) (p4est_t*, void*);
  void (*pre_refine_callback) (p4est_t*, void*);
  p4est_replace_t refine_replace_callback_fcn_ptr;
  p4est_replace_t balance_replace_callback_fcn_ptr;  
  p4est_iter_volume_t iter_volume;
  void* hp_amr_scheme_data;

} hp_amr_scheme_t;


typedef struct {
  
  double* data;
  int* refinement_log;
  int refinement_log_stride;
  void* hp_amr_scheme_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase;  
  p4est_replace_t refine_replace_callback_fcn_ptr;
  p4est_replace_t balance_replace_callback_fcn_ptr;
  estimator_stats_t* estimator_stats;
  int elements_marked_for_hrefine;
  int elements_marked_for_prefine;
  
} hp_amr_data_t;

void
hp_amr
(
 p4est_t* p4est,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double** data_to_hp_refine,
 estimator_stats_t* stats,
 hp_amr_scheme_t* scheme
);

void
hp_amr_store_balance_changes
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 int num_outgoing,
 p4est_quadrant_t * outgoing[],
 int num_incoming,
 p4est_quadrant_t * incoming[]
);

void
hp_amr_set_deg_on_refine
(
 p4est_t* p4est,
 p4est_topidx_t which_tree,
 int num_outgoing,
 p4est_quadrant_t* outgoing[],
 int num_incoming,
 p4est_quadrant_t* incoming[]
);

void
hp_amr_save_to_vtk
(
 p4est_t* p4est,
 char* filename,
 int save_local_est
);



#endif
