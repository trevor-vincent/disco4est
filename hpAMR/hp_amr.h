#ifndef HP_AMR_H
#define HP_AMR_H 

#include "../pXest/pXest.h"
#include "../dGMath/dgmath.h"
#include "../Estimators/estimator_stats.h"


/**
 * @file   hp_amr.c
 * @author Trevor Vincent <tvincent@cita.utoronto.ca>
 * @date   Mon Nov  2 01:08:28 2015
 * 
 * @brief  Here are the rules for the refinement_log
 * array :
 *
 * This is the interface for running hp-amr loops.
 * 
 * Suppose the element currently has degree p, then:
 *
 * if refinement_log[i] = p' with p' < 0
 * then we h-refine in element i with the absolute
 * value of p as the degree of the parent
 * and children octants
 *
 * if refinement_log[i] = p' with p' > 0 then we
 * p-refine if p' > p and we p-coarsen if p' < p
 * We do nothing if p' = p.
 *
 */

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
 hp_amr_scheme_t* scheme,
 int curved
);



#endif
