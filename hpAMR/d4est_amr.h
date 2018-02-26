#ifndef D4EST_AMR_H
#define D4EST_AMR_H

#include <pXest.h>
#include <d4est_operators.h>
#include <d4est_estimator_stats.h>

typedef enum {AMR_NOT_SET,
              AMR_SMOOTH_PRED,
              AMR_CUSTOM,
              AMR_UNIFORM_H,
              AMR_UNIFORM_P,
              AMR_RANDOM_H,
              AMR_RANDOM_HP} d4est_amr_scheme_type_t;

/**
 * @file   amr.c
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

typedef struct d4est_amr_scheme d4est_amr_scheme_t;

struct d4est_amr_scheme {
  d4est_amr_scheme_type_t amr_scheme_type;

  void (*post_balance_callback) (p4est_t*, void*);
  void (*pre_refine_callback) (p4est_t*, void*);
  p4est_replace_t refine_replace_callback_fcn_ptr;
  p4est_replace_t balance_replace_callback_fcn_ptr;
  p4est_iter_volume_t mark_elements;
  void (*destroy)(d4est_amr_scheme_t*);
  void* amr_scheme_data;
};

typedef struct {

  /* externally set */
  int num_of_amr_steps;
  int max_degree;

  /* internal use */
  int mpirank;
  int initial_log_size;
  int balance_log_size;
  int* refinement_log;
  int* initial_log;
  int* balance_log;

  d4est_amr_scheme_t* scheme;
  d4est_operators_t* d4est_ops;
  d4est_estimator_stats_t** d4est_estimator_stats;
  double* d4est_estimator;
  
} d4est_amr_t;
/* This file was automatically generated.  Do not edit! */
void d4est_amr_destroy(d4est_amr_t *d4est_amr);
void d4est_amr_step(p4est_t *p4est,p4est_ghost_t **ghost,d4est_element_data_t **ghost_data,d4est_operators_t *d4est_ops,d4est_amr_t *d4est_amr,double **field,double *d4est_estimator,d4est_estimator_stats_t **stats);
d4est_amr_t *d4est_amr_custom_init(p4est_t *p4est,int max_degree,int num_of_amr_steps,void(*d4est_amr_custom_mark_elements)(p4est_iter_volume_info_t *,void *),void *user);
d4est_amr_t *d4est_amr_init_random_hp(p4est_t *p4est,int max_degree,int num_of_amr_steps);
d4est_amr_t *d4est_amr_init_uniform_p(p4est_t *p4est,int max_degree,int num_of_amr_steps);
d4est_amr_t *d4est_amr_init_uniform_h(p4est_t *p4est,int max_degree,int num_of_amr_steps);
d4est_amr_t *d4est_amr_init(p4est_t *p4est,const char *input_file,void *scheme_data);
void d4est_amr_input(const char *input_file,d4est_amr_t *d4est_amr);
 
#endif
