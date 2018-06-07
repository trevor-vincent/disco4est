#ifndef MULTIGRID_PROFILER_BASIC_H
#define MULTIGRID_PROFILER_BASIC_H 

#include <pXest.h>
#include <time.h>
#include <multigrid.h>

typedef struct  {
  
  /* temporaries */
  time_t begin;
  time_t end;
  
  double time_start_pre_v;
  double time_downv_smooth;
  double time_downv_coarsen;
  double time_downv_balance;
  double time_downv_restriction;
  double time_coarse_solve;
  double time_upv_refine;
  double time_upv_smooth;
  double time_post_v_end;
  double time_total;
  
} multigrid_profiler_basic_data_t;


/* This file was automatically generated.  Do not edit! */
/* This file was automatically generated.  Do not edit! */
multigrid_profiler_t *multigrid_profiler_basic_init();
void multigrid_profiler_basic_destroy(multigrid_profiler_t *profiler);


#endif
