#ifndef D4EST_ESTIMATOR_STATS_H
#define D4EST_ESTIMATOR_STATS_H 

#include <pXest.h>
#include <d4est_element_data.h>
#define MAX_BINS 50

typedef struct {

  int mpirank;
  int tree;
  
  double total;
  double mean;
  double std;
  double max;
  double min;
  double p5;
  double p10;
  double p15;
  double p20;
  double p25;
  int sample_size;
  
} d4est_estimator_stats_t;

void d4est_estimator_stats_print(d4est_estimator_stats_t *stats);
double d4est_estimator_stats_get_percentile(d4est_estimator_stats_t *stats,int percentile);
void d4est_estimator_stats_compute_max_percentiles_across_proc(d4est_estimator_stats_t **stats,int num_bins);
double d4est_estimator_stats_compute_per_tree(p4est_t *p4est,d4est_estimator_stats_t **stats,int curved);
double d4est_estimator_stats_compute_per_bin(p4est_t *p4est,d4est_estimator_stats_t **stats,int num_bins,int(*in_bin)(d4est_element_data_t *,int));
void d4est_estimator_stats_compute_stats(d4est_estimator_stats_t *stats,double *eta2,int sample_size,double total_eta2);
void d4est_estimator_stats_compute(p4est_t *p4est,d4est_estimator_stats_t *stats,int curved);

#endif
