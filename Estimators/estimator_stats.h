#ifndef ESTIMATOR_STATS_H
#define ESTIMATOR_STATS_H 

#include <pXest.h>
#include <curved_element_data.h>
#define MAX_BINS 50

typedef struct {

  int mpirank;
  int tree;
  
  double total;
  double mean;
  double std;
  double max;
  double min;
  /* double median; */
  /* double Q1; */
  /* double Q3; */
  /* double IQR; */
  /* double customQ; */
  double p5;
  double p10;
  double p15;
  double p20;
  double p25;
  
  /* bin size from scott's rule */
  /* double bin_size_scott; */
  /* double bin_size_freedman; */
  /* int num_bins_scott; */
  /* int num_bins_freedman; */
  int sample_size;
  
} estimator_stats_t;

/* This file was automatically generated.  Do not edit! */
void estimator_stats_write_to_file(p4est_t *p4est,estimator_stats_t *stats,const char *file_name_stats,const char *file_name_estimator_prefix,int curved,int mpi_rank,int print_gnuplot);
void estimator_stats_print(estimator_stats_t *stats);
double estimator_stats_get_percentile(estimator_stats_t *stats,int percentile);
double estimator_stats_compute_per_tree(p4est_t *p4est,estimator_stats_t **stats,int curved);
void estimator_stats_compute_stats(estimator_stats_t *stats,double *eta2,int sample_size,double total_eta2);
void estimator_stats_compute(p4est_t *p4est,estimator_stats_t *stats,int curved);
void
estimator_stats_compute_max_percentiles_across_proc
(
 estimator_stats_t** stats,
 int num_trees
);


double
estimator_stats_compute_per_bin
(
 p4est_t* p4est,
 estimator_stats_t** stats,
 int num_bins,
 int(*in_bin)(curved_element_data_t*,int)
);

#endif
