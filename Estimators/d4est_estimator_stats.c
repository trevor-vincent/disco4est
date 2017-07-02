#define _GNU_SOURCE
#include <d4est_estimator_stats.h>
#include <d4est_element_data.h>
#include <util.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <sc_reduce.h>

void d4est_estimator_stats_compute
(
 p4est_t* p4est,
 d4est_estimator_stats_t* stats,
 int curved
)
{
  double* eta2 = P4EST_ALLOC(double, p4est->local_num_quadrants);
  
  double total_eta2 = 0.;
  int k = 0;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q, ++k) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        eta2[k] = (((d4est_element_data_t*)(quad->p.user_data))->local_estimator);
        total_eta2 += eta2[k];
      }
    }

  stats->tree = -1;
  stats->mpirank = p4est->mpirank;
  d4est_estimator_stats_compute_stats
    (
     stats,
     eta2,
     p4est->local_num_quadrants,
     total_eta2
    );
  
  P4EST_FREE(eta2);
}

double
d4est_estimator_stats_compute_per_bin
(
 p4est_t* p4est,
 d4est_estimator_stats_t** stats,
 int num_bins,
 int(*in_bin)(d4est_element_data_t*,int)
)
{
  double* eta2 = P4EST_ALLOC_ZERO(double, p4est->local_num_quadrants);
  double local_eta2 = 0.;

  for (int b = 0; b < num_bins; b++){
    double total_eta2_per_bin = 0.;
    int bsize = 0;
    
    for (p4est_topidx_t tt = p4est->first_local_tree;
         tt <= p4est->last_local_tree;
         ++tt)
      {
        p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
        sc_array_t* tquadrants = &tree->quadrants;
        int Q = (p4est_locidx_t) tquadrants->elem_count;

        for (int q = 0; q < Q; ++q) {
          p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
          if(in_bin((d4est_element_data_t*)(quad->p.user_data),b)){
            eta2[bsize] = (((d4est_element_data_t*)(quad->p.user_data))->local_estimator);
            total_eta2_per_bin += eta2[bsize];
            local_eta2 += eta2[bsize];
            bsize++;
          }
        }

        stats[b]->tree = -1;
        stats[b]->mpirank = p4est->mpirank;
        d4est_estimator_stats_compute_stats
          (
           stats[b],
           eta2,
           bsize,
           total_eta2_per_bin
          );

      }
  }
  
  P4EST_FREE(eta2);
  return local_eta2;
}

double d4est_estimator_stats_compute_per_tree(p4est_t* p4est, d4est_estimator_stats_t** stats, int curved){

  double* eta2 = P4EST_ALLOC_ZERO(double, p4est->local_num_quadrants);
  double local_eta2 = 0.;
  
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      double total_eta2_per_tree = 0.;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        if(curved)
          eta2[q] = (((d4est_element_data_t*)(quad->p.user_data))->local_estimator);
        else{
          mpi_abort("Soon to be deprecated\n");
        }
        total_eta2_per_tree += eta2[q];
        local_eta2 += eta2[q];
      }

      stats[tt]->tree = tt;
      stats[tt]->mpirank = p4est->mpirank;
      d4est_estimator_stats_compute_stats
        (
         stats[tt],
         eta2,
         Q,
         total_eta2_per_tree
        );

    }
  
  P4EST_FREE(eta2);
  return local_eta2;
}


void
d4est_estimator_stats_compute_stats
(
 d4est_estimator_stats_t* stats,
 double* eta2,
 int sample_size,
 double total_eta2
)
{
  util_sort_double(eta2, sample_size);

  if (sample_size > 0){

  stats->total = total_eta2;
  stats->mean = total_eta2/(double)sample_size;
  stats->std = 0.;
  stats->max = eta2[sample_size-1];
  stats->min = eta2[0];
  stats->sample_size = sample_size;
  
  for (int i = 0; i < sample_size; i++){
    stats->std += (eta2[i] - stats->mean)*(eta2[i] - stats->mean);
  }
  if (stats->sample_size != 1)
    stats->std = stats->std/(double)(stats->sample_size-1.);
  else
    stats->std = 0.;
  stats->std = sqrt(stats->std);
 
  stats->p5 =  eta2[(int)(((double)sample_size)*.95)];
  stats->p10 = eta2[(int)(((double)sample_size)*.9)];
  stats->p15 = eta2[(int)(((double)sample_size)*.85)];
  stats->p20 = eta2[(int)(((double)sample_size)*.8)];
  stats->p25 = eta2[(int)(((double)sample_size)*.75)];
  
  }
  else {
    stats->total = 0;
    stats->mean = -1;
    stats->max = -1;
    stats->min = -1;
    stats->std = -1;
    stats->p5 = -1;
    stats->p10 = -1;
    stats->p15 = -1;
    stats->p20 = -1;
  }

}

void
d4est_estimator_stats_compute_max_percentiles_across_proc
(
 d4est_estimator_stats_t** stats,
 int num_bins
){
  mpi_assert(num_bins < MAX_BINS);
  double local_percentiles [MAX_BINS*5];
  double global_percentile_maxes [MAX_BINS*5];

  for (int i = 0; i < num_bins; i++){
    local_percentiles[i*5 + 0] = stats[i]->p5;
    local_percentiles[i*5 + 1] = stats[i]->p10;
    local_percentiles[i*5 + 2] = stats[i]->p15;
    local_percentiles[i*5 + 3] = stats[i]->p20;
    local_percentiles[i*5 + 4] = stats[i]->p25;
  }
    
  sc_allreduce
    (
     &local_percentiles,
     &global_percentile_maxes,
     num_bins*5,
     sc_MPI_DOUBLE,
     sc_MPI_MAX,
     sc_MPI_COMM_WORLD
    );


  for (int i = 0; i < num_bins; i++){
    stats[i]->p5 = global_percentile_maxes[i*5 + 0];
    stats[i]->p10 = global_percentile_maxes[i*5 + 1];
    stats[i]->p15 = global_percentile_maxes[i*5 + 2];
    stats[i]->p20 = global_percentile_maxes[i*5 + 3];
    stats[i]->p25 = global_percentile_maxes[i*5 + 4];
  }
  
}


double
d4est_estimator_stats_get_percentile
(
 d4est_estimator_stats_t* stats,
 int percentile
){

    if (percentile == 5)
      return stats->p5;
    else if (percentile == 10)
      return stats->p10;
    else if (percentile == 15)
      return stats->p15;
    else if (percentile == 20)
      return stats->p20;
    else if (percentile == 25)
      return stats->p25;
    else{
      mpi_abort("[D4EST_ERROR]: Not a supported percentile");
      return -1;
    }
}


void d4est_estimator_stats_print(d4est_estimator_stats_t* stats){
  printf("mpirank = %d\n", stats->mpirank);
  if(stats->tree == -1){
    printf("tree = ALL TREES\n");
  }
  else {
    printf("tree = %d\n", stats->tree);
  }
  printf("sample_size = %d\n", stats->sample_size);   
  printf("mean = %.25f\n", stats->mean);
  printf("std = %.25f\n", stats->std);
  printf("max = %.25f\n", stats->max);
  printf("min = %.25f\n", stats->min);
  printf("p5 = %.25f\n", stats->p5);
  printf("p10 = %.25f\n", stats->p10);
  printf("p15 = %.25f\n", stats->p15);
  printf("p20 = %.25f\n", stats->p20);
}

