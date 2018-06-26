#define _GNU_SOURCE
#include <d4est_estimator_stats.h>
#include <d4est_element_data.h>
#include <d4est_util.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <sc_reduce.h>

void d4est_estimator_stats_compute
(
 p4est_t* p4est,
 double* estimator,
 d4est_estimator_stats_t* stats,
 int compute_percentile,
 int compute_mean,
 int compute_max
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
        eta2[k] = estimator[k];
        total_eta2 += eta2[k];
      }
    }

  stats->tree = -1;
  stats->mpirank = p4est->mpirank;
  d4est_estimator_stats_compute_aux
    (
     p4est,
     stats,
     eta2,
     p4est->local_num_quadrants,
     total_eta2,
     compute_percentile,
     compute_mean,
     compute_max
    );
  
  P4EST_FREE(eta2);
}

/* double */
/* d4est_estimator_stats_compute_per_bin */
/* ( */
/*  p4est_t* p4est, */
/*  double* estimator, */
/*  d4est_estimator_stats_t* stats, */
/*  int num_bins, */
/*  int(*in_bin)(d4est_element_data_t*,int) */
/* ) */
/* { */
/*   double* eta2 = P4EST_ALLOC_ZERO(double, p4est->local_num_quadrants); */
/*   double local_eta2 = 0.; */

/*   for (int b = 0; b < num_bins; b++){ */
/*     double total_eta2_per_bin = 0.; */
/*     int bsize = 0; */

/*     int k = 0; */
/*     for (p4est_topidx_t tt = p4est->first_local_tree; */
/*          tt <= p4est->last_local_tree; */
/*          ++tt) */
/*       { */
/*         p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
/*         sc_array_t* tquadrants = &tree->quadrants; */
/*         int Q = (p4est_locidx_t) tquadrants->elem_count; */

/*         for (int q = 0; q < Q; ++q) { */
/*           p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q); */
/*           if(in_bin((d4est_element_data_t*)(quad->p.user_data),b)){ */
/*             eta2[bsize] = estimator[k]; */
/*             total_eta2_per_bin += eta2[bsize]; */
/*             local_eta2 += eta2[bsize]; */
/*             bsize++; */
/*           }     */
/*           k++; */
/*         } */

/*         stats[b].tree = -1; */
/*         stats[b].mpirank = p4est->mpirank; */
/*         d4est_estimator_stats_compute_aux */
/*           ( */
/*            &stats[b], */
/*            eta2, */
/*            bsize, */
/*            total_eta2_per_bin */
/*           ); */

/*       } */
/*   } */
  
/*   P4EST_FREE(eta2); */
/*   return local_eta2; */
/* } */

/* double d4est_estimator_stats_compute_per_tree */
/* ( */
/*  p4est_t* p4est, */
/*  double* estimator, */
/*  d4est_estimator_stats_t** stats */
/* ) */
/* { */
/*   double* eta2 = P4EST_ALLOC_ZERO(double, p4est->local_num_quadrants); */
/*   double local_eta2 = 0.; */

/*   int k = 0; */
/*   for (p4est_topidx_t tt = p4est->first_local_tree; */
/*        tt <= p4est->last_local_tree; */
/*        ++tt) */
/*     { */
/*       p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
/*       sc_array_t* tquadrants = &tree->quadrants; */
/*       int Q = (p4est_locidx_t) tquadrants->elem_count; */
/*       double total_eta2_per_tree = 0.; */
/*       for (int q = 0; q < Q; ++q) { */
/*         p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q); */
/*         eta2[q] = estimator[k]; */
/*         total_eta2_per_tree += eta2[q]; */
/*         local_eta2 += eta2[q]; */
/*         k++; */
/*       } */

/*       stats[tt]->tree = tt; */
/*       stats[tt]->mpirank = p4est->mpirank; */
/*       d4est_estimator_stats_compute_aux */
/*         ( */
/*          stats[tt], */
/*          eta2, */
/*          Q, */
/*          total_eta2_per_tree */
/*         ); */

/*     } */
  
/*   P4EST_FREE(eta2); */
/*   return local_eta2; */
/* } */


void
d4est_estimator_stats_compute_aux
(
 p4est_t* p4est,
 d4est_estimator_stats_t* stats,
 double* eta2,
 int local_size,
 double total_eta2,
 int compute_percentile,
 int compute_mean,
 int compute_max
)
{
  d4est_util_sort_double(eta2, local_size);

  stats->estimator_total = 0;
  stats->estimator_mean = -1;
  stats->estimator_max = -1;
  stats->estimator_at_percentile = -1;
  stats->percentile = -1;

  if (p4est->mpisize == 1){

    if (local_size > 0){
      stats->estimator_total = total_eta2;
      stats->estimator_mean = total_eta2/(double)local_size;
      stats->estimator_max = eta2[local_size-1];
      stats->estimator_sample_size = local_size;
      stats->estimator_at_percentile = eta2[(int)(((double)local_size)*(1.-((double)compute_percentile/100.0)))];
      /* stats->estimator_at_percentile =  eta2[(int)(((double)local_size)*.9)]; */

      /* printf("(int)(((double)local_size)*(1.-((double)compute_percentile/100.0))) = %d\n, (1.-((double)compute_percentile/100.0)) = %f\n, local_size = %d\n, compute_percentile = %d\n", (int)(((double)local_size)*(1.-((double)compute_percentile/100.0))),(1.-((double)compute_percentile/100.0)), local_size, compute_percentile); */
    }

    /* printf("estimator at percentile = %.15f\n",stats->estimator_at_percentile); */
    /* stats->min = eta2[0]; */  
    /* for (int i = 0; i < sample_size; i++){ */
    /* stats->std += (eta2[i] - stats->mean)*(eta2[i] - stats->mean); */
    /* } */
    /* if (stats->sample_size != 1) */
    /* stats->std = stats->std/(double)(stats->sample_size-1.); */
    /* else */
    /* stats->std = 0.; */
    /* stats->std = sqrt(stats->std); */
 
    /* stats->p5 =  eta2[(int)(((double)local_size)*.95)]; */
    /* stats->p10 = eta2[(int)(((double)local_size)*.9)]; */
    /* stats->p15 = eta2[(int)(((double)local_size)*.85)]; */
    /* stats->p20 = eta2[(int)(((double)local_size)*.8)]; */
    /* stats->p25 = eta2[(int)(((double)local_size)*.75)];  */
  }
  else {
    stats->estimator_sample_size = p4est->global_num_quadrants;
    d4est_util_parallel_sort(p4est->mpicomm, eta2, local_size);
    if (compute_percentile){
      stats->estimator_at_percentile = d4est_estimator_stats_get_global_percentile_parallel
                                       (
                                        p4est,
                                        compute_percentile,
                                        eta2
                                       );
    }
    if (compute_mean){
      sc_allreduce(&total_eta2, &stats->estimator_total, 1, sc_MPI_DOUBLE, sc_MPI_SUM, p4est->mpicomm);
      stats->estimator_mean = stats->estimator_total/p4est->global_num_quadrants;
    }
    if (compute_max){
      sc_allreduce(&eta2[local_size-1], &stats->estimator_max, 1, sc_MPI_DOUBLE, sc_MPI_MAX, p4est->mpicomm);
    }
    if (!compute_mean && !compute_mean && !compute_percentile){
      D4EST_ABORT("Don't call estimator stats if you're not computing anything");
    }
  }

}

/* void */
/* d4est_estimator_stats_compute_max_percentiles_across_proc */
/* ( */
/*  d4est_estimator_stats_t* stats, */
/*  int num_bins */
/* ){ */
/*   D4EST_ASSERT(num_bins < MAX_BINS); */
/*   double local_percentiles [MAX_BINS*5]; */
/*   double global_percentile_maxes [MAX_BINS*5]; */

/*   for (int i = 0; i < num_bins; i++){ */
/*     local_percentiles[i*5 + 0] = stats[i].p5; */
/*     local_percentiles[i*5 + 1] = stats[i].p10; */
/*     local_percentiles[i*5 + 2] = stats[i].p15; */
/*     local_percentiles[i*5 + 3] = stats[i].p20; */
/*     local_percentiles[i*5 + 4] = stats[i].p25; */
/*   } */
    
/*   sc_allreduce */
/*     ( */
/*      &local_percentiles, */
/*      &global_percentile_maxes, */
/*      num_bins*5, */
/*      sc_MPI_DOUBLE, */
/*      sc_MPI_MAX, */
/*      sc_MPI_COMM_WORLD */
/*     ); */

/*   for (int i = 0; i < num_bins; i++){ */
/*     stats[i].p5 = global_percentile_maxes[i*5 + 0]; */
/*     stats[i].p10 = global_percentile_maxes[i*5 + 1]; */
/*     stats[i].p15 = global_percentile_maxes[i*5 + 2]; */
/*     stats[i].p20 = global_percentile_maxes[i*5 + 3]; */
/*     stats[i].p25 = global_percentile_maxes[i*5 + 4]; */
/*   } */
  
/* } */

/* double */
/* d4est_estimator_stats_get_percentile */
/* ( */
/*  d4est_estimator_stats_t* stats, */
/*  int percentile */
/* ) */
/* { */
/*   if (percentile == 5) */
/*     return stats->p5; */
/*   else if (percentile == 10) */
/*     return stats->p10; */
/*   else if (percentile == 15) */
/*     return stats->p15; */
/*   else if (percentile == 20) */
/*     return stats->p20; */
/*   else if (percentile == 25) */
/*     return stats->p25; */
/*   else{ */
/*     D4EST_ABORT("[D4EST_ERROR]: Not a supported percentile"); */
/*     return -1; */
/*   } */
/* } */


void d4est_estimator_stats_print(d4est_estimator_stats_t* stats){
  printf("mpirank = %d\n", stats->mpirank);
  if(stats->tree == -1){
    printf("tree = ALL TREES\n");
  }
  else {
    printf("tree = %d\n", stats->tree);
  }
  printf("sample_size = %d\n", stats->estimator_sample_size);   
  printf("total = %.25f\n", stats->estimator_total);   
  printf("mean = %.25f\n", stats->estimator_mean);
  printf("max = %.25f\n", stats->estimator_max);
}

double
d4est_estimator_stats_get_global_percentile_parallel
(
 p4est_t* p4est,
 int percentile,
 double* estimator
)
{
  int global_id_from_end = (int)((double)(p4est->global_num_quadrants*percentile)/100.0);
  double estimator_at_percentile = -1;
  int stride = 0;
  int temp_rank = p4est->mpisize - 1;
  while(!(estimator_at_percentile > 0)){

    if (p4est->mpirank == temp_rank){
      if (global_id_from_end <= p4est->local_num_quadrants + stride){
        if (global_id_from_end != stride){
          estimator_at_percentile = estimator[p4est->local_num_quadrants - (global_id_from_end - stride)];
        }
        else {
          D4EST_ABORT("global_id_from_end == stride, should never happen");
          /* estimator_at_percentile = estimator[local_num_quadrants - 1];  */
        }
      }
      else {
        stride += p4est->local_num_quadrants;
      }
    }
    /* sc_MPI_Bcast(&estimator_at_percentile, sizeof(double), sc_MPI_DOUBLE, temp_rank, p4est->mpicomm); */
    /* sc_MPI_Barrier(p4est->mpicomm); */

    double estimator_at_percentile_temp = estimator_at_percentile;
    sc_allreduce(&estimator_at_percentile_temp, &estimator_at_percentile, 1, sc_MPI_DOUBLE, sc_MPI_MAX,
                 p4est->mpicomm);
    
    if (!(estimator_at_percentile > 0)){
      temp_rank--;
      if (temp_rank == -1){
        printf("percentile = %d", percentile);
        D4EST_ABORT("Could not find percentile");
      }
    }
  }

  /* sc_MPI_Barrier(p4est->mpicomm); */
  return estimator_at_percentile;
}

