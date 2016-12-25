#define _GNU_SOURCE
#include "estimator_stats.h"
#include "../ElementData/curved_element_data.h"
#include "../ElementData/element_data.h"
#include "../Utilities/util.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <math.h>

void estimator_stats_compute(p4est_t* p4est, estimator_stats_t* stats, int curved){

  double* eta2 = P4EST_ALLOC(double, p4est->local_num_quadrants);
  
  double total = 0.;
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
        if(curved)
          eta2[k] = (((curved_element_data_t*)(quad->p.user_data))->local_estimator);
        else{
          eta2[k] = (((element_data_t*)(quad->p.user_data))->local_estimator);
        }
        total += eta2[k];
      }
    }

  util_sort_double(eta2, p4est->local_num_quadrants);

  /* if (p4est->mpirank == 0) */
    /* util_print_matrix(eta2, p4est->local_num_quadrants, 1, "eta2 = ", 0); */
  /* check that assumptions hold for the below calculations */
  /* mpi_assert(p4est->local_num_quadrants >= (P4EST_CHILDREN) && eta2[1] >= eta2[0] && eta2[2] >= eta2[1]); */
  double max = eta2[p4est->local_num_quadrants-1];
  double min = eta2[0];
  
  stats->total = total;
  stats->mean = total/(double)p4est->local_num_quadrants;
  stats->median = util_compute_median(eta2, p4est->local_num_quadrants);
  stats->std = 0.;
  stats->max = max;
  stats->min = min;
  stats->sample_size = p4est->local_num_quadrants;
  
  for (int i = 0; i < p4est->local_num_quadrants; i++){
    stats->std += (eta2[i] - stats->mean)*(eta2[i] - stats->mean);
  }
  if (stats->sample_size != 1)
    stats->std = stats->std/(double)(stats->sample_size-1.);
  else
    stats->std = 0.;
  stats->std = sqrt(stats->std);
  stats->bin_size_scott = 3.5*stats->std*(1./pow(stats->sample_size,1./3.));
  stats->num_bins_scott = round((stats->max - stats->min)/stats->bin_size_scott);

  stats->p5 =  eta2[(int)(((double)p4est->local_num_quadrants)*.95)];
  stats->p10 = eta2[(int)(((double)p4est->local_num_quadrants)*.9)];
  stats->p15 = eta2[(int)(((double)p4est->local_num_quadrants)*.85)];
  stats->p20 = eta2[(int)(((double)p4est->local_num_quadrants)*.8)];
    
  /* if even */
  if (stats->sample_size != 1){
    
    if (stats->sample_size % 2 == 0){
      stats->Q1 = util_compute_median(&eta2[0], stats->sample_size/2);
      stats->Q3 = util_compute_median(&eta2[stats->sample_size/2], stats->sample_size/2);
      stats->IQR = stats->Q3 - stats->Q1;
      stats->bin_size_freedman = (2.*stats->IQR)/pow((double)stats->sample_size, 1./3.);
    }
    else {
      stats->Q1 = util_compute_median(&eta2[0], (stats->sample_size-1)/2);
      stats->Q3 = util_compute_median(&eta2[(stats->sample_size+1)/2], stats->sample_size/2);
      stats->IQR = stats->Q3 - stats->Q1;
      stats->bin_size_freedman = (2.*stats->IQR)/pow((double)stats->sample_size, 1./3.);
    }
    stats->num_bins_freedman = round((stats->max - stats->min)/stats->bin_size_freedman);
  }
  else {
    stats->Q1 = 0.;
    stats->Q3 = 0.;
    stats->IQR = 0.;
    stats->bin_size_freedman = 0.;
  }
  P4EST_FREE(eta2);
}

double
estimator_stats_get_percentile
(
 estimator_stats_t* stats,
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
      return stats->Q3;
    else{
      mpi_abort("[D4EST_ERROR]: Not a supported percentile");
      return -1;
    }
}


void estimator_stats_print(estimator_stats_t* stats, int mpi_rank){
  if (mpi_rank == 0){
    printf("sample_size = %d\n", stats->sample_size);   
    printf("mean = %.25f\n", stats->mean);
    printf("std = %.25f\n", stats->std);
    printf("max = %.25f\n", stats->max);
    printf("min = %.25f\n", stats->min);
    printf("Q1 = %.25f\n", stats->Q1);
    printf("Q3 = %.25f\n", stats->Q3);
    printf("p5 = %.25f\n", stats->p5);
    printf("p10 = %.25f\n", stats->p10);
    printf("p15 = %.25f\n", stats->p15);
    printf("p20 = %.25f\n", stats->p20);
    printf("IQR = %.25f\n", stats->IQR);
    printf("num_bins_scott = %d\n", stats->num_bins_scott);
    printf("num_bins_freedman = %d\n", stats->num_bins_freedman);
    printf("bins_size_scot = %.25f\n", stats->bin_size_scott);
    printf("bins_size_freedman = %.25f\n", stats->bin_size_freedman);
  }
}

void estimator_stats_write_to_file
(
 p4est_t* p4est,
 estimator_stats_t* stats,
 const char* file_name_stats,
 const char* file_name_estimator_prefix,
 int curved,
 int mpi_rank,
 int print_gnuplot
)
{
  if (mpi_rank == 0){
    FILE* f = fopen(file_name_stats, "w");
    if (f == NULL){
      mpi_abort("Error opening file");
    }
    fprintf(f,"sample_size = %d\n", stats->sample_size);   
    fprintf(f,"mean = %.25f\n", stats->mean);
    fprintf(f,"std = %.25f\n", stats->std);
    fprintf(f,"max = %.25f\n", stats->max);
    fprintf(f,"min = %.25f\n", stats->min);
    fprintf(f,"p5 = %.25f\n", stats->p5);
    fprintf(f,"p10 = %.25f\n", stats->p10);
    fprintf(f,"p15 = %.25f\n", stats->p15);
    fprintf(f,"p20 = %.25f\n", stats->p20);
    fprintf(f,"Q1 = %.25f\n", stats->Q1);
    fprintf(f,"Q3 = %.25f\n", stats->Q3); 
    fprintf(f,"IQR = %.25f\n", stats->IQR);
    fprintf(f,"num_bins_scott = %d\n", stats->num_bins_scott);
    fprintf(f,"num_bins_freedman = %d\n", stats->num_bins_freedman);
    fprintf(f,"bins_size_scot = %.25f\n", stats->bin_size_scott);
    fprintf(f,"bins_size_freedman = %.25f\n", stats->bin_size_freedman);
    fclose(f);

    
    if(print_gnuplot){

      char* gnuplot_freedman_save_as = NULL;
      D4EST_ASPRINTF(gnuplot_freedman_save_as, "%s_%s","freedman_gnuplot", file_name_estimator_prefix );
      FILE* file_freedman = fopen(gnuplot_freedman_save_as, "w");  
      if (file_freedman == NULL){
        mpi_abort("Error opening file");
      }
      free(gnuplot_freedman_save_as);
      
      fprintf(file_freedman,"reset \n");
      fprintf(file_freedman,"n=%d #number of intervals\n",stats->num_bins_freedman);
      fprintf(file_freedman,"max=%.25f #max value\n",stats->max);
      fprintf(file_freedman,"min=%.25f #min value\n",stats->min);
      fprintf(file_freedman,"width=%.25f #interval width\n",stats->bin_size_freedman);
      fprintf(file_freedman,"#function used to map a value to the intervals\n");
      fprintf(file_freedman,"hist(x,width)=width*floor(x/width)+width/2.0\n");
      fprintf(file_freedman,"set term png #output terminal and file\n");
      fprintf(file_freedman,"set output \"freedman_histogram.png\"\n");
      fprintf(file_freedman,"set xrange [min:max]\n");
      fprintf(file_freedman,"set yrange [0:]\n");
      fprintf(file_freedman,"#to put an empty boundary around the\n");
      fprintf(file_freedman,"#data inside an autoscaled graph.\n");
      fprintf(file_freedman,"set offset graph 0.05,0.05,0.05,0.0\n");
      fprintf(file_freedman,"set xtics min,(max-min)/5,max\n");
      fprintf(file_freedman,"set boxwidth width*0.9\n");
      fprintf(file_freedman,"set style fill solid 0.5 #fillstyle\n");
      fprintf(file_freedman,"set tics out nomirror\n");
      fprintf(file_freedman,"set xlabel \"eta2\"\n");
      fprintf(file_freedman,"set ylabel \"Frequency\"\n");
      fprintf(file_freedman,"#count and plot\n");
      fprintf(file_freedman,"plot \"%s\" u (hist($1,width)):(1.0) smooth freq w boxes lc rgb \"green\" notitle\n",file_name_estimator_prefix);
      fclose(file_freedman);
  

      char* gnuplot_scott_save_as = NULL;
      D4EST_ASPRINTF(gnuplot_scott_save_as, "%s_%s", "scott_gnuplot", file_name_estimator_prefix);
      FILE* file_scott = fopen(gnuplot_scott_save_as, "w");  
      if (file_scott == NULL){
        mpi_abort("Error opening file");
      }
      free(gnuplot_scott_save_as);
      
      fprintf(file_scott,"reset \n");
      fprintf(file_scott,"n=%d #number of intervals\n",stats->num_bins_scott);
      fprintf(file_scott,"max=%.25f #max value\n",stats->max);
      fprintf(file_scott,"min=%.25f #min value\n",stats->min);
      fprintf(file_scott,"width=%.25f #interval width\n",stats->bin_size_scott);
      fprintf(file_scott,"#function used to map a value to the intervals\n");
      fprintf(file_scott,"hist(x,width)=width*floor(x/width)+width/2.0\n");
      fprintf(file_scott,"set term png #output terminal and file\n");
      fprintf(file_scott,"set output \"scott_histogram.png\"\n");
      fprintf(file_scott,"set xrange [min:max]\n");
      fprintf(file_scott,"set yrange [0:]\n");
      fprintf(file_scott,"#to put an empty boundary around the\n");
      fprintf(file_scott,"#data inside an autoscaled graph.\n");
      fprintf(file_scott,"set offset graph 0.05,0.05,0.05,0.0\n");
      fprintf(file_scott,"set xtics min,(max-min)/5,max\n");
      fprintf(file_scott,"set boxwidth width*0.9\n");
      fprintf(file_scott,"set style fill solid 0.5 #fillstyle\n");
      fprintf(file_scott,"set tics out nomirror\n");
      fprintf(file_scott,"set xlabel \"eta2\"\n");
      fprintf(file_scott,"set ylabel \"Frequency\"\n");
      fprintf(file_scott,"#count and plot\n");
      fprintf(file_scott,"plot \"%s\" u (hist($1,width)):(1.0) smooth freq w boxes lc rgb \"green\" notitle\n",file_name_estimator_prefix);
      fclose(file_scott);
    }
  }
  
    
  char* rank_save_as = NULL;
  D4EST_ASPRINTF(rank_save_as, "%s_%d", file_name_estimator_prefix, mpi_rank);
  FILE* file_est = fopen(rank_save_as, "w");
  free(rank_save_as);

  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        if(curved)
          fprintf(file_est, "%.20f\n", ((curved_element_data_t*)(quad->p.user_data))->local_estimator);
        else
          fprintf(file_est, "%.20f\n", ((element_data_t*)(quad->p.user_data))->local_estimator);
      }
    }
  
  fclose(file_est);
}
