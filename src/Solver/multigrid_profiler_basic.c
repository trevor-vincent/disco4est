#include <multigrid_profiler_basic.h>
#include <time.h>
#include <sc_reduce.h>

void
multigrid_profiler_basic_destroy
(
 multigrid_profiler_t* profiler
)
{
  P4EST_FREE(profiler->user);
  P4EST_FREE(profiler);
}

static
void
multigrid_profiler_basic_print(multigrid_profiler_t* profiler, int mpi_rank)
{
  
 multigrid_profiler_basic_data_t* profiler_data = profiler->user;

 profiler_data->time_total = profiler_data->time_start_pre_v +
                        profiler_data->time_start_pre_v +
                        profiler_data->time_downv_smooth +
                        profiler_data->time_downv_coarsen +
                        profiler_data->time_downv_balance +
                        profiler_data->time_downv_restriction +
                        profiler_data->time_coarse_solve +
                        profiler_data->time_upv_refine +
                        profiler_data->time_upv_smooth +
                        profiler_data->time_post_v_end +
                        profiler_data->time_total;
 
 double times_local [] =
   {
     profiler_data->time_start_pre_v,
     profiler_data->time_downv_smooth,
     profiler_data->time_downv_coarsen,
     profiler_data->time_downv_balance,
     profiler_data->time_downv_restriction,
     profiler_data->time_coarse_solve,
     profiler_data->time_upv_refine,
     profiler_data->time_upv_smooth,
     profiler_data->time_post_v_end,
     profiler_data->time_total
   };

 double times_global [10];
  
 sc_reduce
   (
    &times_local[0],
    &times_global[0],
    10,
    sc_MPI_DOUBLE,
    sc_MPI_SUM,
    0,
    sc_MPI_COMM_WORLD
   );    


 zlog_category_t *c_default = zlog_get_category("d4est_multigrid_profiler");
 if (mpi_rank == 0){
   zlog_info(c_default, "Multigrid profiler report");
   zlog_info(c_default, "*************************");
   zlog_info(c_default, "time_start_pre_v = %f, fraction = %f percent", times_global[0], 100*times_global[0]/times_global[9]);
   zlog_info(c_default, "time_downv_smooth = %f, fraction = %f percent", times_global[1], 100*times_global[1]/times_global[9]);
   zlog_info(c_default, "time_downv_coarsen = %f, fraction = %f percent", times_global[2], 100*times_global[2]/times_global[9]);
   zlog_info(c_default, "time_downv_balance = %f, fraction = %f percent", times_global[3], 100*times_global[3]/times_global[9]);
   zlog_info(c_default, "time_downv_restriction = %f, fraction = %f percent", times_global[4], 100*times_global[4]/times_global[9]);
   zlog_info(c_default, "time_coarse_solve = %f, fraction = %f percent", times_global[5], 100*times_global[5]/times_global[9]);
   zlog_info(c_default, "time_upv_refine = %f, fraction = %f percent", times_global[6], 100*times_global[6]/times_global[9]);
   zlog_info(c_default, "time_upv_smooth = %f, fraction = %f percent", times_global[7], 100*times_global[7]/times_global[9]);
   zlog_info(c_default, "time_post_v_end = %f, fraction = %f percent", times_global[8], 100*times_global[8]/times_global[9]);
   zlog_info(c_default, "time_total = %f, fraction = %f percent", times_global[9], 100*times_global[9]/times_global[9]);
   zlog_info(c_default, "*************************");
 }
}


static
void
multigrid_profiler_basic_update
(
 p4est_t* p4est,
 int level,
 d4est_elliptic_data_t* vecs
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_profiler_t* profiler = mg_data->profiler;
  multigrid_profiler_basic_data_t* profiler_data = profiler->user;
  
  if(mg_data->mg_state == START){
    profiler_data->time_start_pre_v = 0.;
    profiler_data->time_downv_smooth = 0.;
    profiler_data->time_downv_coarsen = 0.;
    profiler_data->time_downv_balance = 0.;
    profiler_data->time_downv_restriction = 0.;
    profiler_data->time_coarse_solve = 0.;
    profiler_data->time_upv_refine = 0.;
    profiler_data->time_upv_smooth = 0.;
    profiler_data->time_post_v_end = 0.;
    profiler_data->time_total = 0.;    
    profiler_data->begin = clock();
  }
  else if (mg_data->mg_state == PRE_V){
    profiler_data->end = clock();
    profiler_data->time_start_pre_v += (double)(profiler_data->end - profiler_data->begin)/(CLOCKS_PER_SEC);
  }
  else if (mg_data->mg_state == DOWNV_PRE_SMOOTH){
    profiler_data->begin = clock();    
  }
  else if (mg_data->mg_state == DOWNV_POST_SMOOTH){
    profiler_data->end = clock();
    profiler_data->time_downv_smooth += (double)(profiler_data->end - profiler_data->begin)/(CLOCKS_PER_SEC);
  }
  else if (mg_data->mg_state == DOWNV_PRE_COARSEN){
    profiler_data->begin = clock();
  }
  else if (mg_data->mg_state == DOWNV_POST_COARSEN){
    profiler_data->end = clock();
    profiler_data->time_downv_coarsen += (double)(profiler_data->end - profiler_data->begin)/(CLOCKS_PER_SEC);
  }  
  else if (mg_data->mg_state == DOWNV_PRE_BALANCE){
    profiler_data->begin = clock();
  }
  else if (mg_data->mg_state == DOWNV_POST_BALANCE){
    profiler_data->end = clock();
    profiler_data->time_downv_balance += (double)(profiler_data->end - profiler_data->begin)/(CLOCKS_PER_SEC);
  }
  else if (mg_data->mg_state == DOWNV_PRE_RESTRICTION){
    profiler_data->begin = clock();
  }
  else if (mg_data->mg_state == DOWNV_POST_RESTRICTION){
    profiler_data->end = clock();
    profiler_data->time_downv_restriction += (double)(profiler_data->end - profiler_data->begin)/(CLOCKS_PER_SEC);
  }    
  else if (mg_data->mg_state == COARSE_PRE_SOLVE){
    profiler_data->begin = clock();
  }
  else if (mg_data->mg_state == COARSE_POST_SOLVE){
    profiler_data->end = clock();
    profiler_data->time_coarse_solve += (double)(profiler_data->end - profiler_data->begin)/(CLOCKS_PER_SEC);    
  }
  else if (mg_data->mg_state == UPV_PRE_REFINE){
    profiler_data->begin = clock();
  }
  else if (mg_data->mg_state == UPV_POST_REFINE){
    profiler_data->end = clock();
    profiler_data->time_upv_refine += (double)(profiler_data->end - profiler_data->begin)/(CLOCKS_PER_SEC);        
  }
  else if (mg_data->mg_state == UPV_PRE_SMOOTH){
    profiler_data->begin = clock();
  }
  else if (mg_data->mg_state == UPV_POST_SMOOTH){
    profiler_data->end = clock();
    profiler_data->time_upv_smooth += (double)(profiler_data->end - profiler_data->begin)/(CLOCKS_PER_SEC);
    profiler_data->begin = clock();
  }
  else if (mg_data->mg_state == POST_V){
    profiler_data->begin = clock();
  }
  else if (mg_data->mg_state == END){
    profiler_data->end = clock();
    profiler_data->time_post_v_end += (double)(profiler_data->end - profiler_data->begin)/(CLOCKS_PER_SEC);
    multigrid_profiler_basic_print(profiler, p4est->mpirank);
  }  
  else {
    return;
  }
}

multigrid_profiler_t*
multigrid_profiler_basic_init
(
)
{
  multigrid_profiler_t* profiler = P4EST_ALLOC(multigrid_profiler_t, 1);
  multigrid_profiler_basic_data_t* profiler_data = P4EST_ALLOC(multigrid_profiler_basic_data_t, 1);
  profiler->update = multigrid_profiler_basic_update;
  profiler->user = profiler_data;
  return profiler;
}
