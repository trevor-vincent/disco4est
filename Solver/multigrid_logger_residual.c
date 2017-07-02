#include <multigrid_logger_residual.h>

void
multigrid_logger_residual_destroy
(
 multigrid_logger_t* logger
)
{
  P4EST_FREE(logger);
}

static
void
multigrid_logger_residual_update
(
 p4est_t* p4est,
 int level,
 d4est_elliptic_problem_data_t* vecs
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  double r2_i_global = mg_data->vcycle_r2_global_current;
  double old_r2_i_global = mg_data->vcycle_r2_global_last;
  int v = mg_data->vcycle_num_finished;
  
  if(mg_data->mg_state == START || mg_data->mg_state == POST_RESIDUAL_UPDATE){
    printf("[MG_SOLVER]: %d %.30f %f\n", v, sqrt(r2_i_global), sqrt(r2_i_global/old_r2_i_global));
  }  
  else {
    return;
  }
}

multigrid_logger_t*
multigrid_logger_residual_init
(
)
{
  multigrid_logger_t* logger = P4EST_ALLOC(multigrid_logger_t, 1);
  logger->update = multigrid_logger_residual_update;
  logger->user = NULL;
  return logger;
}
