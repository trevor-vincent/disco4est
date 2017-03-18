#include "sc_reduce.h"
#include "../pXest/pXest.h"
#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../dGMath/dgmath.h"
#include "../Solver/multigrid.h"
#include "../Solver/cg_eigs.h"
#include "../hpAMR/hp_amr.h"
#include "../Solver/multigrid_callbacks.h"
#include <ini.h>
#include <util.h>

static void
multigrid_update_components
(
 p4est_t* p4est,
 int level,
 problem_data_t* data
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  if (mg_data->user_callbacks != NULL){
    if (mg_data->user_callbacks->update != NULL){
      mg_data->user_callbacks->update(p4est, level, data);
    }
  }

  if (mg_data->logger != NULL){
    if (mg_data->logger->update != NULL){
      mg_data->logger->update(p4est, level, data);
    }
  }

  if (mg_data->bottom_solver != NULL){
    if (mg_data->bottom_solver->update != NULL){
      mg_data->bottom_solver->update(p4est, level, data);
    }
  }   

  if (mg_data->smoother != NULL){
    if (mg_data->smoother->update != NULL){
      mg_data->smoother->update(p4est, level, data);
    }
  }

  if (mg_data->elem_data_updater != NULL){
    if (mg_data->elem_data_updater->update != NULL){
      mg_data->elem_data_updater->update(p4est, level, data);
    }
  }
  else {
    mpi_abort("You should set elem_data_updater in multigrid\n");
  }
  
}

static
int multigrid_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  multigrid_data_t* mg_data = ((multigrid_data_t*)user);
  
  if (util_match_couple(section,"multigrid",name,"vcycle_imax")) {
    mpi_assert(mg_data->vcycle_imax == -1);
    mg_data->vcycle_imax = atoi(value);
  }
  else if (util_match_couple(section,"multigrid",name,"vcycle_rtol")) {
    mpi_assert(mg_data->vcycle_rtol == -1);
    mg_data->vcycle_rtol = atof(value);
  }
  else if (util_match_couple(section,"multigrid",name,"vcycle_atol")) {
    mpi_assert(mg_data->vcycle_atol == -1);
    mg_data->vcycle_atol = atof(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

multigrid_data_t*
multigrid_data_init
(
 p4est_t* p4est,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 int num_of_levels,
 multigrid_smoother_t* smoother,
 multigrid_bottom_solver_t* bottom_solver,
 multigrid_logger_t* logger,
 multigrid_user_callbacks_t* user_callbacks,
 multigrid_element_data_updater_t* updater,
 const char* input_file
)
{
  multigrid_data_t* mg_data = P4EST_ALLOC(multigrid_data_t, 1);
  mpi_assert(num_of_levels > 0);

  mg_data->dgmath_jit_dbase = dgmath_jit_dbase;
  mg_data->num_of_levels = num_of_levels;
  mg_data->vcycle_atol = -1;
  mg_data->vcycle_rtol = -1;
  mg_data->vcycle_imax = -1;
  mg_data->Ae_at0 = NULL;
  mg_data->err_at0 = NULL;
  mg_data->rres_at0 = NULL;
  mg_data->res_at0 = NULL;
  mg_data->elements_on_level_of_multigrid = NULL;
  mg_data->elements_on_level_of_surrogate_multigrid = NULL;
  mg_data->nodes_on_level_of_multigrid = NULL;
  mg_data->nodes_on_level_of_surrogate_multigrid = NULL;
  
  if (ini_parse(input_file, multigrid_input_handler, mg_data) < 0) {
    mpi_abort("Can't load input file");
  }
  if(mg_data->vcycle_atol == -1){
    mpi_abort("[D4EST_ERROR]: vcycle_atol not set in multigrid input");
  }
  if(mg_data->vcycle_rtol == -1){
    mpi_abort("[D4EST_ERROR]: vcycle_rtol not set in multigrid input");
  }
  if(mg_data->vcycle_imax == -1){
    mpi_abort("[D4EST_ERROR]: vcycle_imax not set in multigrid input");
  }  
  if(p4est->mpirank == 0){
    printf("[D4EST_INFO]: Multigrid Parameters\n");
    printf("[D4EST_INFO]: vcycle imax = %d\n", mg_data->vcycle_imax);
    printf("[D4EST_INFO]: vcycle rtol = %.25f\n", mg_data->vcycle_rtol);
    printf("[D4EST_INFO]: vcycle atol = %.25f\n", mg_data->vcycle_atol);
  }

  mg_data->smoother = smoother;
  mg_data->logger = logger;
  mg_data->bottom_solver = bottom_solver;
  mg_data->user_callbacks = user_callbacks;
  mg_data->elem_data_updater = updater;
  
  return mg_data;
}


void
multigrid_data_destroy(multigrid_data_t* mg_data)
{
  P4EST_FREE(mg_data);
}


static void
multigrid_vcycle
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  int level;

  /* start and end levels of multigrid */
  int toplevel = mg_data->num_of_levels - 1;
  int bottomlevel = 0;

  /* double* max_eigs = mg_data->max_eigs;/\* P4EST_ALLOC(double, mg_data->num_of_levels); *\/ */
  
  int* elements_on_level_of_multigrid = P4EST_ALLOC_ZERO(int, mg_data->num_of_levels);
  int* nodes_on_level_of_multigrid = P4EST_ALLOC_ZERO(int, mg_data->num_of_levels);
  int* elements_on_level_of_surrogate_multigrid = P4EST_ALLOC_ZERO(int, mg_data->num_of_levels);
  int* nodes_on_level_of_surrogate_multigrid = P4EST_ALLOC_ZERO(int, mg_data->num_of_levels);

  int total_elements_on_surrogate_multigrid = 0;
  int total_nodes_on_surrogate_multigrid = 0;
  int total_elements_on_multigrid = 0;
  int total_nodes_on_multigrid = 0;

  elements_on_level_of_multigrid[toplevel] = p4est->local_num_quadrants;
  nodes_on_level_of_multigrid[toplevel] = vecs->local_nodes;
  
  total_elements_on_multigrid += elements_on_level_of_multigrid[toplevel];
  total_nodes_on_multigrid += nodes_on_level_of_multigrid[toplevel];
  
  elements_on_level_of_surrogate_multigrid[toplevel] = p4est->local_num_quadrants;
  nodes_on_level_of_surrogate_multigrid[toplevel] = nodes_on_level_of_multigrid[toplevel];
  
  total_elements_on_surrogate_multigrid += elements_on_level_of_surrogate_multigrid[toplevel];
  total_nodes_on_surrogate_multigrid += nodes_on_level_of_surrogate_multigrid[toplevel];


  
  /* FIXME: REALLOCATE AS WE COARSEN */
  multigrid_refine_data_t* coarse_grid_refinement = P4EST_ALLOC
                                                    (
                                                     multigrid_refine_data_t,
                                                     (toplevel-bottomlevel) * p4est->local_num_quadrants
                                                    );
  mg_data->coarse_grid_refinement = coarse_grid_refinement;
  
  /* initialize */
  mg_data->fine_nodes = total_nodes_on_multigrid;
  mg_data->coarse_nodes = total_nodes_on_multigrid;

  double* Ae_at0 = P4EST_ALLOC(double, total_nodes_on_multigrid);
  double* err_at0 = P4EST_ALLOC_ZERO(double, total_nodes_on_multigrid);
  double* res_at0 = P4EST_ALLOC(double, total_nodes_on_multigrid);
  double* rres_at0 = P4EST_ALLOC(double, total_nodes_on_multigrid);

  /* for logging */
  mg_data->Ae_at0 = &(Ae_at0)[0];
  mg_data->err_at0 = &(err_at0)[0];
  mg_data->res_at0 = &(res_at0)[0];
  mg_data->rres_at0 = &(rres_at0)[0];
  mg_data->elements_on_level_of_multigrid = elements_on_level_of_multigrid;
  mg_data->nodes_on_level_of_multigrid = nodes_on_level_of_multigrid;
  mg_data->elements_on_level_of_surrogate_multigrid = elements_on_level_of_surrogate_multigrid;
  mg_data->nodes_on_level_of_surrogate_multigrid = nodes_on_level_of_surrogate_multigrid;


  
  double* u = vecs->u;
  
  problem_data_t vecs_for_smooth;
  problem_data_copy_ptrs(vecs, &vecs_for_smooth);
  problem_data_t vecs_for_bottom_solve;
  problem_data_copy_ptrs(vecs, &vecs_for_bottom_solve);
  
  /* initialize error to zero */
  /* linalg_fill_vec(mg_data->err, 0., mg_data->fine_nodes); */

  /* initialize stride */
  mg_data->stride = 0;


  /**********************************************************/
  /**********************************************************/
  /******************* BEGIN GOING DOWN V *******************/
  /**********************************************************/
  /**********************************************************/  
  /* DEBUG_PRINT_ARR_DBL(vecs->u, vecs->local_nodes); */

  #ifdef DEBUG_INFO
  printf("Level = %d\n", toplevel);
  printf("State = %s\n", "PRE_V");
  printf("elements_on_level = %d\n", elements_on_level_of_multigrid[toplevel]);
  printf("elements_on_surrogate_level = %d\n", elements_on_level_of_surrogate_multigrid[toplevel]);
  printf("nodes_on_level = %d\n", nodes_on_level_of_multigrid[toplevel]);
  printf("nodes_on_surrogate_level = %d\n", nodes_on_level_of_surrogate_multigrid[toplevel]);
#endif    
  mg_data->mg_state = PRE_V; multigrid_update_components(p4est, toplevel, NULL);

  
  int stride_to_fine_data = 0;
  for (level = toplevel; level > bottomlevel; --level){

    int fine_level = level;
    int coarse_level = level-1;

    /**********************************************************/
    /**********************************************************/
    /******************* PRE SMOOTH ***************************/
    /**********************************************************/
    /**********************************************************/  
    
    /* set initial guess for error */
    linalg_fill_vec(&err_at0[stride_to_fine_data], 0., nodes_on_level_of_multigrid[level]);//mg_data->fine_nodes);

    if (level == toplevel){
      vecs_for_smooth.Au = vecs->Au;
      vecs_for_smooth.u = vecs->u; 
      vecs_for_smooth.rhs = vecs->rhs;
      vecs_for_smooth.local_nodes = vecs->local_nodes;
    }
    else{
      vecs_for_smooth.Au = &Ae_at0[stride_to_fine_data];//mg_data->Ae;//[mg_data->fine_nodes];
      vecs_for_smooth.u = &err_at0[stride_to_fine_data];//mg_data->err;//[mg_data->fine_nodes];
      vecs_for_smooth.rhs = &res_at0[stride_to_fine_data];//mg_data->res;//)[mg_data->fine_nodes];
      vecs_for_smooth.local_nodes = nodes_on_level_of_multigrid[level];
    }

    /**********************************************************/
    /**********************************************************/
    /******************* BEGIN SMOOTH *************************/
    /**********************************************************/
    /**********************************************************/  

#ifdef DEBUG_INFO
    printf("Level = %d\n", level);
    printf("State = %s\n", "DOWNV_PRE_SMOOTH");
    printf("elements_on_level = %d\n", elements_on_level_of_multigrid[level]);
    printf("elements_on_surrogate_level = %d\n", elements_on_level_of_surrogate_multigrid[level]);
    printf("nodes_on_level = %d\n", nodes_on_level_of_multigrid[level]);
    printf("nodes_on_surrogate_level = %d\n", nodes_on_level_of_surrogate_multigrid[level]);
#endif
    mg_data->mg_state = DOWNV_PRE_SMOOTH; multigrid_update_components(p4est, level, &vecs_for_smooth);
    
    mg_data->smoother->smooth
      (
       p4est,
       &vecs_for_smooth,
       fcns,
       &rres_at0[stride_to_fine_data],
       fine_level
      );
    
#ifdef DEBUG_INFO
    printf("Level = %d\n", level);
    printf("State = %s\n", "DOWNV_POST_SMOOTH");
    printf("elements_on_level = %d\n", elements_on_level_of_multigrid[level]);
    printf("elements_on_surrogate_level = %d\n", elements_on_level_of_surrogate_multigrid[level]);
    printf("nodes_on_level = %d\n", nodes_on_level_of_multigrid[level]);
    printf("nodes_on_surrogate_level = %d\n", nodes_on_level_of_surrogate_multigrid[level]);
#endif    
    mg_data->mg_state = DOWNV_POST_SMOOTH; multigrid_update_components(p4est, level, &vecs_for_smooth);
    
    /**********************************************************/
    /**********************************************************/
    /********************* END SMOOTH *************************/
    /**********************************************************/
    /**********************************************************/  


    /**********************************************************/
    /**********************************************************/
    /******************* BEGIN COARSEN ************************/
    /**********************************************************/
    /**********************************************************/  
#ifdef DEBUG_INFO
    printf("Level = %d\n", level);
    printf("State = %s\n", "DOWNV_PRE_COARSEN");
    printf("elements_on_level = %d\n", elements_on_level_of_multigrid[level]);
    printf("elements_on_surrogate_level = %d\n", elements_on_level_of_surrogate_multigrid[level]);
    printf("nodes_on_level = %d\n", nodes_on_level_of_multigrid[level]);
    printf("nodes_on_surrogate_level = %d\n", nodes_on_level_of_surrogate_multigrid[level]);
#endif
    mg_data->mg_state = DOWNV_PRE_COARSEN; multigrid_update_components(p4est, level, NULL);
    
    /* increments the stride */
    p4est_coarsen_ext (p4est,
                       0,
                       1,
                       multigrid_coarsen,
                       multigrid_coarsen_init,
                       NULL
                      );

#ifdef DEBUG_INFO
    printf("Level = %d\n", level);
    printf("State = %s\n", "DOWNV_POST_COARSEN");
    printf("elements_on_level = %d\n", elements_on_level_of_multigrid[level]);
    printf("elements_on_surrogate_level = %d\n", elements_on_level_of_surrogate_multigrid[level]);
    printf("nodes_on_level = %d\n", nodes_on_level_of_multigrid[level]);
    printf("nodes_on_surrogate_level = %d\n", nodes_on_level_of_surrogate_multigrid[level]);
#endif
    
    mg_data->mg_state = DOWNV_POST_COARSEN; multigrid_update_components(p4est, level, NULL);

    /**********************************************************/
    /**********************************************************/
    /******************* END COARSEN **************************/
    /**********************************************************/
    /**********************************************************/  
    
    /* p4est has changed, so update the element data, this does not change the stride */
    /* element_data_init(p4est, -1); */

    /* update surrogate info */
    nodes_on_level_of_surrogate_multigrid[level-1] = mg_data->elem_data_updater->get_local_nodes(p4est);
    elements_on_level_of_surrogate_multigrid[level-1] = p4est->local_num_quadrants;

    total_elements_on_surrogate_multigrid += elements_on_level_of_surrogate_multigrid[level-1];
    total_nodes_on_surrogate_multigrid += nodes_on_level_of_surrogate_multigrid[level-1];
        
    /* decrements the stride */
    mg_data->stride -= elements_on_level_of_surrogate_multigrid[level-1];

    /**********************************************************/
    /**********************************************************/
    /******************* BEGIN BALANCE ************************/
    /**********************************************************/
    /**********************************************************/  
    mg_data->mg_state = DOWNV_PRE_BALANCE; multigrid_update_components(p4est, level, NULL);   
    
    /* does not change the stride */
    p4est_balance_ext (p4est, P4EST_CONNECT_FACE, NULL, multigrid_store_balance_changes);

    mg_data->mg_state = DOWNV_POST_BALANCE; multigrid_update_components(p4est, level, NULL);
    /**********************************************************/
    /**********************************************************/
    /******************* END BALANCE **************************/
    /**********************************************************/
    /**********************************************************/  
    
    /* element_data_init(p4est, -1); */

    /* /\* update ghost data *\/ */
    /* P4EST_FREE(*ghost_data); */
    /* p4est_ghost_destroy (*ghost); */
    /* *ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE); */
    /* *ghost_data = P4EST_ALLOC (element_data_t, */
    /*                            (*ghost)->ghosts.elem_count); */

    /* update coarse grid info */
    elements_on_level_of_multigrid[level-1] = p4est->local_num_quadrants;
    total_elements_on_multigrid += elements_on_level_of_multigrid[level-1];
    nodes_on_level_of_multigrid[level-1] = mg_data->elem_data_updater->get_local_nodes(p4est);
    total_nodes_on_multigrid += nodes_on_level_of_multigrid[level-1];

    mg_data->fine_nodes = nodes_on_level_of_multigrid[level];
    mg_data->coarse_nodes = nodes_on_level_of_multigrid[level - 1];

    Ae_at0 = P4EST_REALLOC(Ae_at0, double, total_nodes_on_multigrid);
    err_at0 = P4EST_REALLOC(err_at0, double, total_nodes_on_multigrid);
    res_at0 = P4EST_REALLOC(res_at0, double, total_nodes_on_multigrid);
    rres_at0 = P4EST_REALLOC(rres_at0, double, total_nodes_on_multigrid);


    /* always zero these before restriction or prolongation */
    mg_data->coarse_stride = 0;
    mg_data->fine_stride = 0;
    mg_data->temp_stride = 0;
    mg_data->intergrid_ptr = &rres_at0[stride_to_fine_data];//&(mg_data->rres)[0];

#ifdef DEBUG_INFO
    printf("Level = %d\n", level);
    printf("State = %s\n", "DOWNV_PRE_RESTRICTION");
    printf("elements_on_level = %d\n", elements_on_level_of_multigrid[level]);
    printf("elements_on_surrogate_level = %d\n", elements_on_level_of_surrogate_multigrid[level]);
    printf("nodes_on_level = %d\n", nodes_on_level_of_multigrid[level]);
    printf("nodes_on_surrogate_level = %d\n", nodes_on_level_of_surrogate_multigrid[level]);
#endif
    
    mg_data->mg_state = DOWNV_PRE_RESTRICTION; multigrid_update_components(p4est, level, NULL);

    
    
    p4est_iterate(
                  p4est,
                  NULL,
                  NULL,
                  multigrid_apply_restriction,
                  NULL,
#if (P4EST_DIM)==3
                  NULL,
#endif
                  NULL
    );

#ifdef DEBUG_INFO
    printf("Level = %d\n", level);
    printf("State = %s\n", "DOWNV_POST_RESTRICTION");
    printf("elements_on_level = %d\n", elements_on_level_of_multigrid[level]);
    printf("elements_on_surrogate_level = %d\n", elements_on_level_of_surrogate_multigrid[level]);
    printf("nodes_on_level = %d\n", nodes_on_level_of_multigrid[level]);
    printf("nodes_on_surrogate_level = %d\n", nodes_on_level_of_surrogate_multigrid[level]);    
#endif
    
    mg_data->mg_state = DOWNV_POST_RESTRICTION; multigrid_update_components(p4est, level, NULL);   

    linalg_copy_1st_to_2nd
      (
       &(rres_at0)[stride_to_fine_data + mg_data->fine_nodes],
       &(res_at0)[stride_to_fine_data + mg_data->fine_nodes],
       nodes_on_level_of_multigrid[coarse_level]//mg_data->coarse_nodes
      );

    stride_to_fine_data += mg_data->fine_nodes;
  }

  /**********************************************************/
  /**********************************************************/
  /******************* END GOING DOWN V *********************/
  /**********************************************************/
  /**********************************************************/

  

  /**********************************************************/
  /**********************************************************/
  /******************* BEGIN COARSE SOLVE *******************/
  /**********************************************************/
  /**********************************************************/

  /* set initial guess for error */
  linalg_fill_vec(&err_at0[stride_to_fine_data], 0., nodes_on_level_of_multigrid[bottomlevel]);
  
  /* vecs_for_bottom_solve.Au = mg_data->Ae;//[mg_data->fine_nodes]; */
  vecs_for_bottom_solve.Au = &Ae_at0[stride_to_fine_data];//[mg_data->fine_nodes];
  vecs_for_bottom_solve.u = &err_at0[stride_to_fine_data];//mg_data->err;//[mg_data->fine_nodes];
  vecs_for_bottom_solve.rhs = &res_at0[stride_to_fine_data];//mg_data->res;//)[mg_data->fine_nodes];
  vecs_for_bottom_solve.local_nodes = nodes_on_level_of_multigrid[bottomlevel];


  mg_data->mg_state = COARSE_PRE_SOLVE; multigrid_update_components(p4est, level, &vecs_for_bottom_solve);   
  
  mg_data->bottom_solver->solve
    (
     p4est,
     &vecs_for_bottom_solve,
     fcns,
     /* &(mg_data->rres[mg_data->fine_nodes]), */
     &rres_at0[stride_to_fine_data]//[mg_data->fine_nodes]),
    );

  mg_data->mg_state = COARSE_POST_SOLVE; multigrid_update_components(p4est, level, &vecs_for_bottom_solve);   

  /* util_print_matrix(&rres_at0[stride_to_fine_data],mg_data->coarse_nodes,1,"rres coarse solve= ", 0); */
  
  /**********************************************************/
  /**********************************************************/
  /******************* END COARSE SOLVE *********************/
  /**********************************************************/
  /**********************************************************/
    
  mg_data->stride = total_elements_on_surrogate_multigrid;
  mg_data->stride -= elements_on_level_of_surrogate_multigrid[toplevel]; /* subtract fine grid */
  
  /**********************************************************/
  /**********************************************************/
  /******************* BEGIN GOING UP THE V *****************/
  /**********************************************************/
  /**********************************************************/
  for (level = bottomlevel; level < toplevel; ++level){

    int fine_level = level + 1;
    int coarse_level = level;

    mg_data->fine_nodes = nodes_on_level_of_multigrid[level+1];
    mg_data->coarse_nodes = nodes_on_level_of_multigrid[level];
    stride_to_fine_data -= mg_data->fine_nodes;
    mg_data->stride -= elements_on_level_of_surrogate_multigrid[level];
    mg_data->fine_stride = 0;
    mg_data->coarse_stride = 0;
    mg_data->temp_stride = 0;

    
    linalg_copy_1st_to_2nd(&err_at0[stride_to_fine_data + mg_data->fine_nodes],//&(mg_data->err)[mg_data->fine_nodes],
                           &rres_at0[stride_to_fine_data + mg_data->fine_nodes],//&(mg_data->rres)[mg_data->fine_nodes],
                           nodes_on_level_of_multigrid[coarse_level]);

    mg_data->intergrid_ptr = &rres_at0[stride_to_fine_data];/* mg_data->rres */;


    /**********************************************************/
    /**********************************************************/
    /******************* BEGIN REFINE/PROLONGATION ************/
    /**********************************************************/
    /**********************************************************/

    
    mg_data->mg_state = UPV_PRE_REFINE; multigrid_update_components(p4est, level, NULL);   
    
    /* increments the stride */
    p4est_refine_ext(p4est,
                     0,
                     P4EST_QMAXLEVEL,
                     multigrid_refine_and_apply_prolongation,
                     NULL,
                     multigrid_refine_and_apply_prolongation_replace
                    );


    mg_data->mg_state = UPV_POST_REFINE; multigrid_update_components(p4est, level, NULL);   


    /**********************************************************/
    /**********************************************************/
    /******************* END REFINE/PROLONGATION **************/
    /**********************************************************/
    /**********************************************************/


    
    /* P4EST_FREE(*ghost_data); */
    /* p4est_ghost_destroy (*ghost); */
    /* *ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE); */
    /* *ghost_data = P4EST_ALLOC (element_data_t, */
    /*                            (*ghost)->ghosts.elem_count); */

    /* element_data_init(p4est, -1); */
    
    /* decrements the stride */
    mg_data->stride -= elements_on_level_of_surrogate_multigrid[level];

    /* If the grid is not the finiest */
    if(fine_level != toplevel){
      vecs_for_smooth.Au = &Ae_at0[stride_to_fine_data];//mg_data->Ae;
      vecs_for_smooth.u = &err_at0[stride_to_fine_data];//mg_data->err;
      vecs_for_smooth.rhs = &res_at0[stride_to_fine_data];//mg_data->res;
      vecs_for_smooth.local_nodes = mg_data->fine_nodes;
    }
    else {
      vecs_for_smooth.Au = vecs->Au;//&Ae_at0[stride_to_fine_data];//mg_data->Ae;
      vecs_for_smooth.u = vecs->u;//&err_at0[stride_to_fine_data];//mg_data->err;
      vecs_for_smooth.rhs = vecs->rhs;//&res_at0[stride_to_fine_data];//mg_data->res;
      vecs_for_smooth.local_nodes = mg_data->fine_nodes;
    }
     
    linalg_vec_axpy(1.0, &rres_at0[stride_to_fine_data], vecs_for_smooth.u, mg_data->fine_nodes);



    /**********************************************************/
    /**********************************************************/
    /******************* BEGIN SMOOTH *************************/
    /**********************************************************/
    /**********************************************************/
    

    mg_data->mg_state = UPV_PRE_SMOOTH; multigrid_update_components(p4est, level, &vecs_for_smooth);   
    
    mg_data->smoother->smooth
      (
       p4est,
       &vecs_for_smooth,
       fcns,
       &rres_at0[stride_to_fine_data],//mg_data->rres,
       fine_level
      );

    mg_data->mg_state = UPV_POST_SMOOTH; multigrid_update_components(p4est, level, &vecs_for_smooth);   


    /**********************************************************/
    /**********************************************************/
    /******************* END SMOOTH ***************************/
    /**********************************************************/
    /**********************************************************/

    
  }

  mg_data->mg_state = POST_V; multigrid_update_components(p4est, toplevel, NULL);
  
  for (level = toplevel; level >= bottomlevel; --level){
    printf("Level %d, Number of elements %d, Number of nodes %d\n",
           level,
           elements_on_level_of_multigrid[level],
           nodes_on_level_of_multigrid[level]
          );
  }


  
  mg_data->vcycle_r2_local_current = linalg_vec_dot(&rres_at0[stride_to_fine_data],
                                           &rres_at0[stride_to_fine_data],
                                           mg_data->fine_nodes);
  
  /**********************************************************/
  /**********************************************************/
  /******************* END GOING UP THE V *******************/
  /**********************************************************/
  /**********************************************************/
  
  P4EST_FREE(elements_on_level_of_multigrid);
  P4EST_FREE(nodes_on_level_of_multigrid);
  P4EST_FREE(elements_on_level_of_surrogate_multigrid);
  P4EST_FREE(nodes_on_level_of_surrogate_multigrid);
  P4EST_FREE(coarse_grid_refinement);
  P4EST_FREE(Ae_at0);
  P4EST_FREE(err_at0);
  P4EST_FREE(res_at0);
  P4EST_FREE(rres_at0);
}


static double
multigrid_compute_residual
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns
){
  
  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  p4est_ghost_t* ghost = *(updater->ghost);
  void* ghost_data = *(updater->ghost_data);
  d4est_geometry_t* d4est_geom = updater->d4est_geom;
  
  if (mg_data->mg_state == START){
    double* Au; 
    double* rhs;
    int local_nodes = vecs->local_nodes;
    Au = vecs->Au;
    rhs = vecs->rhs;
    double* r = P4EST_ALLOC(double, vecs->local_nodes);
    fcns->apply_lhs(p4est, ghost, ghost_data, vecs, mg_data->dgmath_jit_dbase, d4est_geom);  
    linalg_vec_axpyeqz(-1., Au, rhs, r, local_nodes);
    double r2_0_local = linalg_vec_dot(r,r,local_nodes);  
    P4EST_FREE(r);
  
    double r2_0_global = -1;
    sc_allreduce
      (
       &r2_0_local,
       &r2_0_global,
       1,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       sc_MPI_COMM_WORLD
      );  
  
    return r2_0_global;
  }
  else {
    double r2_i_local = mg_data->vcycle_r2_local_current;
    double r2_i_global;
    sc_allreduce
      (
       &r2_i_local,
       &r2_i_global,
       1,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       sc_MPI_COMM_WORLD
      );
    return r2_i_global;
  }
}

void
multigrid_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 multigrid_data_t* mg_data
)
{
  void* tmp_ptr = p4est->user_pointer;
  p4est->user_pointer = mg_data;
  /*
   * Calculate the initial residual
   * for termination condition
   */
  
  mg_data->mg_state = START;
  mg_data->vcycle_r2_global_current = multigrid_compute_residual
                                               (
                                                p4est,
                                                vecs,
                                                fcns
                                               );

  mg_data->vcycle_r2_global_last = mg_data->vcycle_r2_global_current;
  
  mg_data->vcycle_num_finished = 0;
  mg_data->vcycle_r2_global_stoptol = mg_data->vcycle_rtol*mg_data->vcycle_rtol*mg_data->vcycle_r2_global_current
                                      + mg_data->vcycle_atol*mg_data->vcycle_atol;


  multigrid_update_components(p4est, -1, NULL);
  /**
   * START VCYCLING
   * 
   */
  while
    (
     mg_data->vcycle_num_finished <= mg_data->vcycle_imax
     &&
     mg_data->vcycle_r2_global_current > mg_data->vcycle_r2_global_stoptol
    )
    {

      multigrid_vcycle(p4est, vecs, fcns);

      mg_data->vcycle_r2_global_current = multigrid_compute_residual
                                          (
                                           p4est,
                                           vecs,
                                           fcns
                                          );

      mg_data->vcycle_num_finished++;
      mg_data->mg_state = POST_RESIDUAL_UPDATE;
      multigrid_update_components(p4est, -1, NULL);
  
      if (sqrt(mg_data->vcycle_r2_global_current
               /mg_data->vcycle_r2_global_last) >= .99){
        break;
      }
      
      mg_data->vcycle_r2_global_last = mg_data->vcycle_r2_global_current;
    }


  mg_data->mg_state = END;
  multigrid_update_components(p4est, -1, NULL);
  
  p4est->user_pointer = tmp_ptr;
}
