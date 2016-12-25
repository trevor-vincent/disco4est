#include "sc_reduce.h"
#include "../pXest/pXest.h"
#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../dGMath/dgmath.h"
#include "../Solver/multigrid.h"
#include "../Solver/cg_eigs.h"
#include "../hpAMR/hp_amr.h"
#include "../Solver/multigrid_cg_coarse_solver.h"
#include "../Solver/multigrid_cheby_smoother.h"
#include "../Solver/multigrid_callbacks.h"

multigrid_data_t*
multigrid_data_init
(
 int mpi_rank, 
 int num_of_levels,
 int vcycle_iter,
 /* Stopping condition := r2 <= vcycle_rtol*vcycle_rtol*r2initial + vcycle_atol*vcycle_atol */
 double vcycle_rtol,
 double vcycle_atol,
 int smooth_iter,
 /* Iterations of cg to find spectral radius */
 int cg_eigs_iter,
 /* Multiplier for maximum eigenvalue estimate, = 1.0 for no multiplier */
 double max_eig_factor,
 /* Reuse eigenvalue estimates from the last vcycle, =1 is most sensible */
 int max_eig_reuse,
 double lmax_lmin_rat,
 multigrid_coarse_solver_t coarse_solver_type,
 int coarse_iter,
 double coarse_rtol,
 int save_vtk_snapshot,
 int perform_checksum,
 multigrid_log_option_t log_option,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  multigrid_data_t* mg_data = P4EST_ALLOC(multigrid_data_t, 1);
  mg_data->max_eigs = P4EST_ALLOC(double, num_of_levels);
  
  mg_data->user_defined_fields = 0;
  mg_data->mg_prolong_user_callback = NULL;
  mg_data->mg_restrict_user_callback = NULL;
  mg_data->mg_update_user_callback = NULL;
  mg_data->user_ctx = NULL;
  
  mg_data->mpi_rank = mpi_rank;
  mg_data->num_of_levels = num_of_levels;
  mg_data->vcycle_iter = vcycle_iter;
  mg_data->vcycle_rtol = vcycle_rtol;
  mg_data->vcycle_atol = vcycle_atol;
  mg_data->smooth_iter = smooth_iter;
  mg_data->cg_eigs_iter = cg_eigs_iter;
  mg_data->max_eig_factor = max_eig_factor;
  mg_data->max_eig_reuse = max_eig_reuse;
  mg_data->lmax_lmin_rat = lmax_lmin_rat;
  mg_data->coarse_solver_type = coarse_solver_type;
  mg_data->coarse_iter = coarse_iter;
  mg_data->coarse_rtol = coarse_rtol;
  mg_data->save_vtk_snapshot = save_vtk_snapshot;
  mg_data->perform_checksum = perform_checksum;
  mg_data->log_option = log_option;
  mg_data->dgmath_jit_dbase = dgmath_jit_dbase;

  return mg_data;
}

void
multigrid_data_set_user_defined_fields
(
 multigrid_data_t* mg_data,
 multigrid_prolong_user_callback_fcn_t mg_prolong_user_callback,
 multigrid_restrict_user_callback_fcn_t mg_restrict_user_callback,
 multigrid_update_user_callback_fcn_t mg_update_user_callback,
 void* user_ctx
)
{
  mg_data->user_defined_fields = 1;
  mg_data->mg_prolong_user_callback = mg_prolong_user_callback;
  mg_data->mg_restrict_user_callback = mg_restrict_user_callback;
  mg_data->mg_update_user_callback = mg_update_user_callback;
  mg_data->user_ctx = user_ctx;
}

void
multigrid_data_destroy(multigrid_data_t* mg_data)
{
  P4EST_FREE(mg_data->max_eigs);
  P4EST_FREE(mg_data);
}

static void
multigrid_vcycle
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 multigrid_data_t* mg_data,
 p4est_ghost_t** ghost,
 element_data_t** ghost_data
)
{
  int level;
  void* tmpptr = p4est->user_pointer;
  p4est->user_pointer = (void*) mg_data;

  /* start and end levels of multigrid */
  int endlevel = mg_data->num_of_levels - 1;
  int startlevel = 0;

  double* max_eigs = mg_data->max_eigs;/* P4EST_ALLOC(double, mg_data->num_of_levels); */
  
  int* elements_on_level_of_multigrid = P4EST_ALLOC_ZERO(int, mg_data->num_of_levels);
  int* nodes_on_level_of_multigrid = P4EST_ALLOC_ZERO(int, mg_data->num_of_levels);
  int* elements_on_level_of_surrogate_multigrid = P4EST_ALLOC_ZERO(int, mg_data->num_of_levels);
  int* nodes_on_level_of_surrogate_multigrid = P4EST_ALLOC_ZERO(int, mg_data->num_of_levels);

  int total_elements_on_surrogate_multigrid = 0;
  int total_nodes_on_surrogate_multigrid = 0;
  int total_elements_on_multigrid = 0;
  int total_nodes_on_multigrid = 0;
  
  elements_on_level_of_multigrid[0] = p4est->local_num_quadrants;
  nodes_on_level_of_multigrid[0] = vecs->local_nodes;
  
  total_elements_on_multigrid += elements_on_level_of_multigrid[0];
  total_nodes_on_multigrid += nodes_on_level_of_multigrid[0];
  
  elements_on_level_of_surrogate_multigrid[0] = p4est->local_num_quadrants;
  nodes_on_level_of_surrogate_multigrid[0] = nodes_on_level_of_multigrid[0];
  
  total_elements_on_surrogate_multigrid += elements_on_level_of_surrogate_multigrid[0];
  total_nodes_on_surrogate_multigrid += nodes_on_level_of_surrogate_multigrid[0];

  /* FIXME: REALLOCATE AS WE COARSEN */
  multigrid_refine_data_t* coarse_grid_refinement = P4EST_ALLOC(multigrid_refine_data_t, (endlevel-startlevel)* p4est->local_num_quadrants);
  mg_data->coarse_grid_refinement = coarse_grid_refinement;
  
  /* initialize */
  mg_data->fine_nodes = total_nodes_on_multigrid;
  mg_data->coarse_nodes = total_nodes_on_multigrid;

  double* Ae_at0 = P4EST_ALLOC(double, total_nodes_on_multigrid);
  double* err_at0 = P4EST_ALLOC_ZERO(double, total_nodes_on_multigrid);
  double* res_at0 = P4EST_ALLOC(double, total_nodes_on_multigrid);
  double* rres_at0 = P4EST_ALLOC(double, total_nodes_on_multigrid);

  mg_data->Ae = &(Ae_at0)[0];
  mg_data->err = &(err_at0)[0];
  mg_data->res = &(res_at0)[0];
  mg_data->rres = &(rres_at0)[0];
  
  double* u = vecs->u;
  
  problem_data_t vecs_for_cheby_smooth;
  /* these shouldn't change */
  vecs_for_cheby_smooth.scalar_flux_fcn_data = vecs->scalar_flux_fcn_data;
  vecs_for_cheby_smooth.vector_flux_fcn_data = vecs->vector_flux_fcn_data;

  problem_data_t vecs_for_cg_coarse_solve;
  /* these shouldn't change */
  vecs_for_cg_coarse_solve.scalar_flux_fcn_data = vecs->scalar_flux_fcn_data;
  vecs_for_cg_coarse_solve.vector_flux_fcn_data = vecs->vector_flux_fcn_data;
  
  /* these shouldn't change */
  multigrid_cheby_params_t cheby_params;    
  if (mg_data->cg_eigs_iter > 0)
    cg_eigs
      (
       p4est,
       vecs,
       fcns,
       *ghost,
       *ghost_data,
       mg_data->dgmath_jit_dbase,
       mg_data->cg_eigs_iter,
       &max_eigs[endlevel]
      );
  
  max_eigs[endlevel] *= mg_data->max_eig_factor;
  cheby_params.lmin = max_eigs[endlevel]/mg_data->lmax_lmin_rat;
  cheby_params.lmax = max_eigs[endlevel];
  cheby_params.iter = mg_data->smooth_iter;

  if (mg_data->log_option == DEBUG){
    printf("[MULTIGRID_DEBUG]: level %d\n", 0);
    printf("[MULTIGRID_DEBUG]: PREV_PRE_SMOOTH \n");
    printf("[MULTIGRID_DEBUG]: elements on level = %d \n", elements_on_level_of_multigrid[0]);
    printf("[MULTIGRID_DEBUG]: nodes on level = %d \n", nodes_on_level_of_multigrid[0]);
    printf("[MULTIGRID_DEBUG]: max_eigs[%d] = %f\n", endlevel, max_eigs[endlevel]);
    printf("[MULTIGRID_DEBUG]: lmin = %f\n", cheby_params.lmin);
    printf("[MULTIGRID_DEBUG]: lmax = %f\n", cheby_params.lmax);
    printf("[MULTIGRID_DEBUG]: iter = %d\n", cheby_params.iter);    
  }
  
  /* these will change */
  vecs_for_cheby_smooth.Au = vecs->Au;
  vecs_for_cheby_smooth.u = vecs->u; 
  vecs_for_cheby_smooth.rhs = vecs->rhs;
  vecs_for_cheby_smooth.local_nodes = vecs->local_nodes;

  int level_start_fix_this_tho = -2;
  mg_data->mg_state = PREV_PRE_SMOOTH;
  if (mg_data->user_defined_fields && mg_data->mg_update_user_callback != NULL){
    mg_data->mg_update_user_callback(p4est, level_start_fix_this_tho, &vecs_for_cheby_smooth);
  }

  
  multigrid_cheby_smoother
    (
     p4est,
     vecs,
     fcns,
     *ghost,
     *ghost_data,
     mg_data->res,
     &cheby_params
    );


  int level_final_fix_this_tho = -2;
  mg_data->mg_state = PREV_POST_SMOOTH;
  if (mg_data->user_defined_fields && mg_data->mg_update_user_callback != NULL){
    mg_data->mg_update_user_callback(p4est, level_final_fix_this_tho, NULL);
  }
  
  /* initialize error to zero */
  linalg_fill_vec(mg_data->err, 0., mg_data->fine_nodes);
   
  /* we have residual on first level now, store it in temp rr */
  linalg_copy_1st_to_2nd(mg_data->res, mg_data->rres, mg_data->fine_nodes);

  /* initialize stride */
  mg_data->stride = 0;
 
  if (mg_data->save_vtk_snapshot == 1){
    char save_as [500];
    sprintf(save_as, "multigrid_level_%d_preV", endlevel);
    hp_amr_save_to_vtk
      (
       p4est,
       save_as,
       0
      );
  }   
   
  /* Going down the V */
  for (level = endlevel; level > startlevel; --level){

    /* increments the stride */
    p4est_coarsen_ext (p4est,
                       0,
                       1,
                       multigrid_coarsen,
                       multigrid_coarsen_init,
                       NULL
                      );



    /* printf("stride after coarsen = %d\n", mg_data->stride); */
   
    /* does not change the stride */
    element_data_init(p4est, -1);

    mg_data->mg_state = DOWNV_POST_COARSEN;
    if (mg_data->user_defined_fields && mg_data->mg_update_user_callback != NULL){
      mg_data->mg_update_user_callback(p4est, level, NULL);
    }

    /* decrements the stride */
    total_elements_on_surrogate_multigrid += p4est->local_num_quadrants;
    nodes_on_level_of_surrogate_multigrid[endlevel-level+1] = element_data_get_local_nodes(p4est);
    total_nodes_on_surrogate_multigrid += nodes_on_level_of_surrogate_multigrid[endlevel-level+1];
    elements_on_level_of_surrogate_multigrid[endlevel-level+1] = p4est->local_num_quadrants;
    mg_data->stride -= elements_on_level_of_surrogate_multigrid[endlevel-level+1];

  if (mg_data->log_option == DEBUG){
    printf("[MULTIGRID_DEBUG]: level %d\n", endlevel);
    printf("[MULTIGRID_DEBUG]: DOWNV_POST_COARSEN \n");
    printf("[MULTIGRID_DEBUG]: surrogate elements on level = %d \n", elements_on_level_of_surrogate_multigrid[endlevel-level+1]);
    printf("[MULTIGRID_DEBUG]: surrogate nodes on level = %d \n", nodes_on_level_of_surrogate_multigrid[endlevel-level+1]);
    /* printf("[MULTIGRID_DEBUG]: max_eigs[%d] = %f\n", endlevel, max_eigs[endlevel]); */
    /* printf("[MULTIGRID_DEBUG]: lmin = %f\n", cheby_params.lmin); */
    /* printf("[MULTIGRID_DEBUG]: lmax = %f\n", cheby_params.lmax); */
    /* printf("[MULTIGRID_DEBUG]: iter = %f\n", cheby_params.iter);     */
  }

    
    if (mg_data->save_vtk_snapshot == 1){
      char save_as1 [500];
      sprintf(save_as1, "multigrid_level_%d_post_coarsen", level - 1);
      hp_amr_save_to_vtk
        (
         p4est,
         save_as1,
         0
        );
    }
   
    /* does not change the stride */
    p4est_balance_ext (p4est, P4EST_CONNECT_FACE, NULL, multigrid_store_balance_changes);
    element_data_init(p4est, -1);


    mg_data->mg_state = DOWNV_POST_BALANCE;
    if (mg_data->user_defined_fields && mg_data->mg_update_user_callback != NULL){
      mg_data->mg_update_user_callback(p4est, level, NULL);
    }

    
    
    P4EST_FREE(*ghost_data);
    p4est_ghost_destroy (*ghost);
    *ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
    *ghost_data = P4EST_ALLOC (element_data_t,
                               (*ghost)->ghosts.elem_count);
    
    elements_on_level_of_multigrid[endlevel-level+1] = p4est->local_num_quadrants;
    total_elements_on_multigrid += elements_on_level_of_multigrid[endlevel-level+1];
    nodes_on_level_of_multigrid[endlevel-level+1] = element_data_get_local_nodes(p4est);
    total_nodes_on_multigrid += nodes_on_level_of_multigrid[endlevel-level+1];

    mg_data->fine_nodes = nodes_on_level_of_multigrid[endlevel - level];
    mg_data->coarse_nodes = nodes_on_level_of_multigrid[endlevel - level + 1];

    if (mg_data->log_option == DEBUG){
      printf("[MULTIGRID_DEBUG]: level %d\n", endlevel);
      printf("[MULTIGRID_DEBUG]: DOWNV_POST_BALANCE \n");
      printf("[MULTIGRID_DEBUG]: elements on level = %d \n", elements_on_level_of_multigrid[endlevel-level+1]);
      printf("[MULTIGRID_DEBUG]: nodes on level = %d \n", nodes_on_level_of_multigrid[endlevel-level+1]);
      printf("[MULTIGRID_DEBUG]: fine nodes = %d \n", nodes_on_level_of_multigrid[endlevel-level]);
      printf("[MULTIGRID_DEBUG]: coarse nodes = %d \n", nodes_on_level_of_multigrid[endlevel-level] - 1);
      /* printf("[MULTIGRID_DEBUG]: max_eigs[%d] = %f\n", endlevel, max_eigs[endlevel]); */
      /* printf("[MULTIGRID_DEBUG]: lmin = %f\n", cheby_params.lmin); */
      /* printf("[MULTIGRID_DEBUG]: lmax = %f\n", cheby_params.lmax); */
      /* printf("[MULTIGRID_DEBUG]: iter = %f\n", cheby_params.iter);     */
    }
    
    if (mg_data->save_vtk_snapshot == 1){
      char save_as2 [500];
      sprintf(save_as2, "multigrid_level_%d_post_balance", level - 1);
      hp_amr_save_to_vtk
        (
         p4est,
         save_as2,
         0
        );
    }
    Ae_at0 = P4EST_REALLOC(Ae_at0, double, total_nodes_on_multigrid);
    err_at0 = P4EST_REALLOC(err_at0, double, total_nodes_on_multigrid);
    res_at0 = P4EST_REALLOC(res_at0, double, total_nodes_on_multigrid);
    rres_at0 = P4EST_REALLOC(rres_at0, double, total_nodes_on_multigrid);

    mg_data->Ae = &Ae_at0[total_nodes_on_multigrid - mg_data->fine_nodes - mg_data->coarse_nodes];
    mg_data->err = &err_at0[total_nodes_on_multigrid - mg_data->fine_nodes - mg_data->coarse_nodes];
    mg_data->res = &res_at0[total_nodes_on_multigrid - mg_data->fine_nodes - mg_data->coarse_nodes];
    mg_data->rres = &rres_at0[total_nodes_on_multigrid - mg_data->fine_nodes - mg_data->coarse_nodes];
    
    /* set initial guess for mg_data->error */
    linalg_fill_vec(&mg_data->err[mg_data->fine_nodes], 0., mg_data->coarse_nodes);

    /* always zero these before restriction or prolongation */
    mg_data->coarse_stride = 0;
    mg_data->fine_stride = 0;
    mg_data->temp_stride = 0;
    mg_data->intergrid_ptr = &(mg_data->rres)[0];

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
    

    linalg_copy_1st_to_2nd
      (
       &(mg_data->rres)[mg_data->fine_nodes],
       &(mg_data->res)[mg_data->fine_nodes],
       mg_data->coarse_nodes
      );
    
    /* double debug_res_after_restrict = linalg_vec_dot(&(mg_data->res)[mg_data->fine_nodes],  &(mg_data->res)[mg_data->fine_nodes], mg_data->coarse_nodes); */

    /* printf("%d, debug_res_after_restrict = %.20f\n", level, debug_res_after_restrict); */
    
    /* PROBABLY CAN BE REMOVED */
    element_data_init(p4est, -1);

    /* relax solution  */
    if (level > startlevel + 1){

      
      vecs_for_cheby_smooth.Au = &mg_data->Ae[mg_data->fine_nodes];
      vecs_for_cheby_smooth.u = &mg_data->err[mg_data->fine_nodes];
      vecs_for_cheby_smooth.rhs = &(mg_data->res)[mg_data->fine_nodes];
      vecs_for_cheby_smooth.local_nodes = mg_data->coarse_nodes;
      cheby_params.iter = mg_data->smooth_iter;
      
      if (mg_data->cg_eigs_iter > 0)
        cg_eigs
          (
           p4est,
           &vecs_for_cheby_smooth,
           fcns,
           *ghost,
           *ghost_data,
           mg_data->dgmath_jit_dbase,
           mg_data->cg_eigs_iter,
           &max_eigs[level - 1]
          );

      max_eigs[level - 1] *= mg_data->max_eig_factor;
      cheby_params.lmin = max_eigs[level - 1]/mg_data->lmax_lmin_rat;
      cheby_params.lmax = max_eigs[level - 1];
      cheby_params.iter = mg_data->smooth_iter;

      mg_data->mg_state = DOWNV_PRE_SMOOTH;
      if (mg_data->user_defined_fields && mg_data->mg_update_user_callback != NULL){
        mg_data->mg_update_user_callback(p4est, level, &vecs_for_cheby_smooth);
      }
     
      if (mg_data->log_option == DEBUG){
        printf("[MULTIGRID_DEBUG]: level %d\n", endlevel);
        printf("[MULTIGRID_DEBUG]: PREV_PRE_SMOOTH \n");
        printf("[MULTIGRID_DEBUG]: max_eig = %f\n", max_eigs[level - 1]);
        printf("[MULTIGRID_DEBUG]: lmin = %f\n", cheby_params.lmin);
        printf("[MULTIGRID_DEBUG]: lmax = %f\n", cheby_params.lmax);
        printf("[MULTIGRID_DEBUG]: iter = %d\n", cheby_params.iter);    
      }
      
      multigrid_cheby_smoother
        (
         p4est,
         &vecs_for_cheby_smooth,
         fcns,
         *ghost,
         *ghost_data,
         &(mg_data->rres)[mg_data->fine_nodes],
         &cheby_params
        );

      mg_data->mg_state = DOWNV_POST_SMOOTH;
      if (mg_data->user_defined_fields && mg_data->mg_update_user_callback != NULL){
        mg_data->mg_update_user_callback(p4est, level, NULL);
      }
      
    }

    else{

      mg_data->mg_state = COARSE_PRE_SOLVE;
      if (mg_data->user_defined_fields && mg_data->mg_update_user_callback != NULL){
        mg_data->mg_update_user_callback(p4est, level, NULL);
      }
      

      if (mg_data->coarse_solver_type == CG){
        vecs_for_cg_coarse_solve.Au = &mg_data->Ae[mg_data->fine_nodes];
        vecs_for_cg_coarse_solve.u = &mg_data->err[mg_data->fine_nodes];
        vecs_for_cg_coarse_solve.rhs = &(mg_data->res)[mg_data->fine_nodes];
        vecs_for_cg_coarse_solve.local_nodes = mg_data->coarse_nodes;

        multigrid_cg_params_t cg_params;
        cg_params.iter = mg_data->coarse_iter;
        cg_params.rtol = mg_data->coarse_rtol;
        cg_params.mpi_rank = mg_data->mpi_rank;

        multigrid_cg_coarse_solver
          (
           p4est,
           &vecs_for_cg_coarse_solve,
           fcns,
           *ghost,
           *ghost_data,
           &(mg_data->rres[mg_data->fine_nodes]),
           &cg_params
          );
      }
      else {
        vecs_for_cheby_smooth.Au = &mg_data->Ae[mg_data->fine_nodes];
        vecs_for_cheby_smooth.u = &mg_data->err[mg_data->fine_nodes];
        vecs_for_cheby_smooth.rhs = &(mg_data->res)[mg_data->fine_nodes];
        vecs_for_cheby_smooth.local_nodes = mg_data->coarse_nodes;

        if (mg_data->cg_eigs_iter > 0)
          cg_eigs
            (
             p4est,
             &vecs_for_cheby_smooth,
             fcns,
             *ghost,
             *ghost_data,
             mg_data->dgmath_jit_dbase,
             mg_data->cg_eigs_iter,
             &max_eigs[level - 1]
            );

        max_eigs[level - 1] *= mg_data->max_eig_factor;
        cheby_params.lmin = max_eigs[level - 1]/mg_data->lmax_lmin_rat;
        cheby_params.lmax = max_eigs[level - 1];
        cheby_params.iter = mg_data->coarse_iter;
        multigrid_cheby_smoother
          (
           p4est,
           &vecs_for_cheby_smooth,
           fcns,
           *ghost,
           *ghost_data,
           &(mg_data->rres)[mg_data->fine_nodes],
           &cheby_params
          );
      }


      mg_data->mg_state = COARSE_POST_SOLVE;
      if (mg_data->user_defined_fields && mg_data->mg_update_user_callback != NULL){
        mg_data->mg_update_user_callback(p4est, level, NULL);
      }
      
    }

    mg_data->Ae += mg_data->fine_nodes;
    mg_data->err += mg_data->fine_nodes;
    mg_data->res += mg_data->fine_nodes;
    mg_data->rres += mg_data->fine_nodes;
  }

   mg_data->stride = total_elements_on_surrogate_multigrid;
   mg_data->stride -= elements_on_level_of_surrogate_multigrid[0]; /* subtract fine grid */
  
  /* Going up the V */
  for (level = startlevel; level < endlevel; ++level){


    mg_data->mg_state = UPV_PRE_REFINE;
    if (mg_data->user_defined_fields && mg_data->mg_update_user_callback != NULL){
      mg_data->mg_update_user_callback(p4est, level, NULL);
    }
    
    mg_data->fine_nodes = nodes_on_level_of_multigrid[endlevel-(level+1)];
    mg_data->coarse_nodes = nodes_on_level_of_multigrid[endlevel-level];

    mg_data->err -= mg_data->fine_nodes;
    mg_data->Ae -= mg_data->fine_nodes;
    mg_data->rres -= mg_data->fine_nodes;
    mg_data->res -= mg_data->fine_nodes;
    
    mg_data->stride -= elements_on_level_of_surrogate_multigrid[endlevel-level];
    mg_data->fine_stride = 0;
    mg_data->coarse_stride = 0;
    mg_data->temp_stride = 0;
    linalg_copy_1st_to_2nd(&(mg_data->err)[mg_data->fine_nodes], &(mg_data->rres)[mg_data->fine_nodes], mg_data->coarse_nodes);

    mg_data->intergrid_ptr = mg_data->rres;
    
    /* increments the stride */
    p4est_refine_ext(p4est,
                     0,
                     P4EST_QMAXLEVEL,
                     multigrid_refine_and_apply_prolongation,
                     NULL,
                     multigrid_refine_and_apply_prolongation_replace
                    );

    P4EST_FREE(*ghost_data);
    p4est_ghost_destroy (*ghost);
    *ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
    *ghost_data = P4EST_ALLOC (element_data_t,
                               (*ghost)->ghosts.elem_count);
    
    if (mg_data->save_vtk_snapshot == 1){
      char save_as3 [500];
      sprintf(save_as3, "multigrid_level_%d_post_refine", level + 1);
      hp_amr_save_to_vtk
        (
         p4est,
         save_as3,
         0
        );
    }

    /* decrements the stride */
    mg_data->stride -= elements_on_level_of_surrogate_multigrid[endlevel - level];

    /* PROBABLY NOT NEEDED */
    element_data_init(p4est, -1);    
    linalg_vec_axpy(1.0, mg_data->rres, mg_data->err, mg_data->fine_nodes);

    /* If the grid is not the finiest */
    if(level != endlevel-1){

      
      vecs_for_cheby_smooth.Au = mg_data->Ae;
      vecs_for_cheby_smooth.u = mg_data->err;
      vecs_for_cheby_smooth.rhs = mg_data->res;
      vecs_for_cheby_smooth.local_nodes = mg_data->fine_nodes;



      if (mg_data->cg_eigs_iter > 0)
        cg_eigs
          (
           p4est,
           &vecs_for_cheby_smooth,
           fcns,
           *ghost,
           *ghost_data,
           mg_data->dgmath_jit_dbase,
           mg_data->cg_eigs_iter,
           &max_eigs[level + 1]
          );

      max_eigs[level+1] *= mg_data->max_eig_factor;
      cheby_params.lmin = max_eigs[level + 1]/mg_data->lmax_lmin_rat;
      cheby_params.lmax = max_eigs[level + 1];
      cheby_params.iter = mg_data->smooth_iter;

      mg_data->mg_state = UPV_PRE_SMOOTH;
      if (mg_data->user_defined_fields && mg_data->mg_update_user_callback != NULL){
        mg_data->mg_update_user_callback(p4est, level, &vecs_for_cheby_smooth);
      }
      
      multigrid_cheby_smoother
        (
         p4est,
         &vecs_for_cheby_smooth,
         fcns,
         *ghost,
         *ghost_data,
         mg_data->rres,
         &cheby_params
        );
    }

    /* Finest grid smooth */
    else{
      
      mg_data->mg_state = POSTV_PRE_SMOOTH;
      if (mg_data->user_defined_fields && mg_data->mg_update_user_callback != NULL){
        mg_data->mg_update_user_callback(p4est, level, vecs);
      }
      linalg_vec_axpy(1.0, mg_data->err, u, mg_data->fine_nodes);

      if (mg_data->cg_eigs_iter > 0)
        cg_eigs
          (
           p4est,
           vecs,
           fcns,
           *ghost,
           *ghost_data,
           mg_data->dgmath_jit_dbase,
           mg_data->cg_eigs_iter,
           &max_eigs[endlevel]
          );

      max_eigs[endlevel] *= mg_data->max_eig_factor;
      cheby_params.lmin = max_eigs[endlevel]/mg_data->lmax_lmin_rat;
      cheby_params.lmax = max_eigs[endlevel];
      cheby_params.iter = mg_data->smooth_iter;
      multigrid_cheby_smoother
        (
         p4est,
         vecs,
         fcns,
         *ghost,
         *ghost_data,
         mg_data->res,
         &cheby_params
        );
      
      mg_data->vcycle_r2local = linalg_vec_dot(mg_data->res,
                                               mg_data->res,
                                               mg_data->fine_nodes);
      /* p4est_ghost_destroy (ghost); */
      mg_data->mg_state = POSTV_POST_SMOOTH;
      if (mg_data->user_defined_fields && mg_data->mg_update_user_callback != NULL){
        mg_data->mg_update_user_callback(p4est, level, NULL);
      }

      
    }
  }

  /* P4EST_FREE(ghost_data); */
  P4EST_FREE(elements_on_level_of_multigrid);
  P4EST_FREE(nodes_on_level_of_multigrid);
  P4EST_FREE(elements_on_level_of_surrogate_multigrid);
  P4EST_FREE(nodes_on_level_of_surrogate_multigrid);
  P4EST_FREE(coarse_grid_refinement);
  P4EST_FREE(Ae_at0);
  P4EST_FREE(err_at0);
  P4EST_FREE(res_at0);
  P4EST_FREE(rres_at0);
  /* P4EST_FREE(max_eigs); */
  p4est->user_pointer = tmpptr;
  
}

void
multigrid_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 multigrid_data_t* mg_data,
 p4est_ghost_t** ghost,
 element_data_t** ghost_data
)
{

  /*
   * Calculate the initial residual
   * for termination condition
   */
  double* Au; 
  double* rhs;
  int local_nodes = vecs->local_nodes;
  Au = vecs->Au;
  rhs = vecs->rhs;
  double* r = P4EST_ALLOC(double, vecs->local_nodes);
  fcns->apply_lhs(p4est, *ghost, *ghost_data, vecs, mg_data->dgmath_jit_dbase);  
  linalg_vec_axpyeqz(-1., Au, rhs, r, local_nodes);
  double r2_0_local = linalg_vec_dot(r,r,local_nodes);  
  P4EST_FREE(r);
  
  double r2_0_global;
  sc_allreduce
    (
     &r2_0_local,
     &r2_0_global,
     1,
     sc_MPI_DOUBLE,
     sc_MPI_SUM,
     sc_MPI_COMM_WORLD
    );  
  
  double r2_i_local; 
  double r2_i_global = r2_0_global;
  double old_r2_i_global = r2_i_global;
  /*
   * End Termination residual 
   * calculation
   */


  double* err = NULL;
  double* u_analytic = NULL; 
  double err_2_temp, err_2, old_err_2 = -1.;
  int v = 0;
  
  if (mg_data->log_option == ERR_LOG || mg_data->log_option == ERR_AND_EIG_LOG){
    mpi_abort("Needs vecs->u_analytic_fcn to be fixed");
    err = P4EST_ALLOC(double, vecs->local_nodes);
    u_analytic = P4EST_ALLOC(double, vecs->local_nodes);
    /* element_data_init_node_vec( */
    /*                            p4est, */
    /*                            u_analytic, */
    /*                            vecs->u_analytic_fcn, */
    /*                            mg_data->dgmath_jit_dbase */
    /*                           ); */
    linalg_vec_axpyeqz(-1., vecs->u, u_analytic, err, local_nodes);
    err_2_temp = (element_data_compute_l2_norm_sqr_no_local(p4est, err, mg_data->dgmath_jit_dbase));
    sc_reduce
      (
       &err_2_temp,
       &err_2,
       1,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       0,
       sc_MPI_COMM_WORLD
      );
    old_err_2 = err_2;
    if (mg_data->mpi_rank == 0)
      printf("[MG_SOLVER]: %d %.20f %.20f %f %f\n", v, sqrt(err_2), sqrt(r2_i_global), sqrt(err_2/old_err_2), sqrt(r2_i_global/old_r2_i_global));
  }
  else {
    if (mg_data->mpi_rank == 0)
      printf("[MG_SOLVER]: %d %.20f %f\n", v, sqrt(r2_i_global), sqrt(r2_i_global/old_r2_i_global));
  }

  double r2_stop_atol = mg_data->vcycle_rtol*mg_data->vcycle_rtol*r2_0_global + mg_data->vcycle_atol*mg_data->vcycle_atol;

  /**
   * START VCYCLING
   * 
   */
  while
    (
     v <= mg_data->vcycle_iter
     &&
     r2_i_global > r2_stop_atol
    )
    {
      
    if (mg_data->perform_checksum == 1){
      int checksum_b = p4est_checksum(p4est);
      printf("[MG_SOLVER]: Checksum before VCYCLE: %d\n", checksum_b);
    }

    if ( v != 0 && mg_data->max_eig_reuse == 1){
      mg_data->cg_eigs_iter = -1.;
      mg_data->max_eig_factor = 1.;
    }
    multigrid_vcycle(p4est, vecs, fcns, mg_data, ghost, ghost_data);

    if (mg_data->perform_checksum == 1){
      int checksum_a = p4est_checksum(p4est);
      printf("[MG_SOLVER]: Checksum after VCYCLE: %d\n", checksum_a);
    }
    
    /* Calculate for termination check*/
    r2_i_local = mg_data->vcycle_r2local;
    sc_allreduce
      (
       &r2_i_local,
       &r2_i_global,
       1,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       sc_MPI_COMM_WORLD
      );

    double gamma = sqrt(r2_i_global/old_r2_i_global);
    /* print out error/residual diagnostics */
    if (mg_data->log_option == ERR_LOG || mg_data->log_option == ERR_AND_EIG_LOG)
      {
        linalg_vec_axpyeqz(-1., vecs->u, u_analytic, err, local_nodes);
        err_2_temp = (element_data_compute_l2_norm_sqr_no_local(p4est, err,mg_data->dgmath_jit_dbase));
        sc_reduce
          (
           &err_2_temp,
           &err_2,
           1,
           sc_MPI_DOUBLE,
           sc_MPI_SUM,
           0,
           sc_MPI_COMM_WORLD
          );
        if (mg_data->mpi_rank == 0)
          printf("[MG_SOLVER]: %d %.20f %.20f %f %f %f\n",v, sqrt(err_2), sqrt(r2_i_global), sqrt(err_2/old_err_2), sqrt(r2_i_global/old_r2_i_global), pow((r2_i_global/r2_0_global), 1./(v + 1)));
        old_err_2 = err_2;
      }
    else {
      if (mg_data->mpi_rank == 0)
        printf("[MG_SOLVER]: %d %.20f %.25f %.25f\n", v, sqrt(r2_i_global), sqrt(r2_i_global/old_r2_i_global), pow((r2_i_global/r2_0_global), 1./(v + 1)));
    }

    if (mg_data->log_option == RES_AND_EIG_LOG || mg_data->log_option == ERR_AND_EIG_LOG){
      int i;
      if(mg_data->mpi_rank == 0)
        for (i = 0; i < mg_data->num_of_levels; i++){
          printf("[MG_SOLVER]: LEV %d MAX_EIG %.20f\n", i, mg_data->max_eigs[i]);
        }
    }
    
    if (gamma >= .99){
      break;
    } 
    old_r2_i_global = r2_i_global;
    v++;
  }
  
  if (mg_data->log_option == ERR_LOG || mg_data->log_option == ERR_AND_EIG_LOG){
    P4EST_FREE(u_analytic);
    P4EST_FREE(err);
  }

  mg_data->final_vcycles = v;
  /* P4EST_FREE(mg_data->max_eigs); */
}
