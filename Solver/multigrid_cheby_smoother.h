#ifndef MULTIGRID_CHEBY_SMOOTHER_NEW_H
#define MULTIGRID_CHEBY_SMOOTHER_NEW_H 

static void 
multigrid_cheby_smoother_iterate
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 double* r,
 multigrid_cheby_params_t* cheby_params
)
{
  const int iter = cheby_params->iter;
  const double lmin = cheby_params->lmin;
  const double lmax = cheby_params->lmax;

  multigrid_data_t* mg_data = (multigrid_data_t*) p4est->user_pointer;
  dgmath_jit_dbase_t* dgmath_jit_dbase = mg_data->dgmath_jit_dbase;
 
  int i;
  double d = (lmax + lmin)*.5;
  double c = (lmax - lmin)*.5;

  int local_nodes = vecs->local_nodes;
  double* Au = vecs->Au;
  double* u = vecs->u;
  double* rhs = vecs->rhs;
  
  double* p;
  double alpha,beta;
  p = P4EST_ALLOC(double, local_nodes);
  
  linalg_fill_vec(p, 0., local_nodes);
  for (i = 0; i < iter; i++){
    /* calculate residual r = rhs - Au */
    fcns->apply_lhs(p4est, ghost, ghost_data, vecs, dgmath_jit_dbase);
    linalg_copy_1st_to_2nd(Au, r, local_nodes);
    linalg_vec_xpby(rhs, -1., r, local_nodes);
    
    if (i == 0)
      alpha = 1./d;
    else if (i == 1)
      alpha = 2.*d/(2*d*d - c*c);
    else
      alpha = 1./(d-(alpha*c*c/4.));

    beta = alpha*d - 1.;
   
    linalg_vec_scale(alpha,r,local_nodes);
    linalg_vec_xpby(&r[0], beta, &p[0], local_nodes);   
    linalg_vec_axpy(1., p, u, local_nodes);
  }

  /* calculate the residual */
  fcns->apply_lhs(p4est, ghost, ghost_data, vecs, dgmath_jit_dbase);
  linalg_copy_1st_to_2nd(Au, r, local_nodes);
  linalg_vec_xpby(rhs, -1., r, local_nodes);

  P4EST_FREE(p);
  /* P4EST_FREE(ghost_data); */
  /* p4est_ghost_destroy (ghost); */
}

void
multigrid_cheby_smoother
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 double* r,
 multigrid_data_t* mg_data,
 int level
)
{
  multigrid_cheby_params_t cheby_params;
  
  if (mg_data->solve_for_eigs){
    
      cg_eigs
        (
         p4est,
         vecs,
         fcns,
         ghost,
         ghost_data,
         mg_data->dgmath_jit_dbase,
         mg_data->cg_eigs_iter,
         &mg_data->max_eigs[level],
         mg_data->cg_eigs_use_zero_vec_as_initial
        );

      mg_data->max_eigs[level] *= mg_data->max_eig_factor;
  }
  
  cheby_params.lmin = mg_data->max_eigs[level]/mg_data->lmax_lmin_rat;
  cheby_params.lmax = mg_data->max_eigs[level];
  cheby_params.iter = mg_data->smooth_iter;

  /* if (mg_data->log_option == DEBUG){ */
  /* printf("[MG_CHEBY_SMOOTHER]: Lev %d Max_eig %.25f Iter %d\n", level, mg_data->max_eigs[level], mg_data->smooth_iter); */
  /* }   */
  
  multigrid_cheby_smoother_iterate
    (
     p4est,
     vecs,
     fcns,
     ghost,
     ghost_data,
     r,
     &cheby_params
    );
}

#endif
