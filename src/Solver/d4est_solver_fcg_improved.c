/**
 * @file   d4est_solver_fcg_improved.c
 * @author  <tvincent@cita.utoronto.ca>
 * 
 * @brief Based off of "A massively parallel solver for discrete Poisson-like problems" by Yvan Notay and Artem Napov
 * 
 * 
 */

#define _GNU_SOURCE
#include <d4est_solver_fcg_improved.h>
#include <d4est_util.h>
#include <ini.h>
#include <d4est_linalg.h>
#include <d4est_krylov_pc.h>
#include <d4est_checkpoint.h>
#include <sc_reduce.h>

static
int d4est_solver_fcg_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_solver_fcg_params_t* pconfig = (d4est_solver_fcg_params_t*)user;
  if (d4est_util_match_couple(section,pconfig->input_section,name,"atol")) {
    D4EST_ASSERT(pconfig->atol == -1);
    pconfig->atol = atof(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"rtol")) {
    D4EST_ASSERT(pconfig->rtol == -1);
    pconfig->rtol = atof(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"checkpoint_every_n_krylov_its")) {
    D4EST_ASSERT(pconfig->checkpoint_every_n_krylov_its == 0);
    /* D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1); */
    pconfig->checkpoint_every_n_krylov_its = atoi(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"imax")) {
    D4EST_ASSERT(pconfig->imax == -1);
    D4EST_ASSERT(atoi(value) >= 0);
    pconfig->imax = atoi(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"do_not_use_preconditioner")) {
    D4EST_ASSERT(pconfig->do_not_use_preconditioner == 0);
    D4EST_ASSERT(atoi(value) == 1 || atoi(value) == 0);
    pconfig->do_not_use_preconditioner = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void
d4est_solver_fcg_input
(
 p4est_t* p4est,
 const char* input_file,
 const char* input_section,
 d4est_solver_fcg_params_t* input
)
{

  input->imax = -1;
  input->rtol = -1;
  input->atol = -1;
  input->do_not_use_preconditioner = 0;
  input->checkpoint_every_n_krylov_its = 0;
  
  D4EST_ASSERT(sizeof(input->input_section) <= 50);
  snprintf (input->input_section, sizeof(input->input_section), "%s", input_section);
  
  if (ini_parse(input_file, d4est_solver_fcg_input_handler, input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, input->imax, -1);
  D4EST_CHECK_INPUT(input_section, input->atol, -1);
  D4EST_CHECK_INPUT(input_section, input->rtol, -1);
    
  if(p4est->mpirank == 0){
    zlog_category_t *c_default = zlog_get_category("d4est_solver_fcg");
    zlog_debug(c_default,"imax = %d\n",input->imax);
    zlog_debug(c_default,"atol = %.15f\n",input->atol);
    zlog_debug(c_default,"rtol = %.15f\n",input->rtol);
    zlog_debug(c_default,"do not use preconditioner = %d\n", input->do_not_use_preconditioner);
    zlog_debug(c_default,"checkpoint_every_n_krylov_its = %d\n", input->checkpoint_every_n_krylov_its);
  }
  
}

void
d4est_solver_fcg_solve
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 d4est_ghost_t** ghost,
 d4est_ghost_data_t** ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* fcg_params,
 d4est_krylov_pc_t* d4est_krylov_pc,
 int amr_level,
 int newton_iteration /* set to -1 if not using newton */
)
{
  d4est_solver_fcg_params_t* params = fcg_params;
  zlog_category_t *c_default = zlog_get_category("d4est_solver_fcg");
  int local_nodes = vecs->local_nodes;
  double* Au = vecs->Au;
  double* u = vecs->u;
  double* rhs = vecs->rhs;

  d4est_elliptic_data_t vecs_aux;
  d4est_elliptic_data_copy_ptrs(vecs,&vecs_aux);

  const int imax = params->imax;
  const double atol = params->atol;
  const double rtol = params->rtol;
  /* const int m_max = params->vi; */
  
  double* r_k = P4EST_ALLOC(double, local_nodes);
  double* w_k = P4EST_ALLOC(double, local_nodes);
  double* d_k = P4EST_ALLOC(double, local_nodes);
  double* v_k = P4EST_ALLOC(double, local_nodes);
  double* q_k = P4EST_ALLOC(double, local_nodes);


  /* BEGIN: INITIALIZE PRECOND CTX */
  /***********************/
  krylov_ctx_t krylov_ctx;
  krylov_ctx.p4est = p4est;
  krylov_ctx.vecs = vecs;
  krylov_ctx.fcns = fcns;
  krylov_ctx.ghost = ghost;
  krylov_ctx.ghost_data = ghost_data;
  krylov_ctx.d4est_ops = d4est_ops;
  krylov_ctx.d4est_geom = d4est_geom;
  krylov_ctx.d4est_quad = d4est_quad;
  krylov_ctx.d4est_factors = d4est_factors;
  krylov_ctx.checkpoint_every_n_krylov_its = params->checkpoint_every_n_krylov_its;
  krylov_ctx.last_krylov_checkpoint_it = 0;
  krylov_ctx.amr_level = amr_level;
  krylov_ctx.time_start = clock();
  krylov_ctx.using_newton = (newton_iteration != -1) ? 1 : 0;
  krylov_ctx.newton_iteration = newton_iteration;
  /* krylov_ctx.ksp = &ksp; */
  /* END: INITIALIZE PRECOND CTX */
  /***********************/


  /* BEGIN: COMPUTE TOLERANCE */
  /****************************/
  
  /* Compute Au_0 */
  d4est_elliptic_eqns_apply_lhs
    (
     p4est,
     *ghost,
     *ghost_data,
     fcns,
     vecs,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors
    );

  /* r0 = RHS - Au_0*/
  d4est_linalg_vec_axpyeqz(-1., Au, rhs, r_k, local_nodes);
  
  double r0_norm_global;
  double r0_norm_local = d4est_linalg_vec_dot(r_k,r_k,local_nodes);

  sc_allreduce
    (
     &r0_norm_local,
     &r0_norm_global,
     1,
     sc_MPI_DOUBLE,
     sc_MPI_SUM,
     sc_MPI_COMM_WORLD
    );

  double tol = params->atol + params->rtol*sqrt(r0_norm_global);
  /***********************/
  /* END: COMPUTE TOLERANCE */
  
  double local_dots [4];
  double global_dots [4];

  double alpha_k;
  double beta_k;
  double rho_k;
  double gamma_k;
  
  for (int k = 0; k < imax; k++) {
    /* v_i = B(r_i) */
    if(d4est_krylov_pc != NULL && params->do_not_use_preconditioner == 0){
      d4est_krylov_pc->pc_ctx = &krylov_ctx;
      if (d4est_krylov_pc->pc_setup != NULL){
        d4est_krylov_pc->pc_setup(d4est_krylov_pc);
      }
      d4est_krylov_pc->pc_apply(d4est_krylov_pc, r_k, v_k);
    }
    else {
      /* Identity preconditioner */
      d4est_util_copy_1st_to_2nd(r_k, v_k, local_nodes);
    }

    /* w_i = Av */
    vecs_aux.Au = w_k;
    vecs_aux.u = v_k;
    
    d4est_elliptic_eqns_apply_lhs
      (
       p4est,
       *ghost,
       *ghost_data,
       fcns,
       &vecs_aux,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors
      );

    /* alpha_k = v_k . r_k */
    local_dots[0] = d4est_linalg_vec_dot(v_k, r_k, local_nodes);
    /* beta_k = v_k . w_k */
    local_dots[1] = d4est_linalg_vec_dot(v_k, w_k, local_nodes);

    if (k > 0){
      local_dots[2] = d4est_linalg_vec_dot(v_k, q_k, local_nodes);
      local_dots[3] = d4est_linalg_vec_dot(r_k, r_k, local_nodes);
    }
    
    sc_allreduce
      (
       &local_dots[0],
       &global_dots[0],
       (k > 0) ? 4 : 2,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       sc_MPI_COMM_WORLD
      );

    alpha_k = global_dots[0];
    beta_k = global_dots[1];
    
    if (k > 0){
      /* gamma_k = v_k . q_{k-1} */
      gamma_k = global_dots[2];
      /* d_k = v_k - (gamma_k/rho_{k-1})d_{k-1} */
      d4est_linalg_vec_axpyeqz(-gamma_k/rho_k, d_k, v_k, d_k, local_nodes);
      /* q_k = w_k - (gamma_k/rho_{k-1})q_{k-1} */
      d4est_linalg_vec_axpyeqz(-gamma_k/rho_k, q_k, w_k, q_k, local_nodes);
      /* rho_k = beta_k - gamma_k^2/rho_{k-1} */
      rho_k = beta_k - (gamma_k*gamma_k)/rho_k;      
    }
    else {
      rho_k = beta_k;

      /* d_k = v_k */
      d4est_util_copy_1st_to_2nd(v_k,d_k, local_nodes);
      
      /* q_k = w_k */
      d4est_util_copy_1st_to_2nd(w_k, q_k, local_nodes);
    }

    /* u_{k+1} = u_k + (alpha_k/rho_k)d_k */
    d4est_linalg_vec_axpyeqz(alpha_k/rho_k, d_k, u, u, local_nodes);
    /* r_{k+1} = r_k - (alpha_k/rho_k)q_k */
    d4est_linalg_vec_axpyeqz(-alpha_k/rho_k, q_k, r_k, r_k, local_nodes);

    if (k > 0 && sqrt(global_dots[3]) <= tol){
      break;
    }
    else
      {
        zlog_category_t* c_default = zlog_get_category("d4est_solver_fcg");
        zlog_category_t* its_output = zlog_get_category("d4est_solver_fcg_iteration_info");

        if (p4est->mpirank == 0){
          double duration_seconds = ((double)(clock() - krylov_ctx.time_start)) / CLOCKS_PER_SEC;
          if (newton_iteration < 0){
          zlog_info(its_output, "ELEM NODES AMR_IT KRY_IT KRY_NORM TIME: %d %d %d %d %.25f %f",
                    krylov_ctx.p4est->global_num_quadrants,
                    krylov_ctx.d4est_factors->global_nodes,
                    krylov_ctx.amr_level,
                    k,
                    (k > 0) ? sqrt(global_dots[3]) : sqrt(r0_norm_global),
                    duration_seconds);
          }
          else {
          zlog_info(its_output, "ELEM NODES AMR_IT NEWT_IT KRY_IT KSP_NORM TIME: %d %d %d %d %d %.25f %f",
                    krylov_ctx.p4est->global_num_quadrants,
                    krylov_ctx.d4est_factors->global_nodes,
                    krylov_ctx.amr_level,
                    newton_iteration,
                    k,
                    (k > 0) ? sqrt(global_dots[3]) : sqrt(r0_norm_global),
                    duration_seconds);
          }
        }  
        if (krylov_ctx.checkpoint_every_n_krylov_its > 0 && krylov_ctx.amr_level >= 0){
          if ( (k - krylov_ctx.last_krylov_checkpoint_it) >= krylov_ctx.checkpoint_every_n_krylov_its ){
            char* output = NULL;
            asprintf(&output,"checkpoint_krylov_%d", k);

            /* d4est_elliptic_data_t* vecs = /\* petsc_ctx-> *\/vecs; */
            d4est_checkpoint_save
              (
               amr_level,
               output,
               p4est,
               NULL,
               NULL,
               (const char * []){"u", NULL},
               (hid_t []){H5T_NATIVE_DOUBLE},
               (int []){vecs->local_nodes},
               (void* []){vecs->u}
              );

            krylov_ctx.last_krylov_checkpoint_it = k;
            free(output);
          }
        }
      }
    
  }


  P4EST_FREE(r_k);
  P4EST_FREE(w_k);
  P4EST_FREE(d_k);
  P4EST_FREE(v_k);
  P4EST_FREE(q_k);
}
