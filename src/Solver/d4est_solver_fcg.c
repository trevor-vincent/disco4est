/**
 * @file   d4est_solver_fcg.c
 * 
 * @brief Based off of arxiv 1511.07226
 * 
 * 
 */


#include <d4est_solver_fcg.h>
#include <d4est_util.h>
#include <ini.h>
#include <d4est_linalg.h>
#include <krylov_pc.h>
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
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"monitor")) {
    D4EST_ASSERT(pconfig->monitor == -1);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
    pconfig->monitor = atoi(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"iter")) {
    D4EST_ASSERT(pconfig->iter == -1);
    D4EST_ASSERT(atoi(value) >= 0);
    pconfig->iter = atoi(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"vi")) {
    D4EST_ASSERT(pconfig->vi == -1);
    D4EST_ASSERT(atoi(value) >= 1);
    pconfig->vi = atoi(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"precond_flag")) {
    D4EST_ASSERT(pconfig->precond_flag == -1);
    D4EST_ASSERT(atoi(value) == 1 || atoi(value) == 0);
    pconfig->precond_flag = atoi(value);
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
 const char* printf_prefix,
 d4est_solver_fcg_params_t* input
)
{
  input->monitor = -1;
  input->iter = -1;
  input->rtol = -1;
  input->atol = -1;
  input->vi = -1;
  input->precond_flag = -1;

  D4EST_ASSERT(sizeof(input->input_section) <= 50);
  snprintf (input->input_section, sizeof(input->input_section), "%s", input_section);
  
  if (ini_parse(input_file, d4est_solver_fcg_input_handler, input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, input->iter, -1);
  D4EST_CHECK_INPUT(input_section, input->monitor, -1);
  D4EST_CHECK_INPUT(input_section, input->atol, -1);
  D4EST_CHECK_INPUT(input_section, input->rtol, -1);
  D4EST_CHECK_INPUT(input_section, input->vi, -1);
  D4EST_CHECK_INPUT(input_section, input->precond_flag, -1);
    
  if(p4est->mpirank == 0){
    zlog_category_t *c_default = zlog_get_category("d4est_solver_fcg");
    zlog_debug(c_default,"%s: iter = %d\n",printf_prefix, input->iter);
    zlog_debug(c_default,"%s: atol = %.15f\n",printf_prefix, input->atol);
    zlog_debug(c_default,"%s: rtol = %.15f\n",printf_prefix, input->rtol);
    zlog_debug(c_default,"%s: vi = %d\n",printf_prefix, input->vi);
    zlog_debug(c_default,"%s: monitor = %d\n",printf_prefix, input->monitor);
    zlog_debug(c_default,"%s: use preconditioner = %d\n",printf_prefix, input->precond_flag);
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
 krylov_pc_t* krylov_pc
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

  const int imax = params->iter;
  const double atol = params->atol;
  const double rtol = params->rtol;
  const int m_max = params->vi;
  
  double* r = P4EST_ALLOC(double, local_nodes);
  double* w = P4EST_ALLOC(double, local_nodes);
  double** d_storage = P4EST_ALLOC(double*, m_max);
  double** Ad_storage = P4EST_ALLOC(double*, m_max);
  for (int m = 0; m < m_max; m++){
    d_storage[m] = P4EST_ALLOC(double, local_nodes);
    Ad_storage[m] = P4EST_ALLOC(double, local_nodes);
  }
  double* d_dot_Ad_storage = P4EST_ALLOC(double, m_max);
  double* local_w_i_dot_Ad_k = P4EST_ALLOC(double, m_max);
  double* global_w_i_dot_Ad_k = P4EST_ALLOC(double, m_max);
  
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

  krylov_ctx_t krylov_ctx;
  if (krylov_pc != NULL && params->precond_flag == 1){
    krylov_ctx.p4est = p4est;
    krylov_ctx.vecs = vecs;
    krylov_ctx.fcns = fcns;
    krylov_ctx.ghost = ghost;
    krylov_ctx.ghost_data = ghost_data;
    krylov_ctx.d4est_ops = d4est_ops;
    krylov_ctx.d4est_geom = d4est_geom;
    krylov_ctx.d4est_quad = d4est_quad;
    krylov_ctx.d4est_factors = d4est_factors;
  }  

  /* Compute residual r = RHS - Au_0*/
  d4est_linalg_vec_axpyeqz(-1.,Au,rhs, r, local_nodes);

  double r0_norm_global;
  double r0_norm_local = d4est_linalg_vec_dot(r,r,local_nodes);

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
  
  for (int i = 0; i < imax; i++) {

    /* w_i = B(r_i) */
    if(krylov_pc != NULL && params->precond_flag == 1){
      krylov_pc->pc_ctx = &krylov_ctx;
      if (krylov_pc->pc_setup != NULL){
        krylov_pc->pc_setup(krylov_pc);
      }
      krylov_pc->pc_apply(krylov_pc, r, w);
    }
    else {
      /* Identity preconditioner */
      d4est_util_copy_1st_to_2nd(r, w, local_nodes);
    }

    if (i != 0){

      int mi = d4est_util_max_int(1, i % (m_max + 1));
      D4EST_ASSERT(mi <= m_max);
      
      for (int k = (i-mi); k < i-1; k++){
        int kmod = k % m_max;
        local_w_i_dot_Ad_k[kmod] = d4est_linalg_vec_dot(w, Ad_storage[kmod], local_nodes);
      }

      sc_allreduce
        (
         local_w_i_dot_Ad_k,
         global_w_i_dot_Ad_k,
         mi,
         sc_MPI_DOUBLE,
         sc_MPI_SUM,
         sc_MPI_COMM_WORLD
        );

      for (int n = 0; n < local_nodes; n++){
        double sum = 0.;
        for (int k = (i-mi); k < i-1; k++){
          int kmod = k % m_max;
          sum += d_storage[kmod][n]*global_w_i_dot_Ad_k[kmod]/d_dot_Ad_storage[kmod];
        }
        w[n] -= sum;
      }
    }
  
    double* Adi = Ad_storage[i % m_max];
    double* di = d_storage[i % m_max];
    d4est_util_copy_1st_to_2nd(w,di,local_nodes);
    
    vecs_aux.Au = Adi;
    vecs_aux.u = di;
    
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
    
    double di_dot_ri_local = d4est_linalg_vec_dot(di,r,local_nodes);
    double di_dot_Adi_local = d4est_linalg_vec_dot(di,Adi,local_nodes);
    double ri_dot_ri_local = d4est_linalg_vec_dot(r,r,local_nodes);
    double local_dots [] = {di_dot_ri_local,di_dot_Adi_local,ri_dot_ri_local};
    double global_dots [3];

    sc_allreduce
      (
       local_dots,
       global_dots,
       3,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       sc_MPI_COMM_WORLD
      );
    
    double alpha = global_dots[0]/global_dots[1];
    if (params->monitor == 1){
      zlog_info(c_default, "iter = %d, r = %.25f\n", i, sqrt(global_dots[2]));
    }
    if (sqrt(global_dots[2]) < tol){
      break;
    }
    
    d4est_linalg_vec_axpyeqz(alpha, di, u, u, local_nodes);
    d4est_linalg_vec_axpyeqz(-alpha, Adi, r, r, local_nodes);
    d_dot_Ad_storage[i % m_max] = global_dots[1];
  }

  P4EST_FREE(r);
  P4EST_FREE(w);
  for (int m = 0; m < m_max; m++){
    P4EST_FREE(d_storage[m]);
    P4EST_FREE(Ad_storage[m]);
  }
  P4EST_FREE(d_storage);
  P4EST_FREE(Ad_storage);
  P4EST_FREE(d_dot_Ad_storage);
  P4EST_FREE(local_w_i_dot_Ad_k);
  P4EST_FREE(global_w_i_dot_Ad_k);
}

