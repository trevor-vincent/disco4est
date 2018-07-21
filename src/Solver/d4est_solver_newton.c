#include <pXest.h>
#include <d4est_linalg.h>
#include <d4est_solver_newton.h>
#include <d4est_util.h>
#include <ini.h>
#include <d4est_krylov_pc.h>
#include <krylov_petsc.h>
#include <sc_reduce.h>
#include <zlog.h>
#include <d4est_solver_krylov.h>

static
int d4est_solver_newton_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_solver_newton_params_t* pconfig = (d4est_solver_newton_params_t*)user;
  
  if (d4est_util_match_couple(section,"d4est_solver_newton",name,"atol")) {
    D4EST_ASSERT(pconfig->atol == -1);
    pconfig->atol = atof(value);
  }
  else if (d4est_util_match_couple(section,"d4est_solver_newton",name,"rtol")) {
    D4EST_ASSERT(pconfig->rtol == -1);
    pconfig->rtol = atof(value);
  }
  else if (d4est_util_match_couple(section,"d4est_solver_newton",name,"imax")) {
    D4EST_ASSERT(pconfig->imax == -1);
    pconfig->imax = atoi(value);
  }
  else if (d4est_util_match_couple(section,"d4est_solver_newton",name,"monitor")) {
    D4EST_ASSERT(pconfig->monitor == -1);
    pconfig->monitor = atoi(value);
  }
  else if (d4est_util_match_couple(section,"d4est_solver_newton",name,"imin")) {
    D4EST_ASSERT(pconfig->imin == -1);
    pconfig->imin = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


d4est_solver_newton_params_t
d4est_solver_newton_input
(
 p4est_t* p4est,
 const char* input_file
)
{
  d4est_solver_newton_params_t input;
  input.atol = -1;
  input.rtol = -1;
  input.imin = -1;
  input.imax = -1;
  input.monitor = -1;
  
  if (ini_parse(input_file, d4est_solver_newton_input_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("d4est_solver_newton", input.atol, -1);
  D4EST_CHECK_INPUT("d4est_solver_newton", input.rtol, -1);
  D4EST_CHECK_INPUT("d4est_solver_newton", input.imin, -1);
  D4EST_CHECK_INPUT("d4est_solver_newton", input.imax, -1);
  D4EST_CHECK_INPUT("d4est_solver_newton", input.monitor, -1);


   
  if(p4est->mpirank == 0){
    zlog_category_t *c_default = zlog_get_category("solver_newton");
    zlog_debug(c_default, "atol = %f", input.atol);
    zlog_debug(c_default, "rtol = %f", input.rtol);
    zlog_debug(c_default, "imin = %d", input.imin);
    zlog_debug(c_default, "imax = %d", input.imax);
    zlog_debug(c_default, "monitor = %d", input.monitor);
  }
  
  return input;
}


/**
 * Make sure Jacobian has zeroed boundary conditions and build residual has normal
 * boundary condtions
 *
 * @param p4est
 * @param vecs
 * @param fcns
 * @param params
 * @param ghost
 * @param ghost_data
 * @param d4est_ops
 * @param krylov_solve
 *
 * @return
 */
int
d4est_solver_newton_solve
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
 d4est_solver_newton_params_t* nr_params,
 d4est_solver_krylov_fcn_t krylov_fcn,
 void* krylov_fcn_params,
 d4est_krylov_pc_t* d4est_krylov_pc
)
{
  zlog_category_t *c_default = zlog_get_category("solver_newton");
  if(p4est->mpirank == 0)
    zlog_info(c_default, "Performing Newton solve...");

  int ierr = 0;
  int local_nodes = vecs->local_nodes;
  int n = local_nodes;
  d4est_elliptic_data_t vecs_for_linsolve;
  d4est_elliptic_data_t vecs_for_res_build;

  double* xt = P4EST_ALLOC(double, local_nodes);
  double* ft = P4EST_ALLOC(double, local_nodes);
  double* f0 = P4EST_ALLOC(double, local_nodes);
  double* step = P4EST_ALLOC(double, local_nodes);
  double* Jstep = P4EST_ALLOC(double, local_nodes);
  
  /* these don't change */
  double* x = vecs->u;
  d4est_elliptic_data_copy_ptrs(vecs, &vecs_for_linsolve);
  d4est_elliptic_data_copy_ptrs(vecs, &vecs_for_res_build);
    
  double atol = nr_params->atol;
  double rtol = nr_params->rtol;
  int maxit = nr_params->imax;
  int minit = nr_params->imin;
  
  vecs_for_res_build.rhs = vecs->rhs;
  vecs_for_res_build.Au = f0;
  vecs_for_res_build.u = x;
  
  /* build initial residual vector */
  d4est_elliptic_eqns_build_residual
    (
     p4est,
     *ghost,
     *ghost_data,
     fcns,
     &vecs_for_res_build,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors
    );
  
  double fnrm = d4est_linalg_vec_dot(f0,f0,n);
  double fnrm_global;
  
  sc_allreduce
    (
     &fnrm,
     &fnrm_global,
     1,
     sc_MPI_DOUBLE,
     sc_MPI_SUM,
     sc_MPI_COMM_WORLD
    );

  fnrm = sqrt(fnrm_global);
  
  double fnrmo = 1.;
  double stop_tol = atol + rtol*fnrm;
  int itc = 0;


  if (p4est->mpirank == 0 && nr_params->monitor){
    zlog_debug(c_default, "ITER %03d INITIAL FNRM  %.30f", itc, fnrm);
  }
  
  while((fnrm > stop_tol || itc < minit) && (itc < maxit)){

    /* double ratio = fnrm/fnrmo; */
    fnrmo = fnrm;
    itc++;
    d4est_linalg_vec_scale(-1., f0, n);
    
    vecs_for_linsolve.u0 = x;
    vecs_for_linsolve.rhs = f0;
    vecs_for_linsolve.u = step;
    vecs_for_linsolve.Au = Jstep;

    /* set initial guess */
    d4est_util_fill_array(vecs_for_linsolve.u, 0., n);

    krylov_fcn
      (
       p4est,
       &vecs_for_linsolve,
       fcns,
       ghost,
       ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       krylov_fcn_params,
       d4est_krylov_pc
      );


    /* xt = x + lambda*step */
    d4est_linalg_vec_axpyeqz(1.0, step, x, xt, n);

    /* calculate new residual vector */
    vecs_for_res_build.u = xt;
    vecs_for_res_build.Au = ft;
    d4est_elliptic_eqns_build_residual
      (
       p4est,
       *ghost,
       *ghost_data,
       fcns,
       &vecs_for_res_build,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors
      );

    /* flip pointers */
    double* tmp = x;
    x = xt;
    xt = tmp;
    tmp = f0;
    f0 = ft;
    ft = tmp;

    /* calculate new residual norm */
    fnrm = d4est_linalg_vec_dot(f0,f0,n);

    sc_allreduce
      (
       &fnrm,
       &fnrm_global,
       1,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       sc_MPI_COMM_WORLD
      );

    fnrm = sqrt(fnrm_global);

    if (p4est->mpirank == 0 && nr_params->monitor){
      zlog_debug(c_default, "ITER %03d PRE-FNRM %.15e POST-FNRM  %.15e" ,itc, fnrmo,  fnrm);
    }
    
  }
  
  if (fnrm > stop_tol && ierr == 0){
    ierr = 1;
  }

  if(vecs->u != x){
    d4est_util_copy_1st_to_2nd(x, vecs->u, vecs->local_nodes);
    P4EST_FREE(x);
  }
  else {
    P4EST_FREE(xt);
  }

  P4EST_FREE(Jstep);
  P4EST_FREE(step);
  P4EST_FREE(f0);
  P4EST_FREE(ft);

  if (p4est->mpirank == 0)
    zlog_info(c_default, "Newton solve complete.");

  return ierr;
}
