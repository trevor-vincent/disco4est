#include "../pXest/pXest.h"
#include "../LinearAlgebra/d4est_linalg.h"
#include "../Utilities/util.h"
#include "../Solver/matrix_sym_tester.h"
#include <ini.h>
#include <krylov_pc.h>
#include <krylov_petsc.h>
#include "sc_reduce.h"

typedef struct {

  double rtol;
  double atol;
  int imax;
  int imin;
  int monitor;
  
} newton_d4est_params_t;


static
int newton_d4est_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  newton_d4est_params_t* pconfig = (newton_d4est_params_t*)user;
  
  if (util_match_couple(section,"newton_d4est",name,"atol")) {
    D4EST_ASSERT(pconfig->atol == -1);
    pconfig->atol = atof(value);
  }
  else if (util_match_couple(section,"newton_d4est",name,"rtol")) {
    D4EST_ASSERT(pconfig->rtol == -1);
    pconfig->rtol = atof(value);
  }
  else if (util_match_couple(section,"newton_d4est",name,"imax")) {
    D4EST_ASSERT(pconfig->imax == -1);
    pconfig->imax = atoi(value);
  }
  else if (util_match_couple(section,"newton_d4est",name,"monitor")) {
    D4EST_ASSERT(pconfig->monitor == -1);
    pconfig->monitor = atoi(value);
  }
  else if (util_match_couple(section,"newton_d4est",name,"imin")) {
    D4EST_ASSERT(pconfig->imin == -1);
    pconfig->imin = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;  
}


newton_d4est_params_t
newton_d4est_input
(
 p4est_t* p4est,
 const char* input_file
)
{
  newton_d4est_params_t input;
  input.atol = -1;
  input.rtol = -1;
  input.imin = -1;
  input.imax = -1;
  input.monitor = -1;
  
  if (ini_parse(input_file, newton_d4est_input_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("newton", input.atol, -1);
  D4EST_CHECK_INPUT("newton", input.rtol, -1);
  D4EST_CHECK_INPUT("newton", input.imin, -1);
  D4EST_CHECK_INPUT("newton", input.imax, -1);
  D4EST_CHECK_INPUT("newton", input.monitor, -1);


   
  if(p4est->mpirank == 0){
    printf("[NEWTON_D4EST]: atol = %f\n", input.atol);
    printf("[NEWTON_D4EST]: rtol = %f\n", input.rtol);
    printf("[NEWTON_D4EST]: imin = %d\n", input.imin);
    printf("[NEWTON_D4EST]: imax = %d", input.imax);
    printf("[NEWTON_D4EST]: monitor = %d\n", input.monitor);
  }
  
  return input;
}

typedef
void
(*krylov_solver_fcn_t)
(
 p4est_t*,
 d4est_elliptic_data_t*,
 d4est_elliptic_eqns_t*,
 double*,
 int
);


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
newton_d4est_solve
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 p4est_ghost_t** ghost,
 void** ghost_data, 
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 const char* input_file,
 krylov_pc_t* krylov_pc
)
{

  krylov_petsc_params_t petsc_params;
  krylov_petsc_input(p4est,
                     input_file,
                     "krylov_petsc",
                     "[KRYLOV_PETSC_FOR_NEWTON_D4EST]",
                     &petsc_params);      
  
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
  
  /******** external parameters ********/
  newton_d4est_params_t nr_params = newton_d4est_input(p4est,input_file);
  
  double atol = nr_params.atol;
  double rtol = nr_params.rtol;
  int maxit = nr_params.imax;
  int minit = nr_params.imin;
  
  vecs_for_res_build.rhs = vecs->rhs;
  vecs_for_res_build.Au = f0;
  vecs_for_res_build.u = x;
  
  /* build initial residual vector */
  fcns->build_residual(p4est, *ghost, *ghost_data, &vecs_for_res_build, d4est_ops, d4est_geom, d4est_quad);
  
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
  /* printf("[NEWTON_SOLVER]: INITIAL FNRM = %f\n", fnrm); */
  
  double fnrmo = 1.;
  double stop_tol = atol + rtol*fnrm;
  int itc = 0;


  if (p4est->mpirank == 0 && nr_params.monitor){
    printf("[NEWTON_D4EST]: ITER %03d PRE-FNRM %.30f POST-FNRM  %.30f\n", itc, fnrmo,  fnrm);
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
    d4est_linalg_fill_vec(vecs_for_linsolve.u, 0., n);

    /* petsc_ctx_t petsc_ctx; */
    /* if(krylov_pc != NULL){ */
    /*   petsc_ctx.p4est = p4est; */
    /*   petsc_ctx.vecs = vecs; */
    /*   petsc_ctx.fcns = fcns; */
    /*   petsc_ctx.ghost = ghost; */
    /*   petsc_ctx.ghost_data = ghost_data; */
    /*   petsc_ctx.d4est_ops = d4est_ops; */
    /*   petsc_ctx.d4est_geom = d4est_geom; */
    /*   krylov_pc->pc_ctx = &petsc_ctx; */
    /*   if(krylov_pc->pc_setup != NULL){ */
    /*     krylov_pc->pc_setup(krylov_pc); */
    /*   } */
    /* } */
    
    krylov_petsc_solve
      (
       p4est,
       &vecs_for_linsolve,
       fcns,
       ghost,
       ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &petsc_params,
       krylov_pc
      );
    
    /* xt = x + lambda*step */
    d4est_linalg_vec_axpyeqz(1.0, step, x, xt, n);

    /* calculate new residual vector */
    vecs_for_res_build.u = xt;
    vecs_for_res_build.Au = ft;
    fcns->build_residual(p4est, *ghost, *ghost_data, &vecs_for_res_build, d4est_ops, d4est_geom, d4est_quad);    

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

    if (p4est->mpirank == 0 && nr_params.monitor){
      printf("[NEWTON_D4EST]: ITER %03d PRE-FNRM %.30f POST-FNRM  %.30f\n" ,itc, fnrmo,  fnrm);
    }
    
  }
  
 /* clean: */
  /* nr_params.final_iter = itc; */
  /* nr_params.final_fnrm = fnrm; */

  if (fnrm > stop_tol && ierr == 0){
    ierr = 1;
  }

  if(vecs->u != x){
    /* vecs->u = x; */
    d4est_linalg_copy_1st_to_2nd(x, vecs->u, vecs->local_nodes);
    P4EST_FREE(x);
  }
  else {
    P4EST_FREE(xt);
  }

  P4EST_FREE(Jstep);
  P4EST_FREE(step);
  P4EST_FREE(f0);
  P4EST_FREE(ft);

  return ierr;
}

