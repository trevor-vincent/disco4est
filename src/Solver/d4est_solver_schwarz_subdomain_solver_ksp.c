#define _GNU_SOURCE
#include <pXest.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_solver_schwarz_geometric_data.h>
#include <d4est_solver_schwarz_subdomain_solver_ksp.h>
#include <d4est_solver_schwarz_operators.h>
#include <d4est_solver_schwarz_apply_lhs.h>
#include <d4est_linalg.h>
#include <zlog.h>
#include <ini.h>
#include <time.h>

static
int d4est_solver_schwarz_subdomain_solver_ksp_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_solver_schwarz_subdomain_solver_ksp_data_t* pconfig = user;
  const char* input_section = pconfig->input_section;
  
  if (d4est_util_match_couple(section,input_section,name,"subdomain_atol")) {
    D4EST_ASSERT(pconfig->subdomain_atol == -1);
    pconfig->subdomain_atol = atof(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"subdomain_rtol")) {
    D4EST_ASSERT(pconfig->subdomain_rtol == -1);
    pconfig->subdomain_rtol = atof(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"subdomain_iter")) {
    D4EST_ASSERT(pconfig->subdomain_iter == -1);
    pconfig->subdomain_iter = atoi(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"subdomain_monitor")) {
    D4EST_ASSERT(pconfig->subdomain_monitor == -1);
    pconfig->subdomain_monitor = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void
d4est_solver_schwarz_subdomain_solver_ksp_destroy
(
 void* params
)
{

  d4est_solver_schwarz_subdomain_solver_ksp_data_t* ksp_params
    = params;
  P4EST_FREE(ksp_params); 
}

d4est_solver_schwarz_subdomain_solver_ksp_data_t*
d4est_solver_schwarz_subdomain_solver_ksp_init
(
 p4est_t* p4est,
 const char* input_file,
 const char* input_section
)
{
  d4est_solver_schwarz_subdomain_solver_ksp_data_t* solver_ksp =
    P4EST_ALLOC(d4est_solver_schwarz_subdomain_solver_ksp_data_t, 1);
  
  zlog_category_t* c_default = zlog_get_category("d4est_solver_schwarz_subdomain_solver_ksp");
  solver_ksp->subdomain_iter = -1;
  solver_ksp->subdomain_rtol = -1;
  solver_ksp->subdomain_atol = -1;
  solver_ksp->subdomain_monitor = -1;
  solver_ksp->input_section = input_section;
  
  if(
     ini_parse(input_file,
               d4est_solver_schwarz_subdomain_solver_ksp_input_handler,
               solver_ksp) < 0
  ){
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, solver_ksp->subdomain_iter, -1);
  D4EST_CHECK_INPUT(input_section, solver_ksp->subdomain_rtol, -1);
  D4EST_CHECK_INPUT(input_section, solver_ksp->subdomain_atol, -1);
  D4EST_CHECK_INPUT(input_section, solver_ksp->subdomain_monitor, -1);
  
  if (solver_ksp->subdomain_iter <= 0 ||
      solver_ksp->subdomain_rtol <= 0 ||
      solver_ksp->subdomain_atol <= 0 
     ){
    D4EST_ABORT("Some subdomain solver options are <= 0");
  }

  return solver_ksp;
}

static
PetscErrorCode d4est_solver_schwarz_subdomain_solver_ksp_apply_aij
(
 Mat A,
 Vec x,
 Vec y
)
{
  void           *ctx;
  const double* px;
  double* py;

  /* PetscFunctionBegin; */
  PetscErrorCode ierr = MatShellGetContext( A, &ctx ); CHKERRQ(ierr);
  d4est_solver_schwarz_subdomain_solver_ksp_ctx_t* petsc_ctx = ctx;
  ierr = VecGetArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecGetArray( y, &py ); CHKERRQ(ierr);

  petsc_ctx->apply_lhs->apply_lhs_fcn
    (
     petsc_ctx->p4est,
     petsc_ctx->schwarz_ops->d4est_ops,
     petsc_ctx->d4est_geom,
     petsc_ctx->d4est_quad,
     petsc_ctx->d4est_factors,
     petsc_ctx->ghost,
     petsc_ctx->schwarz_ops,
     petsc_ctx->schwarz_metadata,
     petsc_ctx->schwarz_geometric_data,
     petsc_ctx->subdomain,
     (double*)px,
     py,
     petsc_ctx->apply_lhs->apply_lhs_ctx
    );


  ierr = VecRestoreArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecRestoreArray( y, &py ); CHKERRQ(ierr);
  return ierr;
}



static
PetscErrorCode d4est_solver_schwarz_subdomain_solver_ksp_monitor
(
 KSP ksp,
 PetscInt it,
 PetscReal norm,
 void *ctx
)
{
  d4est_solver_schwarz_subdomain_solver_ksp_ctx_t* petsc_ctx = ctx;
  zlog_category_t* c_default = zlog_get_category("d4est_solver_schwarz_subdomain_solver_ksp");

  if (petsc_ctx->p4est->mpirank == 0){
    double emax, emin;    
    KSPComputeExtremeSingularValues(ksp, &emax, &emin);
    
    zlog_info(c_default, "elems %d nodes %d rank %d sub %d sub_tree %d iter %d r %.15f emax %.15f emin %.15f",
              petsc_ctx->p4est->global_num_quadrants,
              petsc_ctx->d4est_factors->global_nodes,
              petsc_ctx->p4est->mpirank,
              petsc_ctx->subdomain,
              petsc_ctx->schwarz_metadata->subdomain_metadata[petsc_ctx->subdomain].core_tree,
              it,
              norm,
              emax,
              emin);
  }  

  return 0;
}


d4est_solver_schwarz_subdomain_solver_info_t
d4est_solver_schwarz_subdomain_solver_ksp
(
 p4est_t* p4est,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* ghost,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_solver_schwarz_metadata_t* schwarz_metadata,
 d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data,
 d4est_solver_schwarz_apply_lhs_t* apply_lhs,
 double* du_restricted_field_over_subdomain,
 double* rhs_restricted_field_over_subdomain,
 int subdomain,
 void* params
)
{
  d4est_solver_schwarz_subdomain_solver_ksp_data_t* data = params;
  int iter = data->subdomain_iter;
  double atol = data->subdomain_atol;
  double rtol = data->subdomain_rtol;
  int monitor = data->subdomain_monitor;
  
  int nodes
    = schwarz_metadata->subdomain_metadata[subdomain].restricted_nodal_size;

  KSP ksp;
  Vec x,b;
  
  d4est_solver_schwarz_subdomain_solver_ksp_ctx_t petsc_ctx;
  petsc_ctx.p4est = p4est;
  petsc_ctx.d4est_geom = d4est_geom;
  petsc_ctx.d4est_quad = d4est_quad;
  petsc_ctx.d4est_factors = d4est_factors;
  petsc_ctx.ghost = ghost;
  petsc_ctx.schwarz_ops = schwarz_ops;
  petsc_ctx.schwarz_metadata = schwarz_metadata;
  petsc_ctx.schwarz_geometric_data = schwarz_geometric_data;
  petsc_ctx.apply_lhs = apply_lhs;
  petsc_ctx.subdomain = subdomain;  
  petsc_ctx.ksp = &ksp;
  
  double* u = du_restricted_field_over_subdomain;
  double* rhs = rhs_restricted_field_over_subdomain;

  KSPCreate(PETSC_COMM_SELF,&ksp);
  VecCreate(PETSC_COMM_SELF,&x); //CHKERRQ(ierr);
  VecSetSizes(x, nodes, PETSC_DECIDE); //CHKERRQ(ierr);
  VecSetFromOptions(x); //CHKERRQ(ierr);
  VecDuplicate(x,&b); //CHKERRQ(ierr);  

  PetscOptionsClearValue(NULL,"-ksp_monitor");
  PetscOptionsClearValue(NULL,"-ksp_view");
  PetscOptionsClearValue(NULL,"-ksp_converged_reason");
  PetscOptionsClearValue(NULL,"-ksp_monitor_singular_value");
  PetscOptionsClearValue(NULL,"-ksp_type");
  PetscOptionsClearValue(NULL,"-ksp_atol");
  PetscOptionsClearValue(NULL,"-ksp_rtol");
  PetscOptionsClearValue(NULL,"-ksp_max_it");
  PetscOptionsClearValue(NULL,"-ksp_chebyshev_esteig_steps");
  PetscOptionsClearValue(NULL,"-ksp_chebyshev_esteig");
  PetscOptionsClearValue(NULL,"-ksp_chebyshev_esteig_random");
  KSPSetFromOptions(ksp);
  KSPSetTolerances(ksp,rtol,atol,PETSC_DEFAULT,iter);
  KSPSetType(ksp,"fgmres");
  /* Create matrix-free shell for Aij */
  PC pc;
  KSPGetPC(ksp,&pc);
  PCSetType(pc,PCNONE);//CHKERRQ(ierr);


  
  Mat A;
  MatCreateShell
    (
     PETSC_COMM_SELF,
     nodes,
     nodes,
     PETSC_DETERMINE,
     PETSC_DETERMINE,
     (void*)&petsc_ctx,
     &A
    );
  MatShellSetOperation(A,
                       MATOP_MULT,
                       (void(*)())d4est_solver_schwarz_subdomain_solver_ksp_apply_aij
                      );

  KSPSetOperators(ksp,A,A);
  
  double* r = NULL;
  VecPlaceArray(b, rhs);
  VecPlaceArray(x, u);

  if(monitor){
    KSPMonitorSet(ksp, d4est_solver_schwarz_subdomain_solver_ksp_monitor, &petsc_ctx, NULL);
    KSPSetComputeSingularValues(ksp, PETSC_TRUE);
  }
  
  KSPSolve(ksp,b,x);

  double final_res;
  int final_iter;
  KSPGetIterationNumber(ksp, &(final_iter));
  KSPGetResidualNorm(ksp, &(final_res));

  MatDestroy(&A);
  VecResetArray(b);
  VecResetArray(x);
  VecDestroy(&x);
  VecDestroy(&b);
  KSPDestroy(&ksp); 
  
  return (d4est_solver_schwarz_subdomain_solver_info_t){.final_iter = final_iter,
    .final_res = final_res};
}


