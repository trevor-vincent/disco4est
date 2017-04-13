#include <util.h>
#include "../pXest/pXest.h"
#include "../LinearAlgebra/linalg.h"
#include "../Solver/krylov_petsc.h"
#include "petscsnes.h"
#include <krylov_pc.h>
#include <ini.h>

static
int krylov_petsc_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  krylov_petsc_params_t* pconfig = (krylov_petsc_params_t*)user;
  if (util_match_couple(section,pconfig->input_section,name,"ksp_atol")) {
    mpi_assert(pconfig->ksp_atol[0] == '*');
    snprintf (pconfig->ksp_atol, sizeof(pconfig->ksp_atol), "%s", value);
    /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_atol",value); */
  }
  else if (util_match_couple(section,pconfig->input_section,name,"ksp_rtol")) {
    mpi_assert(pconfig->ksp_rtol[0] == '*');
    snprintf (pconfig->ksp_rtol, sizeof(pconfig->ksp_rtol), "%s", value);

    /* pconfig->ksp_rtol = atof(value); */
    /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_rtol",value); */
  }
  else if (util_match_couple(section,pconfig->input_section,name,"ksp_max_it")) {
    mpi_assert(pconfig->ksp_max_it[0] == '*');
    snprintf (pconfig->ksp_max_it, sizeof(pconfig->ksp_max_it), "%s", value);
    mpi_assert(atoi(value) > -1);
    /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_max_it",value); */
  }
  else if (util_match_couple(section,pconfig->input_section,name,"ksp_view")) {
    mpi_assert(pconfig->ksp_view == -1);
    pconfig->ksp_view = atoi(value);
    mpi_assert(atoi(value) == 0 || atoi(value) == 1);
    /* if (atoi(value) == 1){ */
      /* PetscOptionsSetValue(pconfig->petsc_options,"-s_view",""); */
    /* } */
  }  
  else if (util_match_couple(section,pconfig->input_section,name,"ksp_monitor")) {
    mpi_assert(pconfig->ksp_monitor == -1);
    pconfig->ksp_monitor = atoi(value);
    mpi_assert(atoi(value) == 0 || atoi(value) == 1);
    /* if (atoi(value) == 1){ */
      /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_monitor",""); */
    /* } */
  } 
  else if (util_match_couple(section,pconfig->input_section,name,"ksp_converged_reason")) {
    mpi_assert(pconfig->ksp_converged_reason == -1);
    pconfig->ksp_converged_reason = atoi(value);
    mpi_assert(atoi(value) == 0 || atoi(value) == 1);
    /* if (atoi(value) == 1){ */
      /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_converged_reason",""); */
    /* } */
  }    
  else if (util_match_couple(section,pconfig->input_section,name,"ksp_initial_guess_nonzero")) {
    mpi_assert(pconfig->ksp_initial_guess_nonzero == -1);
    pconfig->ksp_initial_guess_nonzero = atoi(value);
    mpi_assert(atoi(value) == 0 || atoi(value) == 1);
    /* if (atoi(value) == 1){ */
      /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_initial_guess_nonzero",""); */
    /* } */
  }  
  else if (util_match_couple(section,pconfig->input_section,name,"ksp_type")) {
    mpi_assert(pconfig->ksp_type[0] == '*');
    snprintf (pconfig->ksp_type, sizeof(pconfig->ksp_type), "%s", value);
    /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_type",value); */
  }
  else if (util_match_couple(section,pconfig->input_section,name,"ksp_chebyshev_esteig_steps")) {
    mpi_assert(pconfig->ksp_chebyshev_esteig_steps[0] == '*');
    snprintf (pconfig->ksp_chebyshev_esteig_steps, sizeof(pconfig->ksp_chebyshev_esteig_steps), "%s", value);
    /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_chebyshev_esteig_steps",value); */
  }
  else if (util_match_couple(section,pconfig->input_section,name,"ksp_chebyshev_esteig")) {
    mpi_assert(pconfig->ksp_chebyshev_esteig[0] == '*');
    snprintf (pconfig->ksp_chebyshev_esteig, sizeof(pconfig->ksp_chebyshev_esteig), "%s", value);
    /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_chebyshev_esteig",value); */
  }
  else if (util_match_couple(section,pconfig->input_section,name,"ksp_chebyshev_esteig_random")) {
    mpi_assert(pconfig->ksp_chebyshev_esteig_random == -1);
    pconfig->ksp_chebyshev_esteig_random = atoi(value);
    mpi_assert(atoi(value) == 0 || atoi(value) == 1);
    /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_chebyshev_esteig_random",value); */
  }  
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void
krylov_petsc_set_options_database_from_params
(
 krylov_petsc_params_t* input
)
{
  if(input->ksp_monitor)
    PetscOptionsSetValue(NULL,"-ksp_monitor","");
  else
    PetscOptionsClearValue(NULL,"-ksp_monitor");

  if(input->ksp_view)
    PetscOptionsSetValue(NULL,"-ksp_view","");
  else
    PetscOptionsClearValue(NULL,"-ksp_view");

  if(input->ksp_converged_reason)
     PetscOptionsSetValue(NULL,"-ksp_converged_reason","");
  else
    PetscOptionsClearValue(NULL,"-ksp_converged_reason");

  if(input->ksp_initial_guess_nonzero)
    PetscOptionsSetValue(NULL,"-ksp_initial_guess_nonzero","");
  else
    PetscOptionsClearValue(NULL,"-ksp_initial_guess_nonzero");

  PetscOptionsClearValue(NULL,"-ksp_type");
  PetscOptionsSetValue(NULL,"-ksp_type",input->ksp_type);
  
  PetscOptionsClearValue(NULL,"-ksp_atol");
  PetscOptionsSetValue(NULL,"-ksp_atol",input->ksp_atol);
  
  PetscOptionsClearValue(NULL,"-ksp_rtol");
  PetscOptionsSetValue(NULL,"-ksp_rtol",input->ksp_rtol);
  
  PetscOptionsClearValue(NULL,"-ksp_max_it");
  PetscOptionsSetValue(NULL,"-ksp_max_it", input->ksp_max_it);

  if(util_match(input->ksp_type,"chebyshev")){
    PetscOptionsClearValue(NULL,"-ksp_chebyshev_esteig_steps");
    PetscOptionsSetValue(NULL,"-ksp_chebyshev_esteig_steps",input->ksp_chebyshev_esteig_steps);
    PetscOptionsClearValue(NULL,"-ksp_chebyshev_esteig");
    PetscOptionsSetValue(NULL,"-ksp_chebyshev_esteig",input->ksp_chebyshev_esteig);
    if(input->ksp_chebyshev_esteig_random)
      PetscOptionsSetValue(NULL,"-ksp_chebyshev_esteig_random","");
    else
      PetscOptionsClearValue(NULL,"-ksp_chebyshev_esteig_random");
  }
}

void
krylov_petsc_input
(
 p4est_t* p4est,
 const char* input_file,
 const char* input_section,
 const char* printf_prefix,
 krylov_petsc_params_t* input
)
{
  input->ksp_view = -1;
  input->ksp_monitor = -1;
  input->ksp_converged_reason = -1;
  input->ksp_initial_guess_nonzero = -1;
  input->ksp_type[0] = '*';
  input->ksp_atol[0] = '*';
  input->ksp_rtol[0] = '*';
  input->ksp_max_it[0] = '*';
  input->ksp_chebyshev_esteig_steps[0] = '*';
  input->ksp_chebyshev_esteig[0] = '*';
  input->ksp_chebyshev_esteig_random = -1;
   
  mpi_assert(sizeof(input->input_section) <= 50);
  snprintf (input->input_section, sizeof(input->input_section), "%s", input_section);
  if (ini_parse(input_file, krylov_petsc_input_handler, input) < 0) {
    mpi_abort("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, input->ksp_view, -1);
  D4EST_CHECK_INPUT(input_section, input->ksp_monitor, -1);
  D4EST_CHECK_INPUT(input_section, input->ksp_converged_reason, -1);
  D4EST_CHECK_INPUT(input_section, input->ksp_initial_guess_nonzero, -1);
  D4EST_CHECK_INPUT(input_section, input->ksp_type[0], '*');
  D4EST_CHECK_INPUT(input_section, input->ksp_atol[0], '*');
  D4EST_CHECK_INPUT(input_section, input->ksp_rtol[0], '*');
  D4EST_CHECK_INPUT(input_section, input->ksp_max_it[0], '*');

  if(util_match(input->ksp_type,"chebyshev")){
    D4EST_CHECK_INPUT(input_section, input->ksp_chebyshev_esteig_steps[0], '*');
    D4EST_CHECK_INPUT(input_section, input->ksp_chebyshev_esteig[0], '*');
    D4EST_CHECK_INPUT(input_section, input->ksp_chebyshev_esteig_random, -1);
  }
    
  if(p4est->mpirank == 0){
    printf("%s: ksp_type = %s\n",printf_prefix, input->ksp_type);
    printf("%s: ksp_view = %d\n",printf_prefix, input->ksp_view);
    printf("%s: ksp_monitor = %d\n",printf_prefix, input->ksp_monitor);
    printf("%s: ksp_atol = %s\n",printf_prefix, input->ksp_atol);
    printf("%s: ksp_rtol = %s\n",printf_prefix, input->ksp_rtol);
    printf("%s: ksp_maxit = %s\n",printf_prefix, input->ksp_max_it);
    printf("%s: ksp_converged_reason = %d\n",printf_prefix, input->ksp_converged_reason);
    printf("%s: ksp_initial_guess_nonzero = %d\n",printf_prefix, input->ksp_initial_guess_nonzero);
    if(util_match(input->ksp_type,"chebyshev")){
      printf("%s: ksp_chebyshev_esteig_steps = %s\n",printf_prefix, input->ksp_chebyshev_esteig_steps);
      printf("%s: ksp_chebyshev_esteig = %s\n",printf_prefix, input->ksp_chebyshev_esteig);
      printf("%s: ksp_chebyshev_esteig_random = %d\n",printf_prefix, input->ksp_chebyshev_esteig_random);
    }
  }

  
}

static
PetscErrorCode krylov_petsc_apply_aij( Mat A, Vec x, Vec y )
{
  void           *ctx;
  PetscErrorCode ierr;

  petsc_ctx_t* petsc_ctx;
  const double* px;
  double* py;

  /* PetscFunctionBegin; */
  ierr = MatShellGetContext( A, &ctx ); CHKERRQ(ierr);  
  petsc_ctx = (petsc_ctx_t *)ctx;
  ierr = VecGetArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecGetArray( y, &py ); CHKERRQ(ierr);

  /* weakeqn_ptrs_t* fcns = petsc_ctx->fcns; */
  p4est_t* p4est = petsc_ctx->p4est;
  p4est_ghost_t* ghost = *petsc_ctx->ghost;
  dgmath_jit_dbase_t* dgmath_jit_dbase = petsc_ctx->dgmath_jit_dbase;

  problem_data_t vecs_for_aij;
  problem_data_copy_ptrs(petsc_ctx->vecs, &vecs_for_aij);

  vecs_for_aij.u = (double*)px;
  vecs_for_aij.Au = py;

  ((weakeqn_ptrs_t*)(petsc_ctx->fcns))->apply_lhs(p4est,
                                                  ghost,
                                                  (*petsc_ctx->ghost_data),
                                                  &vecs_for_aij,
                                                  dgmath_jit_dbase,
                                                  petsc_ctx->d4est_geom);

  ierr = VecRestoreArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecRestoreArray( y, &py ); CHKERRQ(ierr);
  return ierr;
}

krylov_info_t
krylov_petsc_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 void* fcns,
 p4est_ghost_t** ghost,
 void** ghost_data, 
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geom,
 krylov_petsc_params_t* krylov_petsc_params,
 krylov_pc_t* krylov_pc
)
{
  krylov_petsc_set_options_database_from_params(krylov_petsc_params);
  /* PetscOptionsSetFromOptions(petsc_options); */
  /* krylov_petsc_params_t krylov_params = krylov_petsc_input(p4est,input_file); */
  /* krylov_params.pc = krylov_pc; */
  
  krylov_info_t info;
  KSP ksp;
  Vec x,b;
  PC             pc;
  /* double* u_temp; */
  /* double* rhs_temp; */

  petsc_ctx_t petsc_ctx;
  petsc_ctx.p4est = p4est;
  petsc_ctx.vecs = vecs;
  petsc_ctx.fcns = fcns;
  petsc_ctx.ghost = ghost;
  petsc_ctx.ghost_data = ghost_data;
  petsc_ctx.dgmath_jit_dbase = dgmath_jit_dbase;
  petsc_ctx.d4est_geom = d4est_geom;

  int local_nodes = vecs->local_nodes;
  double* u = vecs->u;
  double* rhs = vecs->rhs;

  KSPCreate(PETSC_COMM_WORLD,&ksp);  
  VecCreate(PETSC_COMM_WORLD,&x);//CHKERRQ(ierr);
  VecSetSizes(x, local_nodes, PETSC_DECIDE);//CHKERRQ(ierr);
  VecSetFromOptions(x);//CHKERRQ(ierr);
  VecDuplicate(x,&b);//CHKERRQ(ierr);

  KSPGetPC(ksp,&pc);
  if (krylov_pc != NULL) {
    PCSetType(pc,PCSHELL);//CHKERRQ(ierr);
    krylov_pc->pc_ctx = &petsc_ctx;
    PCShellSetApply(pc, krylov_petsc_pc_apply);//CHKERRQ(ierr);
    if (krylov_pc->pc_setup != NULL){
      PCShellSetSetUp(pc, krylov_petsc_pc_setup);
    }
    PCShellSetContext(pc, krylov_pc);//CHKERRQ(ierr);
  }
  else {
    PCSetType(pc,PCNONE);//CHKERRQ(ierr);
  }

  KSPSetFromOptions(ksp);

  /* Create matrix-free shell for Aij */
  Mat A;
  MatCreateShell
    (
     PETSC_COMM_WORLD,
     local_nodes,
     local_nodes,
     PETSC_DETERMINE,
     PETSC_DETERMINE,
     (void*)&petsc_ctx,
     &A
    ); 
  MatShellSetOperation(A,MATOP_MULT,(void(*)())krylov_petsc_apply_aij);

  /* Set Amat and Pmat, where Pmat is the matrix the Preconditioner needs */
  KSPSetOperators(ksp,A,A);
  VecPlaceArray(b, rhs);
  VecPlaceArray(x, u);
  
  KSPSolve(ksp,b,x);
  
  KSPGetIterationNumber(ksp, &(info.iterations));
  KSPGetResidualNorm(ksp, &(info.residual_norm));
  
  VecResetArray(b);
  VecResetArray(x);
  VecDestroy(&x);
  VecDestroy(&b);
  KSPDestroy(&ksp);

  return info;
}
