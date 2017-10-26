#include <d4est_util.h>
#include <pXest.h>
#include <d4est_linalg.h>
#include <d4est_elliptic_data.h>
#include <krylov_petsc.h>
#include <petscsnes.h>
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
  if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_atol")) {
    D4EST_ASSERT(pconfig->ksp_atol[0] == '*');
    snprintf (pconfig->ksp_atol, sizeof(pconfig->ksp_atol), "%s", value);
    /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_atol",value); */
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_rtol")) {
    D4EST_ASSERT(pconfig->ksp_rtol[0] == '*');
    snprintf (pconfig->ksp_rtol, sizeof(pconfig->ksp_rtol), "%s", value);

    /* pconfig->ksp_rtol = atof(value); */
    /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_rtol",value); */
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_max_it")) {
    D4EST_ASSERT(pconfig->ksp_max_it[0] == '*');
    snprintf (pconfig->ksp_max_it, sizeof(pconfig->ksp_max_it), "%s", value);
    D4EST_ASSERT(atoi(value) > -1);
    /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_max_it",value); */
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_view")) {
    D4EST_ASSERT(pconfig->ksp_view == -1);
    pconfig->ksp_view = atoi(value);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
    /* if (atoi(value) == 1){ */
      /* PetscOptionsSetValue(pconfig->petsc_options,"-s_view",""); */
    /* } */
  }  
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_monitor")) {
    D4EST_ASSERT(pconfig->ksp_monitor == -1);
    pconfig->ksp_monitor = atoi(value);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
    /* if (atoi(value) == 1){ */
      /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_monitor",""); */
    /* } */
  } 
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_converged_reason")) {
    D4EST_ASSERT(pconfig->ksp_converged_reason == -1);
    pconfig->ksp_converged_reason = atoi(value);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
    /* if (atoi(value) == 1){ */
      /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_converged_reason",""); */
    /* } */
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_monitor_singular_value")) {
    pconfig->ksp_monitor_singular_value = atoi(value);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
  }      
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_initial_guess_nonzero")) {
    D4EST_ASSERT(pconfig->ksp_initial_guess_nonzero == -1);
    pconfig->ksp_initial_guess_nonzero = atoi(value);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
    /* if (atoi(value) == 1){ */
      /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_initial_guess_nonzero",""); */
    /* } */
  }  
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_type")) {
    D4EST_ASSERT(pconfig->ksp_type[0] == '*');
    snprintf (pconfig->ksp_type, sizeof(pconfig->ksp_type), "%s", value);
    /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_type",value); */
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_chebyshev_esteig_steps")) {
    D4EST_ASSERT(pconfig->ksp_chebyshev_esteig_steps[0] == '*');
    snprintf (pconfig->ksp_chebyshev_esteig_steps, sizeof(pconfig->ksp_chebyshev_esteig_steps), "%s", value);
    /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_chebyshev_esteig_steps",value); */
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_chebyshev_esteig")) {
    D4EST_ASSERT(pconfig->ksp_chebyshev_esteig[0] == '*');
    snprintf (pconfig->ksp_chebyshev_esteig, sizeof(pconfig->ksp_chebyshev_esteig), "%s", value);
    /* PetscOptionsSetValue(pconfig->petsc_options,"-ksp_chebyshev_esteig",value); */
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"ksp_chebyshev_esteig_random")) {
    D4EST_ASSERT(pconfig->ksp_chebyshev_esteig_random == -1);
    pconfig->ksp_chebyshev_esteig_random = atoi(value);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
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

  if(input->ksp_monitor_singular_value == 1)
    PetscOptionsSetValue(NULL,"-ksp_monitor_singular_value","");
  else
    PetscOptionsClearValue(NULL,"-ksp_monitor_singular_value");
  
  PetscOptionsClearValue(NULL,"-ksp_type");
  PetscOptionsSetValue(NULL,"-ksp_type",input->ksp_type);
  
  PetscOptionsClearValue(NULL,"-ksp_atol");
  PetscOptionsSetValue(NULL,"-ksp_atol",input->ksp_atol);
  
  PetscOptionsClearValue(NULL,"-ksp_rtol");
  PetscOptionsSetValue(NULL,"-ksp_rtol",input->ksp_rtol);
  
  PetscOptionsClearValue(NULL,"-ksp_max_it");
  PetscOptionsSetValue(NULL,"-ksp_max_it", input->ksp_max_it);

  if(d4est_util_match(input->ksp_type,"chebyshev")){
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
  input->ksp_monitor_singular_value = 0;
  
  D4EST_ASSERT(sizeof(input->input_section) <= 50);
  snprintf (input->input_section, sizeof(input->input_section), "%s", input_section);
  if (ini_parse(input_file, krylov_petsc_input_handler, input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, input->ksp_view, -1);
  D4EST_CHECK_INPUT(input_section, input->ksp_monitor, -1);
  D4EST_CHECK_INPUT(input_section, input->ksp_converged_reason, -1);
  D4EST_CHECK_INPUT(input_section, input->ksp_initial_guess_nonzero, -1);
  D4EST_CHECK_INPUT(input_section, input->ksp_type[0], '*');
  D4EST_CHECK_INPUT(input_section, input->ksp_atol[0], '*');
  D4EST_CHECK_INPUT(input_section, input->ksp_rtol[0], '*');
  D4EST_CHECK_INPUT(input_section, input->ksp_max_it[0], '*');

  if(d4est_util_match(input->ksp_type,"chebyshev")){
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
    if(d4est_util_match(input->ksp_type,"chebyshev")){
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

  /* d4est_elliptic_eqns_t* fcns = petsc_ctx->fcns; */
  p4est_t* p4est = petsc_ctx->p4est;
  p4est_ghost_t* ghost = *petsc_ctx->ghost;

  d4est_elliptic_data_t vecs_for_aij;
  d4est_elliptic_data_copy_ptrs(petsc_ctx->vecs, &vecs_for_aij);

  vecs_for_aij.u = (double*)px;
  vecs_for_aij.Au = py;
 
  d4est_elliptic_eqns_apply_lhs
    (
     petsc_ctx->p4est,
     *petsc_ctx->ghost,
     *petsc_ctx->ghost_data,
     petsc_ctx->fcns,
     &vecs_for_aij,
     petsc_ctx->d4est_ops,
     petsc_ctx->d4est_geom,
     petsc_ctx->d4est_quad,
     petsc_ctx->d4est_factors
    );
  
  ierr = VecRestoreArrayRead( x, &px ); CHKERRQ(ierr);
  ierr = VecRestoreArray( y, &py ); CHKERRQ(ierr);
  return ierr;
}

krylov_info_t
krylov_petsc_solve
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 p4est_ghost_t** ghost,
 d4est_element_data_t** ghost_data, 
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_geometry_storage_t* d4est_factors,
 krylov_petsc_params_t* krylov_petsc_params,
 krylov_pc_t* krylov_pc
)
{
  krylov_petsc_set_options_database_from_params(krylov_petsc_params);

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
  petsc_ctx.d4est_ops = d4est_ops;
  petsc_ctx.d4est_geom = d4est_geom;
  petsc_ctx.d4est_quad = d4est_quad;

  int local_nodes = vecs->local_nodes;
  double* u = vecs->u;
  double* rhs = vecs->rhs;

  /* printf("PRE KRYLOV SOLVE REDUCTIONS\n"); */
  /* DEBUG_PRINT_ARR_DBL_SUM(u, local_nodes); */
  /* DEBUG_PRINT_ARR_DBL_SUM(rhs, local_nodes); */
  
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
  /* KSPSetUp(ksp); */
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

  MatDestroy(&A);
  VecResetArray(b);
  VecResetArray(x);
  VecDestroy(&x);
  VecDestroy(&b);
  KSPDestroy(&ksp);

  return info;
}
