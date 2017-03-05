#include <sc_reduce.h>
#include <pXest.h>
#include <util.h>
#include <linalg.h>
#include <curved_element_data.h>
#include <sipg_flux_vector_fcns.h>
#include <curved_Gauss_primal_sipg_flux_fcns.h>
#include <problem.h>
#include <problem_data.h>
#include <problem_weakeqn_ptrs.h>
#include <central_flux_params.h>
#include <curved_bi_estimator.h>
#include <krylov_petsc.h>
#include <matrix_sym_tester.h>
#include <dg_norm.h>
#include <hp_amr.h>
#include <hp_amr_curved_smooth_pred.h>
#include <d4est_geometry.h>
#include <d4est_geometry_sphere.h>
#include <d4est_geometry_disk.h>
#include <curved_poisson_debug_vecs.h>
#include <d4est_vtk.h>
#include <bi_estimator_flux_fcns.h>
#include <newton_petsc.h>
#include <ini.h>
#include <curved_poisson_operator_primal.h>
#include <curved_Gauss_central_flux_vector_fcns.h>
#include <jacobian_tester.h>
#include "./twopuncturesfcns.h"
#include "./twopuncturesfcns_cactus.h"
#include "time.h"
#include "util.h"

/* soon to be in the input files */
static const double pi = 3.1415926535897932384626433832795;


static
int
amr_mark_element
(
 p4est_t* p4est,
 double eta2,
 estimator_stats_t** stats,
 curved_element_data_t* elem_data,
 void* user
)
{
  int elem_bin;

  /* outer shell */
  if (elem_data->tree < 6){
    elem_bin = 0;
  }
  /* inner shell */
  else if(elem_data->tree < 12){
    elem_bin = 1;
  }
  /* center cube */
  else {
    elem_bin = 2;
  }
 
  /* if (elem_data->tree == 12){ */
  double eta2_percentile
    = estimator_stats_get_percentile(stats[elem_bin], 5);

  /* if (elem_data->tree == 12) */
    return (eta2 >= eta2_percentile);
  /* else */
    /* return 0; */
  /* } */
  /* else { */
    /* return 0; */
  /* } */
}

static
gamma_params_t
amr_set_element_gamma
(
 p4est_t* p4est,
 double eta2,
 estimator_stats_t** stats,
 curved_element_data_t* elem_data,
 void* user
)
{
  gamma_params_t gamma_hpn;
  gamma_hpn.gamma_h = 0;
  gamma_hpn.gamma_p = 0;
  gamma_hpn.gamma_n = 0;

  return gamma_hpn;
}

/* static */
/* double psi_fcn */
/* ( */
/*  double x, */
/*  double y, */
/*  double z, */
/*  double u */
/* ) */
/* { */
/*   double sumn_mn_o_2rn = 0.; */
/*   double dxn, dyn, dzn, r; */
/*   int n; */
/*   for (n = 0; n < NUM_PUNCTURES; n++){ */
/*     dxn = x - xyz_bh[n][0]; */
/*     dyn = y - xyz_bh[n][1]; */
/*     dzn = z - xyz_bh[n][2]; */
/*     r = sqrt(dxn*dxn + dyn*dyn + dzn*dzn); */
/*     if (r == 0.){ */
/*       r += puncture_eps; */
/*     } */
/*     sumn_mn_o_2rn += M_bh[n]/(2.*r); */
/*   } */

/*   return 1. + u + sumn_mn_o_2rn; */
/* } */



static int
uni_refine_function
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t *quadrant
)
{
  /* if (which_tree != 12) */
  return 1;
  /* else */
    /* return 0; */
}



typedef struct {

  int num_unifrefs;
  int num_of_amr_levels;

  int deg_R0;
  int deg_integ_R0;
  int deg_R1;
  int deg_integ_R1;
  int deg_R2;
  int deg_integ_R2;
  int deg_offset_for_nonlinear_integ;
  double ip_flux_penalty;
  int use_cactus;
  int count;
  
} problem_input_t;


static
int problem_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  problem_input_t* pconfig = (problem_input_t*)user;
  if (util_match_couple(section,"amr",name,"num_of_amr_levels")) {
    mpi_assert(pconfig->num_of_amr_levels == -1);
    pconfig->num_of_amr_levels = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"amr",name,"num_unifrefs")) {
    mpi_assert(pconfig->num_unifrefs == -1);
    pconfig->num_unifrefs = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"flux",name,"ip_flux_penalty")) {
    mpi_assert(pconfig->ip_flux_penalty == -1);
    pconfig->ip_flux_penalty = atof(value);
    pconfig->count += 1;
  } 
  else if (util_match_couple(section,"problem",name,"deg_R0")) {
    mpi_assert(pconfig->deg_R0 == -1);
    pconfig->deg_R0 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"use_cactus")) {
    mpi_assert(pconfig->use_cactus == -1);
    pconfig->use_cactus = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_integ_R0")) {
    mpi_assert(pconfig->deg_integ_R0 == -1);
    pconfig->deg_integ_R0 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_R1")) {
    mpi_assert(pconfig->deg_R1 == -1);
    pconfig->deg_R1 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_integ_R1")) {
    mpi_assert(pconfig->deg_integ_R1 == -1);
    pconfig->deg_integ_R1 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_R2")) {
    mpi_assert(pconfig->deg_R2 == -1);
    pconfig->deg_R2 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_integ_R2")) {
    mpi_assert(pconfig->deg_integ_R2 == -1);
    pconfig->deg_integ_R2 = atoi(value);
    pconfig->count += 1;
  }  
  else if (util_match_couple(section,"problem",name,"deg_offset_for_nonlinear_integ")) {
    mpi_assert(pconfig->deg_offset_for_nonlinear_integ == -1);
    pconfig->deg_offset_for_nonlinear_integ = atoi(value);
    pconfig->count += 1;
  }  

  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


static
problem_input_t
problem_input
(
 const char* input_file
)
{
  int num_of_options = 11;
  
  problem_input_t input;
  input.num_unifrefs = -1;
  input.num_of_amr_levels = -1;
  input.ip_flux_penalty = -1;
  input.use_cactus = -1;
  input.deg_R0 = -1;
  input.deg_integ_R0 = -1;
  input.deg_R1 = -1;
  input.deg_integ_R1 = -1;
  input.deg_R2 = -1;
  input.deg_integ_R2 = -1; 
  input.deg_offset_for_nonlinear_integ = -1;
  
  input.count = 0;
  
  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  mpi_assert(input.count == num_of_options);
  return input;
}

/* p4est_geometry_t* */
/* problem_build_geom */
/* ( */
/*  p4est_connectivity_t* conn */
/* ) */
/* { */
/*   /\* mpi_assert((P4EST_DIM)==3); *\/ */
/*   p4est_geometry_t* geom; */
/*   problem_input_t input = problem_input("options.input"); */

/*   /\* geom = d4est_geometry_new_compact_sphere(conn, input.R2, input.R1, input.R0, input.w, input.Rinf); *\/ */
/* #if (P4EST_DIM)==3 */
/*   geom = d4est_geometry_new_sphere(conn, input.R2, input.R1, input.R0); */
/* #endif */
/* #if (P4EST_DIM)==2 */
/*   geom = d4est_geometry_new_disk(conn, input.R1, input.R2); */
/* #endif */
  
/*   /\* printf("input.R2 = %.25f\n", input.R2); *\/ */

  
/*   return geom;   */
/* } */

p4est_t*
problem_build_p4est
(
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t* conn,
 p4est_locidx_t min_quadrants,
 int min_level, 
 int fill_uniform
)
{
  return p4est_new_ext
    (
     mpicomm,
     conn,
     min_quadrants,
     min_level,
     fill_uniform,
     sizeof(curved_element_data_t),
     NULL,
     NULL
    );
}

p4est_t*
problem_load_p4est_from_checkpoint
(
 const char* filename,
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t** conn
){
  int autopartition = 1;
  int load_data = 1;
  int broadcasthead = 0;
  
  return p4est_load_ext (filename,
                mpicomm,
                sizeof(curved_element_data_t),
                load_data,
                autopartition,
                broadcasthead,
                NULL,
                conn);
}


int
in_bin_fcn
(
 curved_element_data_t* elem_data,
 int bin
)
{
  int elem_bin;

  /* outer shell */
  if (elem_data->tree < 6){
    elem_bin = 0;
  }
  /* inner shell */
  else if(elem_data->tree < 12){
    elem_bin = 1;
  }
  /* center cube */
  else {
    elem_bin = 2;
  }

  return (elem_bin == bin);
}



void
problem_set_degrees
(
 curved_element_data_t* elem_data,
 void* user_ctx
)
{
  problem_input_t* input = user_ctx;
  /* outer shell */
  if (elem_data->tree < 6){
    elem_data->deg = input->deg_R2;
    elem_data->deg_integ = input->deg_integ_R2;
  }
  /* inner shell */
  else if(elem_data->tree < 12){
    elem_data->deg = input->deg_R1;
    elem_data->deg_integ = input->deg_integ_R1;
  }
  /* center cube */
  else {
    elem_data->deg = input->deg_R0;
    elem_data->deg_integ = input->deg_integ_R0;
  } 
}


void
problem_save_to_vtk
(
 p4est_t* p4est,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* u,
 int level,
 int with_eta,
 double R0
)
{
   int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
   double* eta_array = P4EST_ALLOC(double, p4est->local_num_quadrants);
   int vtk_nodes = 0;
     
    int stride = 0;
    for (p4est_topidx_t tt = p4est->first_local_tree;
         tt <= p4est->last_local_tree;
         ++tt)
      {
        p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
        sc_array_t* tquadrants = &tree->quadrants;
        int Q = (p4est_locidx_t) tquadrants->elem_count;
        for (int q = 0; q < Q; ++q) {
          /* k++; */
          p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
          curved_element_data_t* ed = quad->p.user_data;
          deg_array[stride] = ed->deg;
          eta_array[stride] = ed->local_estimator;
          vtk_nodes = util_int_pow_int(deg_array[stride], (P4EST_DIM))*(P4EST_CHILDREN);
          stride++;
        }
      }


    p4est_connectivity_t* conn_vtk = p8est_connectivity_new_sphere();
    p4est_geometry_t* geom_vtk = d4est_geometry_compactified_sphere_from_param
                                 (
                                  R0,
                                  2.*R0,
                                  3.*R0,
                                  conn_vtk
                                 );


    char sol_save_as [500];
    if (with_eta)
      sprintf(sol_save_as, "%s_hp_amr_level_%d_sols_witheta", "puncture", level);
    else
      sprintf(sol_save_as, "%s_hp_amr_level_%d_sols_noeta", "puncture", level);
    
    d4est_vtk_context_t* vtk_ctx = d4est_vtk_dg_context_new(p4est, dgmath_jit_dbase, sol_save_as);
    d4est_vtk_context_set_geom(vtk_ctx, geom_vtk);
    d4est_vtk_context_set_scale(vtk_ctx, .99);
    d4est_vtk_context_set_deg_array(vtk_ctx, deg_array);
    vtk_ctx = d4est_vtk_write_dg_header(vtk_ctx, dgmath_jit_dbase);    
    vtk_ctx = d4est_vtk_write_dg_point_dataf(vtk_ctx,
                                             1,
                                             0,
                                             "u",
                                             u,
                                             vtk_ctx
                                            );

    if (with_eta){
      vtk_ctx = d4est_vtk_write_dg_cell_dataf
                (
                 vtk_ctx,
                 1,
                 1,
                 1,
                 0,
                 1,
                 1,
                 0,
                 "eta",
                 eta_array,
                 vtk_ctx
                );
    }
    else {
      vtk_ctx = d4est_vtk_write_dg_cell_dataf
                (
                 vtk_ctx,
                 1,
                 1,
                 1,
                 0,
                 1,
                 0,
                 0,
                 vtk_ctx
                );
    }
    
    d4est_vtk_write_footer(vtk_ctx);
    P4EST_FREE(deg_array);
    P4EST_FREE(eta_array);
    p8est_connectivity_destroy(conn_vtk);
    p8est_geometry_destroy(geom_vtk);
}


void
problem_init
(
 const char* input_file,
 p4est_t* p4est,
 d4est_geometry_t* d4est_geom,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 int proc_size,
 sc_MPI_Comm mpicomm
)
{


  problem_input_t input = problem_input(input_file);
  
  int level;
  
  mpi_assert((P4EST_DIM) == 2 || (P4EST_DIM) == 3);
  int world_rank, world_size;
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &world_rank);
  sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &world_size);

  double* Au = P4EST_ALLOC_ZERO(double, 1);
  double* rhs = P4EST_ALLOC_ZERO(double, 1);
  double* u = P4EST_ALLOC_ZERO(double, 1);
  double* u_prev = P4EST_ALLOC_ZERO(double, 1);
  int local_nodes = 1;

  penalty_calc_t bi_u_penalty_fcn = bi_u_prefactor_conforming_maxp_minh;
  penalty_calc_t bi_u_dirichlet_penalty_fcn = bi_u_prefactor_conforming_maxp_minh;
  penalty_calc_t bi_gradu_penalty_fcn = bi_gradu_prefactor_maxp_minh;
  
  ip_flux_params_t ip_flux_params;
  ip_flux_params.ip_flux_penalty_prefactor = input.ip_flux_penalty;
  ip_flux_params.ip_flux_penalty_calculate_fcn = sipg_flux_vector_calc_penalty_maxp2_over_minh;

  central_flux_params_t central_flux_params;
  central_flux_params.central_flux_penalty_prefactor = input.ip_flux_penalty;
  
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  /* create space for storing the ghost data */
  curved_element_data_t* ghost_data = P4EST_ALLOC (curved_element_data_t,
                                                   ghost->ghosts.elem_count);

  p4est_partition(p4est, 0, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FACE, NULL);
  /* geometric_factors_t* geometric_factors = geometric_factors_init(p4est); */


  /* grid_fcn_t boundary_flux_fcn = zero_fcn; */
  twopunctures_params_t tp_params;
  init_random_puncture_data(p4est, &tp_params, input.deg_offset_for_nonlinear_integ);
  /* init_S_puncture_data(p4est, &tp_params, input.deg_offset_for_nonlinear_integ); */

  twopunctures_cactus_params_t tp_cactus_params;
  init_cactus_puncture_data(&tp_cactus_params, input.deg_offset_for_nonlinear_integ);

  
  problem_data_t prob_vecs;
  prob_vecs.rhs = rhs;
  prob_vecs.Au = Au;
  prob_vecs.u = u;
  prob_vecs.local_nodes = local_nodes;

  prob_vecs.curved_scalar_flux_fcn_data = curved_Gauss_primal_sipg_flux_dirichlet_fetch_fcns
                                          (zero_fcn,&ip_flux_params);

  if(input.use_cactus){
    prob_vecs.user = &tp_cactus_params;
  }
  else {
    prob_vecs.user = &tp_params;
  }

  curved_weakeqn_ptrs_t prob_fcns;


  if(input.use_cactus){
    prob_fcns.build_residual = build_residual_cactus;
    prob_fcns.apply_lhs = apply_jac_cactus;
  }
  else {
    prob_fcns.build_residual = build_residual;
    prob_fcns.apply_lhs = apply_jac;
  }
  
  geometric_factors_t* geometric_factors = geometric_factors_init(p4est);

  /* printf("[D4EST_INFO]: Initial number of elements = %d\n", p4est->local_num_quadrants); */

  for (level = 0; level < input.num_unifrefs; ++level){

    if (world_rank == 0)
      printf("[D4EST_INFO]: UNIFORM REFINEMENT LEVEL %d\n", level);

    
      p4est_refine_ext(p4est,
                       0,
                       -1,
                       uni_refine_function,
                       NULL,
                       NULL
                      );

      p4est_partition(p4est, 0, NULL);
      p4est_balance_ext
        (
         p4est,
         P4EST_CONNECT_FACE,
         NULL,
         NULL
        );

      p4est_ghost_destroy(ghost);
      P4EST_FREE(ghost_data);

      ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
      ghost_data = P4EST_ALLOC(curved_element_data_t, ghost->ghosts.elem_count);




      curved_element_data_init_new(p4est,
                                   geometric_factors,
                                   dgmath_jit_dbase,
                                   d4est_geom,
                                   problem_set_degrees,
                                   (void*)&input);





    
      local_nodes = curved_element_data_get_local_nodes(p4est);

      Au = P4EST_REALLOC(Au, double, local_nodes);
      u = P4EST_REALLOC(u, double, local_nodes);
      rhs = P4EST_REALLOC(rhs, double, local_nodes);
      prob_vecs.Au = Au;
      prob_vecs.u = u;
      prob_vecs.rhs = rhs;
      prob_vecs.local_nodes = local_nodes;
      
      /* curved_bi_estimator_compute */
      /*   ( */
      /*    p4est, */
      /*    &prob_vecs, */
      /*    &prob_fcns, */
      /*    bi_u_penalty_fcn, */
      /*    bi_u_dirichlet_penalty_fcn, */
      /*    bi_gradu_penalty_fcn, */
      /*    zero_fcn, */
      /*    ip_flux_params.ip_flux_penalty_prefactor, */
      /*    ghost, */
      /*    ghost_data, */
      /*    dgmath_jit_dbase, */
      /*    d4est_geom */
      /*   ); */

      
      
      
      /* estimator_stats_t stats; */
      /* estimator_stats_compute(p4est, &stats,1); */

      /* /\* if(world_rank == 0) *\/ */
      /* estimator_stats_print(&stats); */
      


  }


  /* linalg_fill_vec(prob_vecs.u, 0., local_nodes); */

  
  d4est_geom->dxdr_method = INTERP_X_ON_LOBATTO;    
  /* curved_element_data_init(p4est, geometric_factors, dgmath_jit_dbase, d4est_geom, degree, input.gauss_integ_deg); */
  curved_element_data_init_new(p4est,
                               geometric_factors,
                               dgmath_jit_dbase,
                               d4est_geom,
                               problem_set_degrees,
                               (void*)&input);


    
    local_nodes = curved_element_data_get_local_nodes(p4est);

    Au = P4EST_REALLOC(Au, double, local_nodes);
    u = P4EST_REALLOC(u, double, local_nodes);
    rhs = P4EST_REALLOC(rhs, double, local_nodes);

    prob_vecs.Au = Au;
    prob_vecs.u = u;
    prob_vecs.rhs = rhs;
    prob_vecs.local_nodes = local_nodes;

    curved_smooth_pred_marker_t amr_marker;
    amr_marker.user = (void*)&input;
    amr_marker.mark_element_fcn = amr_mark_element;
    amr_marker.set_element_gamma_fcn = amr_set_element_gamma;
    amr_marker.name = "puncture_marker";

    hp_amr_scheme_t* scheme =
      hp_amr_curved_smooth_pred_init
      (
       p4est,
       (MAX_DEGREE)-2,
       amr_marker
      );

    /* linalg_fill_vec(prob_vecs.u, 0.001, prob_vecs.local_nodes); */
    linalg_fill_vec(prob_vecs.u, 0, prob_vecs.local_nodes);
    

    
  for (level = 0; level < input.num_of_amr_levels; ++level){

    if (world_rank == 0)
      printf("[D4EST_INFO]: AMR REFINEMENT LEVEL %d\n", level);
    

    curved_bi_estimator_compute
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       bi_u_penalty_fcn,
       bi_u_dirichlet_penalty_fcn,
       bi_gradu_penalty_fcn,
       zero_fcn,
       ip_flux_params.ip_flux_penalty_prefactor,
       ghost,
       ghost_data,
       dgmath_jit_dbase,
       d4est_geom
      );

    
    estimator_stats_t* stats [3];
    for (int i = 0; i < 3; i++){
      stats[i] = P4EST_ALLOC(estimator_stats_t, 1);
    }
    
    double local_eta2 = estimator_stats_compute_per_bin
                        (
                         p4est,
                         &stats[0],
                         3,
                         in_bin_fcn
                        );

    curved_element_data_print_number_of_elements_per_tree(p4est);
    
    estimator_stats_compute_max_percentiles_across_proc
      (
       stats,
       3
      );

    if (world_rank == 0){
      for (int i = 0; i < 3; i++){
        estimator_stats_print(stats[i]);
      }
    }
    

    double R0 = ((d4est_geometry_sphere_attr_t*)(d4est_geom->p4est_geom->user))->R0;
    /* printf("R0 = %.f\n", R0); */

   /*  d4est_geometry_t* d4est_geom_vtk = d4est_geometry */
  /*   p4est_connectivity_t* conn_vtk = p8est_connectivity_new_sphere(); */
  /*   p4est_geometry_t* geom_vtk = p8est_geometry_new_sphere(conn_vtk, R0*3, R0*2, R0); */
  /*   d4est_vtk_context_t* vtk_ctx = d4est_vtk_dg_context_new(p4est, dgmath_jit_dbase, "puncture-sphere"); */
  /*   d4est_vtk_context_set_geom(vtk_ctx, geom_vtk); */
  /*   d4est_vtk_context_set_scale(vtk_ctx, .99); */
  /*   d4est_vtk_context_set_deg_array(vtk_ctx, deg_array); */
  /*   vtk_ctx = d4est_vtk_write_dg_header(vtk_ctx, dgmath_jit_dbase); */
  /*   vtk_ctx = d4est_vtk_write_dg_point_dataf(vtk_ctx, 1, 0, "u",u, vtk_ctx); */
  /*   vtk_ctx = d4est_vtk_write_dg_cell_dataf */
  /*             ( */
  /*              vtk_ctx, */
  /*              1, */
  /*              1, */
  /*              1, */
  /*              0, */
  /*              1, */
  /*              1, */
  /*              0, */
  /*              "eta", */
  /*              eta_array, */
  /*              vtk_ctx */
  /*             ); */


  
  /* d4est_vtk_write_footer(vtk_ctx); */
  /* P4EST_FREE(deg_array); */
  /* P4EST_FREE(eta_array); */
  /* p8est_connectivity_destroy(conn_vtk); */
  /* p8est_geometry_destroy(geom_vtk); */

    problem_save_to_vtk
      (
       p4est,
       dgmath_jit_dbase,
       u,
       level,
       1,
       R0
      );
  
    
    hp_amr(p4est,
           dgmath_jit_dbase,
           &u,
           &stats[0],
           scheme,
           1
          );        

    for (int i = 0; i < 3; i++){
      P4EST_FREE(stats[i]);
    }
    
    p4est_ghost_destroy(ghost);
    P4EST_FREE(ghost_data);

    ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
    ghost_data = P4EST_ALLOC(curved_element_data_t, ghost->ghosts.elem_count);
    
    curved_element_data_init_new(p4est,
                                 geometric_factors,
                                 dgmath_jit_dbase,
                                 d4est_geom,
                                 problem_set_degrees,
                                 (void*)&input);
    
    local_nodes = curved_element_data_get_local_nodes(p4est);

    u_prev = P4EST_REALLOC(u_prev, double, local_nodes);
    Au = P4EST_REALLOC(Au, double, local_nodes);
    prob_vecs.Au = Au;
    prob_vecs.u = u;
    prob_vecs.u0 = u;
    prob_vecs.local_nodes = local_nodes;

    linalg_copy_1st_to_2nd(u, u_prev, local_nodes);

    clock_t begin = 0;
    clock_t end = -1;

    if (world_rank == 0){
      begin = clock();
    }

    /* jacobian_tester */
    /*   ( */
    /*    p4est, */
    /*    ghost, */
    /*    ghost_data, */
    /*    dgmath_jit_dbase, */
    /*    d4est_geom, */
    /*    &prob_fcns, */
    /*    &prob_vecs */
    /*   ); */

    
    newton_petsc_solve
      (
       p4est,
       &prob_vecs,
       (void*)&prob_fcns,
       &ghost,
       (void**)&ghost_data,
       dgmath_jit_dbase,
       d4est_geom,
       input_file,
       NULL,
       NULL,
       NULL
      );

    /* matrix_spd_tester_parallel */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs,  */
    /*    &prob_fcns, */
    /*    ghost, */
    /*    ghost_data, */
    /*    dgmath_jit_dbase, */
    /*    d4est_geom, */
    /*    1, */
    /*    20 */
    /*   ); */

    /* matrix_sym_tester_parallel */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs,  */
    /*    &prob_fcns, */
    /*    ghost, */
    /*    ghost_data, */
    /*    dgmath_jit_dbase, */
    /*    d4est_geom, */
    /*    1, */
    /*    20, */
    /*    .000000001 */
    /*   );  */

    linalg_vec_axpy(-1., prob_vecs.u, u_prev, local_nodes);

    double local_l2_norm_sqr = curved_element_data_compute_l2_norm_sqr
                                (
                                 p4est,
                                 u_prev,
                                 local_nodes,
                                 dgmath_jit_dbase,
                                 DO_NOT_STORE_LOCALLY
                                );
    

    
    double local_nodes_dbl = (double)local_nodes;
    double local_reduce [3];
    double global_reduce [3];

    local_reduce[0] = local_l2_norm_sqr;
    local_reduce[1] = local_nodes_dbl;
    local_reduce[2] = local_eta2;
    
    sc_reduce
      (
       &local_reduce[0],
       &global_reduce[0],
       3,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       0,
       sc_MPI_COMM_WORLD
      );

    double global_l2_norm_sqr = global_reduce[0];
    double global_nodes_dbl = global_reduce[1];
    double global_eta2 = global_reduce[2];

    if (world_rank == 0){
      end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf
      (
       "[HP_AMR]: %d, %d, %d, %.25f, %.25f, %f\n",
       level,
       (int)p4est->global_num_quadrants,
       (int)global_nodes_dbl,
       sqrt(global_eta2),
       sqrt(global_l2_norm_sqr),
       time_spent
      );
    }

  }
         
  geometric_factors_destroy(geometric_factors);

  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }

  P4EST_FREE(Au);
  P4EST_FREE(rhs);
  /* P4EST_FREE(u_analytic); */
  P4EST_FREE(u);
}
