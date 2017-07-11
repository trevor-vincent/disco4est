#include <sc_reduce.h>
#include <pXest.h>
#include <util.h>
#include <d4est_linalg.h>
#include <d4est_element_data.h>
#include <sipg_flux_vector_fcns.h>
#include <curved_gauss_primal_sipg_flux_fcns.h>
#include <problem.h>
#include <problem_data.h>
#include <d4est_elliptic_eqns.h>
#include <central_flux_params.h>
#include <curved_bi_estimator.h>
#include <krylov_petsc.h>
#include <matrix_sym_tester.h>
#include <dg_norm.h>
#include <hp_amr.h>
#include <hp_amr_curved_smooth_pred.h>
#include <d4est_geometry.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_geometry_disk.h>
#include <curved_poisson_debug_vecs.h>
#include <d4est_vtk.h>
#include <bi_estimator_flux_fcns.h>
#include <newton_petsc.h>
#include <ini.h>
#include <d4est_poisson.h>
#include <curved_gauss_central_flux_vector_fcns.h>
#include <jacobian_tester.h>
#include "./twopuncturesfcns.h"
#include "./twopuncturesfcns_cactus.h"
#include "./twopuncturesfcns_spec.h"
#include <multigrid_smoother_cheby_d4est.h>
#include <multigrid_smoother_krylov_petsc.h>
#include <multigrid_bottom_solver_cg_d4est.h>
#include <multigrid_bottom_solver_cheby_d4est.h>
#include <multigrid_bottom_solver_krylov_petsc.h>
#include <multigrid_logger_residual.h>
#include <multigrid_element_data_updater_curved.h>
#include <krylov_pc_multigrid.h>
#include "time.h"
#include "util.h"
#include <newton_d4est.h>
#include <newton_petsc.h>

/* soon to be in the input files */
static const double pi = 3.1415926535897932384626433832795;




static
int
amr_mark_element
(
 p4est_t* p4est,
 double eta2,
 estimator_stats_t** stats,
 d4est_element_data_t* elem_data,
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
 d4est_element_data_t* elem_data,
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
  int deg_quad_R0;
  int deg_R1;
  int deg_quad_R1;
  int deg_R2;
  int deg_quad_R2;
  int deg_offset_for_nonlinear_quad;
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
    D4EST_ASSERT(pconfig->num_of_amr_levels == -1);
    pconfig->num_of_amr_levels = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"amr",name,"num_unifrefs")) {
    D4EST_ASSERT(pconfig->num_unifrefs == -1);
    pconfig->num_unifrefs = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"flux",name,"ip_flux_penalty")) {
    D4EST_ASSERT(pconfig->ip_flux_penalty == -1);
    pconfig->ip_flux_penalty = atof(value);
    pconfig->count += 1;
  } 
  else if (util_match_couple(section,"problem",name,"deg_R0")) {
    D4EST_ASSERT(pconfig->deg_R0 == -1);
    pconfig->deg_R0 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"use_cactus")) {
    D4EST_ASSERT(pconfig->use_cactus == -1);
    pconfig->use_cactus = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_quad_R0")) {
    D4EST_ASSERT(pconfig->deg_quad_R0 == -1);
    pconfig->deg_quad_R0 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_R1")) {
    D4EST_ASSERT(pconfig->deg_R1 == -1);
    pconfig->deg_R1 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_quad_R1")) {
    D4EST_ASSERT(pconfig->deg_quad_R1 == -1);
    pconfig->deg_quad_R1 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_R2")) {
    D4EST_ASSERT(pconfig->deg_R2 == -1);
    pconfig->deg_R2 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_quad_R2")) {
    D4EST_ASSERT(pconfig->deg_quad_R2 == -1);
    pconfig->deg_quad_R2 = atoi(value);
    pconfig->count += 1;
  }  
  else if (util_match_couple(section,"problem",name,"deg_offset_for_nonlinear_quad")) {
    D4EST_ASSERT(pconfig->deg_offset_for_nonlinear_quad == -1);
    pconfig->deg_offset_for_nonlinear_quad = atoi(value);
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
  input.deg_quad_R0 = -1;
  input.deg_R1 = -1;
  input.deg_quad_R1 = -1;
  input.deg_R2 = -1;
  input.deg_quad_R2 = -1; 
  input.deg_offset_for_nonlinear_quad = -1;
  
  input.count = 0;
  
  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_ASSERT(input.count == num_of_options);
  return input;
}

void
set_deg_quad
(
 void* elem_data_tmp,
 void* user_ctx
)
{
  d4est_element_data_t* elem_data = elem_data_tmp;
  problem_input_t* input = user_ctx;
  /* outer shell */
  if (elem_data->tree < 6){
    elem_data->deg_quad = input->deg_quad_R2;
    elem_data->deg_stiffness = input->deg_quad_R2;
  }
  /* inner shell */
  else if(elem_data->tree < 12){
    elem_data->deg_quad = input->deg_quad_R1;
    elem_data->deg_stiffness = input->deg_quad_R1;
  }
  /* center cube */
  else {
    elem_data->deg_quad = input->deg_quad_R0;
    elem_data->deg_stiffness = input->deg_quad_R0;
  } 
}
 
/* p4est_geometry_t* */
/* problem_build_geom */
/* ( */
/*  p4est_connectivity_t* conn */
/* ) */
/* { */
/*   /\* D4EST_ASSERT((P4EST_DIM)==3); *\/ */
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
     sizeof(d4est_element_data_t),
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
                sizeof(d4est_element_data_t),
                load_data,
                autopartition,
                broadcasthead,
                NULL,
                conn);
}


int
in_bin_fcn
(
 d4est_element_data_t* elem_data,
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
 void* elem_data_tmp,
 void* user_ctx
)
{
  d4est_element_data_t* elem_data = elem_data_tmp;
  problem_input_t* input = user_ctx;
  /* outer shell */
  if (elem_data->tree < 6){
    elem_data->deg = input->deg_R2;
    elem_data->deg_quad = input->deg_quad_R2;
    elem_data->deg_stiffness = input->deg_quad_R2;
  }
  /* inner shell */
  else if(elem_data->tree < 12){
    elem_data->deg = input->deg_R1;
    elem_data->deg_quad = input->deg_quad_R1;
        elem_data->deg_stiffness = input->deg_quad_R1;
  }
  /* center cube */
  else {
    elem_data->deg = input->deg_R0;
    elem_data->deg_quad = input->deg_quad_R0;
    elem_data->deg_stiffness = input->deg_quad_R0;
  } 
}


void
problem_save_to_vtk
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 double* u,
 int level,
 int with_eta,
 double R0,
 double R1,
 double R2,
 int compactify_outer_shell,
 int compactify_inner_shell,
 const char* input_file
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
          d4est_element_data_t* ed = quad->p.user_data;
          deg_array[stride] = ed->deg;
          eta_array[stride] = ed->local_estimator;
          vtk_nodes = util_int_pow_int(deg_array[stride], (P4EST_DIM))*(P4EST_CHILDREN);
          stride++;
        }
      }


     d4est_geometry_t* geom_vtk = d4est_geometry_new
                                 (
                                  p4est->mpirank,
                                  input_file
                                 );


    ((d4est_geometry_cubed_sphere_attr_t*)geom_vtk->user)->R2 = R2;
    ((d4est_geometry_cubed_sphere_attr_t*)geom_vtk->user)->R1 = R1;
    ((d4est_geometry_cubed_sphere_attr_t*)geom_vtk->user)->R0 = R0;
    ((d4est_geometry_cubed_sphere_attr_t*)geom_vtk->user)->compactify_outer_shell = compactify_outer_shell;
    ((d4est_geometry_cubed_sphere_attr_t*)geom_vtk->user)->compactify_inner_shell = compactify_inner_shell;

    char sol_save_as [500];
    if (with_eta)
      sprintf(sol_save_as, "%s_hp_amr_level_%d_sols_witheta", "puncture", level);
    else
      sprintf(sol_save_as, "%s_hp_amr_level_%d_sols_noeta", "puncture", level);
    
    d4est_vtk_context_t* vtk_ctx = d4est_vtk_dg_context_new(p4est, d4est_ops, sol_save_as);
    d4est_vtk_context_set_geom(vtk_ctx, geom_vtk);
    d4est_vtk_context_set_scale(vtk_ctx, .99);
    d4est_vtk_context_set_deg_array(vtk_ctx, deg_array);
    vtk_ctx = d4est_vtk_write_dg_header(vtk_ctx, d4est_ops);    
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
    d4est_geometry_destroy(geom_vtk);
}


void
problem_init
(
 const char* input_file,
 p4est_t* p4est,
 d4est_geometry_t* d4est_geom,
 d4est_operators_t* d4est_ops,
 int proc_size,
 sc_MPI_Comm mpicomm
)
{


  problem_input_t input = problem_input(input_file);
  
  int level;
  
  D4EST_ASSERT((P4EST_DIM) == 2 || (P4EST_DIM) == 3);
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
  
  d4est_poisson_flux_sipg_params_t ip_flux_params;
  ip_flux_params.sipg_penalty_prefactor = input.ip_flux_penalty;
  ip_flux_params.sipg_penalty_fcn = sipg_flux_vector_calc_penalty_maxp2_over_minh;

  central_flux_params_t central_flux_params;
  central_flux_params.central_flux_penalty_prefactor = input.ip_flux_penalty;
  
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  /* create space for storing the ghost data */
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);

  p4est_partition(p4est, 0, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FACE, NULL);
  /* d4est_mesh_geometry_storage_t* geometric_factors = geometric_factors_init(p4est); */

  /* d4est_grid_fcn_t boundary_flux_fcn = zero_fcn; */
  twopunctures_params_t tp_params;
  init_twopunctures_data(&tp_params, input.deg_offset_for_nonlinear_quad);
  /* init_S_puncture_data(p4est, &tp_params, input.deg_offset_for_nonlinear_quad); */

  twopunctures_cactus_params_t tp_cactus_params;
  init_cactus_puncture_data(&tp_cactus_params, input.deg_offset_for_nonlinear_quad);

  
  d4est_elliptic_data_t prob_vecs;
  prob_vecs.rhs = rhs;
  prob_vecs.Au = Au;
  prob_vecs.u = u;
  prob_vecs.local_nodes = local_nodes;

  prob_vecs.flux_fcn_data = curved_gauss_primal_sipg_flux_dirichlet_fetch_fcns
                                          (zero_fcn,&ip_flux_params);

  if(input.use_cactus){
    prob_vecs.user = &tp_cactus_params;
  }
  else {
    prob_vecs.user = &tp_params;
  }

  d4est_elliptic_eqns_t prob_fcns;


  /* if(input.use_cactus){ */
  /*   prob_fcns.build_residual = twopunctures_cactus_build_residual; */
  /*   prob_fcns.apply_lhs = twopunctures_cactus_apply_jac; */
  /* } */
  /* else { */
  prob_vecs.user = &tp_params;
  prob_fcns.build_residual = twopunctures_build_residual;
  prob_fcns.apply_lhs = twopunctures_apply_jac_mg;
  /* } */
  
  d4est_mesh_geometry_storage_t* geometric_factors = geometric_factors_init(p4est);

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
      ghost_data = P4EST_ALLOC(d4est_element_data_t, ghost->ghosts.elem_count);




      d4est_element_data_init_new(p4est,
                                   geometric_factors,
                                   d4est_ops,
                                   d4est_geom,
                                   problem_set_degrees,
                                   (void*)&input,
                                   1,
                                  1);





    
      local_nodes = d4est_element_data_get_local_nodes(p4est);

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
      /*    ip_flux_params.sipg_penalty_prefactor, */
      /*    ghost, */
      /*    ghost_data, */
      /*    d4est_ops, */
      /*    d4est_geom */
      /*   ); */

      
      
      
      /* estimator_stats_t stats; */
      /* estimator_stats_compute(p4est, &stats,1); */

      /* /\* if(world_rank == 0) *\/ */
      /* estimator_stats_print(&stats); */
      


  }


  /* d4est_linalg_fill_vec(prob_vecs.u, 0., local_nodes); */

  
  /* d4est_geom->dxdr_method = INTERP_X_ON_LOBATTO;     */
  /* d4est_element_data_init(p4est, geometric_factors, d4est_ops, d4est_geom, degree, input.gauss_quad_deg); */
  d4est_element_data_init_new(p4est,
                               geometric_factors,
                               d4est_ops,
                               d4est_geom,
                               problem_set_degrees,
                               (void*)&input,1,1);


    
    local_nodes = d4est_element_data_get_local_nodes(p4est);

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

    /* d4est_linalg_fill_vec(prob_vecs.u, 0.001, prob_vecs.local_nodes); */
    d4est_linalg_fill_vec(prob_vecs.u, 0, prob_vecs.local_nodes);
    

    
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
       ip_flux_params.sipg_penalty_prefactor,
       ghost,
       ghost_data,
       d4est_ops,
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

    d4est_element_data_print_number_of_elements_per_tree(p4est);
    
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
    

    double R0 = ((d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R0;
    /* printf("R0 = %.f\n", R0); */

   /*  d4est_geometry_t* d4est_geom_vtk = d4est_geometry */
  /*   p4est_connectivity_t* conn_vtk = p8est_connectivity_new_sphere(); */
  /*   p4est_geometry_t* geom_vtk = p8est_geometry_new_sphere(conn_vtk, R0*3, R0*2, R0); */
  /*   d4est_vtk_context_t* vtk_ctx = d4est_vtk_dg_context_new(p4est, d4est_ops, "puncture-sphere"); */
  /*   d4est_vtk_context_set_geom(vtk_ctx, geom_vtk); */
  /*   d4est_vtk_context_set_scale(vtk_ctx, .99); */
  /*   d4est_vtk_context_set_deg_array(vtk_ctx, deg_array); */
  /*   vtk_ctx = d4est_vtk_write_dg_header(vtk_ctx, d4est_ops); */
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
       d4est_ops,
       u,
       level,
       1,
       R0,
       2*R0,
       3*R0,
       0,
       1,
       input_file
      );
  
    
    hp_amr(p4est,
           d4est_ops,
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
    ghost_data = P4EST_ALLOC(d4est_element_data_t, ghost->ghosts.elem_count);
    
    d4est_element_data_init_new(p4est,
                                 geometric_factors,
                                 d4est_ops,
                                 d4est_geom,
                                 problem_set_degrees,
                                 (void*)&input,1,1);
    
    local_nodes = d4est_element_data_get_local_nodes(p4est);

    u_prev = P4EST_REALLOC(u_prev, double, local_nodes);
    Au = P4EST_REALLOC(Au, double, local_nodes);
    prob_vecs.Au = Au;
    prob_vecs.u = u;
    prob_vecs.u0 = u;
    prob_vecs.local_nodes = local_nodes;

    d4est_linalg_copy_1st_to_2nd(u, u_prev, local_nodes);

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
    /*    d4est_ops, */
    /*    d4est_geom, */
    /*    &prob_fcns, */
    /*    &prob_vecs */
    /*   ); */


    /* double* Au_spec = P4EST_ALLOC(double, local_nodes); */
    /* double* Au_me = P4EST_ALLOC(double, local_nodes); */
    /* double* Au_cactus = P4EST_ALLOC(double, local_nodes); */
    /* double* u_test = P4EST_ALLOC_ZERO(double, local_nodes); */
    
    /* d4est_elliptic_data_t prob_vecs_spec; */
    /* d4est_elliptic_data_t prob_vecs_me; */
    /* d4est_elliptic_data_t prob_vecs_cactus; */
    
    /* d4est_elliptic_data_copy_ptrs(&prob_vecs, &prob_vecs_spec); */
    /* d4est_elliptic_data_copy_ptrs(&prob_vecs, &prob_vecs_me); */
    /* d4est_elliptic_data_copy_ptrs(&prob_vecs, &prob_vecs_cactus); */

    /* prob_vecs_spec.Au = Au_spec; */
    /* prob_vecs_spec.u = u_test; */
    /* prob_vecs_spec.user = &tp_params; */

    /* prob_vecs_me.Au = Au_me; */
    /* prob_vecs_me.u = u_test; */
    /* prob_vecs_me.user = &tp_params; */
    
    /* prob_vecs_cactus.Au = Au_cactus; */
    /* prob_vecs_cactus.u = u_test; */
    /* prob_vecs_cactus.user = &tp_cactus_params; */
    
    /* twopunctures_spec_build_residual */
    /*   ( */
    /*    p4est, */
    /*    ghost, */
    /*    ghost_data, */
    /*    &prob_vecs_spec, */
    /*    d4est_ops, */
    /*    d4est_geom */
    /*   ); */

    /* twopunctures_build_residual */
    /*   ( */
    /*    p4est, */
    /*    ghost, */
    /*    ghost_data, */
    /*    &prob_vecs_me, */
    /*    d4est_ops, */
    /*    d4est_geom */
    /*   ); */
    
    /* prob_vecs_cactus.user = &tp_cactus_params; */
    
    /* twopunctures_cactus_build_residual */
    /*   ( */
    /*    p4est, */
    /*    ghost, */
    /*    ghost_data, */
    /*    &prob_vecs_cactus, */
    /*    d4est_ops, */
    /*    d4est_geom */
    /*   ); */

    
    /* DEBUG_PRINT_ARR_DBL_SUM(Au_spec, local_nodes); */
    /* DEBUG_PRINT_ARR_DBL_SUM(Au_me, local_nodes); */
    /* DEBUG_PRINT_ARR_DBL_SUM(Au_cactus, local_nodes); */
    /* DEBUG_PRINT_ARR_DBL_SUM(u_test, local_nodes); */
    

    /* P4EST_FREE(Au_spec); */
    /* P4EST_FREE(Au_me); */
    /* P4EST_FREE(Au_cactus); */
    /* P4EST_FREE(u_test); */
    
    /* newton_petsc_solve */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs, */
    /*    (void*)&prob_fcns, */
    /*    &ghost, */
    /*    (void**)&ghost_data, */
    /*    d4est_ops, */
    /*    d4est_geom, */
    /*    input_file, */
    /*    NULL */
    /*   ); */


    int min_level, max_level;

    multigrid_get_level_range(p4est, &min_level, &max_level);
    printf("[min_level, max_level] = [%d,%d]\n", min_level, max_level);

    /* need to do a reduce on min,max_level before supporting multiple proc */
    /* D4EST_ASSERT(proc_size == 1); */
    int num_of_levels = max_level + 1;
     
    

    multigrid_logger_t* logger = multigrid_logger_residual_init
                                 (
                                 );

    multigrid_element_data_updater_t* updater
      = multigrid_element_data_updater_curved_init
      (
       num_of_levels,
       &ghost,
       &ghost_data,
       geometric_factors,
       d4est_geom,
       set_deg_quad,
       &input
      );

    multigrid_user_callbacks_t* matrix_op_callbacks = multigrid_matrix_operator_init
                                                      (
                                                       p4est,
                                                       num_of_levels,
                                                       d4est_ops,
                                                       d4est_element_data_get_local_matrix_nodes,
                                                       &tp_params
                                                      );
    
    multigrid_data_t* mg_data = multigrid_data_init(p4est,
                                                    d4est_ops,
                                                    num_of_levels,
                                                    logger,
                                                    /* NULL, //matrix_op_callbacks, */
                                                    matrix_op_callbacks,
                                                    updater,
                                                    input_file
                                                   );

    
    prob_vecs.user = matrix_op_callbacks->user;
    prob_fcns.build_residual = twopunctures_build_residual_mg;
    prob_fcns.apply_lhs = twopunctures_apply_jac_mg;
    
    krylov_pc_t* pc = krylov_pc_multigrid_create(mg_data, twopunctures_compute_jac_matrix_operator_for_pc);
    /* krylov_pc_multigrid_setup(pc); */
      
    
    /* krylov_petsc_params_t krylov_petsc_params; */
    /* krylov_petsc_input(p4est, input_file, "krylov_petsc", "[KRYLOV_PETSC]", &krylov_petsc_params);       */

    /* newton_petsc_params_t newton_petsc_params; */
    /* newton_petsc_input(p4est, input_file, "[NEWTON_PETSC]", &newton_petsc_params);        */

    
    /* newton_petsc_solve */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs, */
    /*    (void*)&prob_fcns, */
    /*    &ghost, */
    /*    (void**)&ghost_data, */
    /*    d4est_ops, */
    /*    d4est_geom, */
    /*    &krylov_petsc_params, */
    /*    &newton_petsc_params, */
    /*    pc */
    /*   ); */

    /* newton_petsc_solve */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs, */
    /*    (void*)&prob_fcns, */
    /*    &ghost, */
    /*    (void**)&ghost_data, */
    /*    d4est_ops, */
    /*    d4est_geom, */
    /*    &krylov_petsc_params, */
    /*    &newton_petsc_params, */
    /*    NULL */
    /*   ); */



    /* petsc_ctx_t petsc_ctx; */
    /* petsc_ctx.p4est = p4est; */
    /* petsc_ctx.vecs = &prob_vecs; */
    /* petsc_ctx.fcns = &prob_fcns; */
    /* petsc_ctx.ghost = &ghost; */
    /* petsc_ctx.ghost_data = (void**)&ghost_data; */
    /* petsc_ctx.d4est_ops = d4est_ops; */
    /* petsc_ctx.d4est_geom = d4est_geom; */
    /* pc->pc_ctx = &petsc_ctx; */
    /* pc->pc_setup(pc); */
    
    /* d4est_linalg_fill_vec(prob_vecs.u, 1., prob_vecs.local_nodes); */
    /* twopunctures_apply_jac_mg(p4est, ghost, ghost_data, &prob_vecs, d4est_ops, d4est_geom); */
    /* double sum1 = d4est_linalg_vec_sum(prob_vecs.Au,prob_vecs.local_nodes); */

    /* void* tmp_ptr = prob_vecs.user; */
    /* prob_vecs.user = ((multigrid_matrix_op_t*)tmp_ptr)->user; */
    /* twopunctures_apply_jac(p4est, ghost, ghost_data, &prob_vecs, d4est_ops, d4est_geom); */
    
    /* double sum2 = d4est_linalg_vec_sum(prob_vecs.Au,prob_vecs.local_nodes); */
    /* prob_vecs.user = tmp_ptr; */

    /* printf("sum1, sum2 = %f,%f\n", sum1, sum2); */
    
    newton_d4est_solve
      (
       p4est,
       &prob_vecs,
       (void*)&prob_fcns,
       &ghost,
       (void**)&ghost_data,
       d4est_ops,
       d4est_geom,
       input_file,
       pc
      );

    krylov_pc_multigrid_destroy(pc);
    multigrid_logger_residual_destroy(logger);
    multigrid_element_data_updater_curved_destroy(updater, num_of_levels);
    multigrid_matrix_operator_destroy(matrix_op_callbacks);
    multigrid_data_destroy(mg_data);

    
    
    /* matrix_spd_tester_parallel */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs,  */
    /*    &prob_fcns, */
    /*    ghost, */
    /*    ghost_data, */
    /*    d4est_ops, */
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
    /*    d4est_ops, */
    /*    d4est_geom, */
    /*    1, */
    /*    20, */
    /*    .000000001 */
    /*   );  */

    d4est_linalg_vec_axpy(-1., prob_vecs.u, u_prev, local_nodes);

    double local_l2_norm_sqr = d4est_element_data_compute_l2_norm_sqr
                                (
                                 p4est,
                                 u_prev,
                                 local_nodes,
                                 d4est_ops,
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
  P4EST_FREE(u_prev);
  /* P4EST_FREE(u_analytic); */
  P4EST_FREE(u);
}


 
