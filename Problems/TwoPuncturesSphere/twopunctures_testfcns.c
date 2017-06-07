#include <sc_reduce.h>
#include <pXest.h>
#include <util.h>
#include <d4est_linalg.h>
#include <d4est_element_data.h>
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
#include <d4est_geometry_cubed_sphere.h>
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
#include "./twopuncturesfcns_spec.h"
#include "time.h"
#include "util.h"

/* soon to be in the input files */
static const double pi = 3.1415926535897932384626433832795;

static
void
twopunctures_test_build_residual
(
 p4est_t* p4est,
 twopunctures_params_t* tp_params,
 twopunctures_cactus_params_t* tp_cactus_params,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
 
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data;
        int volume_nodes_Gauss = d4est_operators_get_nodes((P4EST_DIM), ed->deg_quad);
        printf("\n Element %d\n", ed->id);
        for (int i = 0; i < volume_nodes_Gauss; i++){
          double x = ed->xyz_quad[0][i];
          double y = ed->xyz_quad[1][i];
          double z = ed->xyz_quad[2][i];
 

          /* if (ed->id == 103){ */
          double Aij_cactus [3][3];
          double Aij_me [3][3];
          double Aij_spec [3][3];
          BY_Aijofxyz(x,y,z, Aij_cactus, tp_cactus_params);
          twopunctures_spec_compute_multiBH_tildeAij(x,y,z,tp_params->P_bh, tp_params->S_bh, tp_params->C_bh, tp_params->num_punctures, Aij_spec);

          
          for (int k = 0; k < 3; k++)
            for (int j = 0; j < 3; j++){
              Aij_me[k][j] = compute_confAij(k,j,x,y,z,tp_params);
              double error = fabs(Aij_me[k][j] - Aij_cactus[k][j]) + fabs(Aij_me[k][j] - Aij_spec[k][j]);
              /* if (error > 1e-8){ */
                /* printf("Holy fucks batman error is big on this bitch\n"); */
              /* } */
              /* printf("x,y,z,i,j, Aij_me, Aij_cactus, Aij_spec = %.7f %.7f %.7f %d %d %.25f %.25f %.25f %.25f\n", */
                     /* x,y,z,k,j, Aij_me[k][j], Aij_cactus[k][j], Aij_spec[k][j], error); */
            }


          double tp_me_aijsqr = compute_confAij_sqr(x,y,z, tp_params);
          double tp_cactus_aijsqr = BY_KKofxyz(x,y,z, tp_cactus_params);
          double tp_spec_aijsqr = twopunctures_spec_compute_tildeAijsqr(x,y,z, tp_params->P_bh, tp_params->S_bh, tp_params->C_bh, tp_params->num_punctures);
          double error_aijsqr = fabs(tp_me_aijsqr - tp_cactus_aijsqr) + fabs(tp_me_aijsqr - tp_spec_aijsqr);
          /* printf("x,y,z, me_aijsqr, cactus_aijsqr, spec_aijsqr, error = %.6f %.6f %.6f %.25f %.25f %.25f %.25f\n", x,y,z, tp_me_aijsqr, tp_cactus_aijsqr, tp_spec_aijsqr, error_aijsqr); */

         double u = 0.;
          double tp_me_neg_1o8_K2_psi_neg7 =
            twopunctures_neg_1o8_K2_psi_neg7(x,y,z,u,tp_params);
          double tp_spec_neg_1o8_K2_psi_neg7 =
            twopunctures_spec_neg_1o8_K2_psi_neg7(x,y,z,u,tp_params);
          double tp_cactus_neg_1o8_K2_psi_neg7 =
            twopunctures_cactus_neg_1o8_K2_psi_neg7(x,y,z,u,tp_cactus_params);
          /* printf("i,x,y,z, me, spec, cactus = %d %.7f %.7f %.7f %.25f %.25f %.25f\n",i, */
                 /* x,y,z,tp_me_neg_1o8_K2_psi_neg7, tp_spec_neg_1o8_K2_psi_neg7, */
                 /* tp_cactus_neg_1o8_K2_psi_neg7); */

          double tp_me_plus_7o8_K2_psi_neg8 =
            twopunctures_plus_7o8_K2_psi_neg8(x,y,z,u,tp_params);
          double tp_spec_plus_7o8_K2_psi_neg8 =
            twopunctures_spec_plus_7o8_K2_psi_neg8(x,y,z,u,tp_params);
          double tp_cactus_plus_7o8_K2_psi_neg8 =
            twopunctures_cactus_plus_7o8_K2_psi_neg8(x,y,z,u,tp_cactus_params);


          double batman_error = fabs(tp_me_plus_7o8_K2_psi_neg8 - tp_spec_plus_7o8_K2_psi_neg8);
          batman_error += fabs(tp_me_plus_7o8_K2_psi_neg8 - tp_cactus_plus_7o8_K2_psi_neg8);
          batman_error += fabs(tp_me_neg_1o8_K2_psi_neg7  - tp_cactus_neg_1o8_K2_psi_neg7);
          batman_error += fabs(tp_me_neg_1o8_K2_psi_neg7  - tp_spec_neg_1o8_K2_psi_neg7);
          if (batman_error > 1e-16){
            printf("Holy fucks batman error is big on this bitch\n");
          }
          
          /* printf("i,x,y,z, me, spec, cactus = %d %.7f %.7f %.7f %.25f %.25f %.25f\n",i, */
          /*        x,y,z,tp_me_plus_7o8_K2_psi_neg8, tp_spec_plus_7o8_K2_psi_neg8, */
          /*        tp_cactus_plus_7o8_K2_psi_neg8);       */    
          
          /* } */
        }
      }
    }

}


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
  else if (util_match_couple(section,"problem",name,"deg_quad_R0")) {
    mpi_assert(pconfig->deg_quad_R0 == -1);
    pconfig->deg_quad_R0 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_R1")) {
    mpi_assert(pconfig->deg_R1 == -1);
    pconfig->deg_R1 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_quad_R1")) {
    mpi_assert(pconfig->deg_quad_R1 == -1);
    pconfig->deg_quad_R1 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_R2")) {
    mpi_assert(pconfig->deg_R2 == -1);
    pconfig->deg_R2 = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"problem",name,"deg_quad_R2")) {
    mpi_assert(pconfig->deg_quad_R2 == -1);
    pconfig->deg_quad_R2 = atoi(value);
    pconfig->count += 1;
  }  
  else if (util_match_couple(section,"problem",name,"deg_offset_for_nonlinear_quad")) {
    mpi_assert(pconfig->deg_offset_for_nonlinear_quad == -1);
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
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  problem_input_t* input = user_ctx;
  /* outer shell */
  if (elem_data->tree < 6){
    elem_data->deg = input->deg_R2;
    elem_data->deg_quad = input->deg_quad_R2;
  }
  /* inner shell */
  else if(elem_data->tree < 12){
    elem_data->deg = input->deg_R1;
    elem_data->deg_quad = input->deg_quad_R1;
  }
  /* center cube */
  else {
    elem_data->deg = input->deg_R0;
    elem_data->deg_quad = input->deg_quad_R0;
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
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);

  p4est_partition(p4est, 0, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FACE, NULL);
  /* d4est_mesh_geometry_storage_t* geometric_factors = geometric_factors_init(p4est); */


  /* grid_fcn_t boundary_flux_fcn = zero_fcn; */
  twopunctures_params_t tp_params;
  init_twopunctures_data(&tp_params, input.deg_offset_for_nonlinear_quad);
  /* init_S_puncture_data(p4est, &tp_params, input.deg_offset_for_nonlinear_quad); */

  twopunctures_cactus_params_t tp_cactus_params;
  init_cactus_puncture_data(&tp_cactus_params, input.deg_offset_for_nonlinear_quad);

  
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
    prob_fcns.build_residual = twopunctures_cactus_build_residual;
    prob_fcns.apply_lhs = twopunctures_cactus_apply_jac;
  }
  else {
    prob_fcns.build_residual = twopunctures_build_residual;
    prob_fcns.apply_lhs = twopunctures_apply_jac;
  }
  
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
                                   (void*)&input);





    
      local_nodes = d4est_element_data_get_local_nodes(p4est);

      Au = P4EST_REALLOC(Au, double, local_nodes);
      u = P4EST_REALLOC(u, double, local_nodes);
      rhs = P4EST_REALLOC(rhs, double, local_nodes);
      prob_vecs.Au = Au;
      prob_vecs.u = u;
      prob_vecs.rhs = rhs;
      prob_vecs.local_nodes = local_nodes;
     


  }


  /* d4est_linalg_fill_vec(prob_vecs.u, 0., local_nodes); */

  
  d4est_geom->dxdr_method = INTERP_X_ON_LOBATTO;    
  d4est_element_data_init_new(p4est,
                               geometric_factors,
                               d4est_ops,
                               d4est_geom,
                               problem_set_degrees,
                               (void*)&input);


    
    local_nodes = d4est_element_data_get_local_nodes(p4est);

    Au = P4EST_REALLOC(Au, double, local_nodes);
    u = P4EST_REALLOC(u, double, local_nodes);
    rhs = P4EST_REALLOC(rhs, double, local_nodes);

    prob_vecs.Au = Au;
    prob_vecs.u = u;
    prob_vecs.rhs = rhs;
    prob_vecs.local_nodes = local_nodes;

twopunctures_test_build_residual
(
 p4est,
 &tp_params,
 &tp_cactus_params,
 d4est_ops,
 d4est_geom
);

    /* double* Au_spec = P4EST_ALLOC(double, local_nodes); */
    /* double* Au_me = P4EST_ALLOC(double, local_nodes); */
    /* double* Au_cactus = P4EST_ALLOC(double, local_nodes); */
    /* double* u_test = P4EST_ALLOC_ZERO(double, local_nodes); */
    
    /* problem_data_t prob_vecs_spec; */
    /* problem_data_t prob_vecs_me; */
    /* problem_data_t prob_vecs_cactus; */
    
    /* problem_data_copy_ptrs(&prob_vecs, &prob_vecs_spec); */
    /* problem_data_copy_ptrs(&prob_vecs, &prob_vecs_me); */
    /* problem_data_copy_ptrs(&prob_vecs, &prob_vecs_cactus); */

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
