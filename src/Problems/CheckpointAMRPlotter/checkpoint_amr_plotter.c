#include <sc_reduce.h>
#include <pXest.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <problem.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_estimator_bi.h>
#include <d4est_solver_cg.h>
#include <d4est_amr.h>
#include <d4est_amr_smooth_pred.h>
#include <d4est_geometry.h>
#include <d4est_geometry_brick.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_vtk.h>
#include <d4est_norms.h>
#include <d4est_mesh.h>
#include <ini.h>
#include <d4est_element_data.h>
#include <d4est_estimator_stats.h>
#include <d4est_laplacian_with_opt.h>
#include <d4est_laplacian_with_opt_flux_sipg.h>
#include <d4est_solver_newton.h>
#include <d4est_solver_multigrid.h>
#include <d4est_krylov_pc_multigrid.h>
#include <d4est_solver_multigrid_logger_residual.h>
#include <d4est_solver_multigrid_element_data_updater.h>
#include <d4est_solver_multigrid_matrix_operator.h>
#include <d4est_solver_krylov_petsc.h>
#include <d4est_solver_newton_petsc.h>
#include <d4est_util.h>
#include <d4est_h5.h>
#include <d4est_checkpoint.h>
#include <time.h>
#include "two_punctures_cactus_fcns_with_opt.h"


void
problem_init
(
 p4est_t* p4est,
 d4est_ghost_t** d4est_ghost,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_mesh_initial_extents_t* initial_extents,
 const char* input_file,
 sc_MPI_Comm mpicomm
)
{ 
  int initial_nodes = initial_extents->initial_nodes;
  int min_level = initial_extents->initial_nodes;
  int number_of_regions = initial_extents->number_of_regions;
  int* initial_deg_per_region = initial_extents->deg;

  int number_of_intervals_dx = params->number_of_intervals_dx;
  int x0 = params->number_of_intervals_x0;
  int x1 = params->number_of_intervals_x1;
  
  int* num_refinements = P4EST_ALLOC_ZERO(int, 2*(number_of_regions + number_of_intervals_dx));
  int* global_num_refinements = P4EST_ALLOC_ZERO(int, 2*(number_of_regions + number_of_intervals_dx));
    
  for (int i = 0; i < p4est->local_num_quadrants; i++){
    d4est_element_data_t* ed = d4est_factors->element_data[i];
    double* x_on_ed = d4est_factors->xyz[ed->nodal_stride];
    double x_max = x_on_ed[0];
    for (int j = 0; j < volume_nodes_ed; j++){
      x_max = (x_on_ed[j] > x_max) ? x_on_ed[j] : x_max;
    }
    int interval = floor((xmax - x0)/dx + .5);
    int region = ed->region;
    int h_refinements = ed->level - min_level;
    int p_refinements = ed->deg - initial_deg_per_region[region];
    num_refinements[region] += h_refinements;
    num_h_refinements_per_dx[number_of_regions + interval] += h_refinements;
    num_p_refinements_per_region[number_of_regions + number_of_intervals + region] += p_refinements;
    num_p_refinements_per_dx[number_of_regions + number_of_intervals + number_of_regions + interval] += p_refinements;
  }

  sc_reduce(
    &num_refinements,
    &global_num_refinements,
    2*(number_of_regions + number_of_intervals_dx),
    sc_MPI_INT,
    sc_MPI_SUM,
    0,
    sc_MPI_COMM_WORLD
  );

  zlog_category_t* refinements_per_region = zlog_get_category("d4est_amr_diagnostic_per_region");
  for (int i = 0; i < number_of_regions; i++){
    zlog_info("%d %d %d ", i, global_num_refinements[i], global_num_refinements[number_of_regions + number_of_intervals + i]);
  }

  zlog_category_t* refinements_per_interval = zlog_get_category("d4est_amr_diagnostic_per_interval");
  for (int i = 0; i < number_of_intervals; i++){
    zlog_info("%.15f %d %d ",x0 + i*dx, global_num_refinements[number_of_regions + i]
              , global_num_refinements[2*number_of_regions + number_of_intervals + i]);
  }  
  /* sc_reduce to zero node */
  /* print to file on zero nodewith zlog */
  /* add to two punctures file and lorentzian file */
}
