#define _GNU_SOURCE
#include <pXest.h>
#include <d4est_solver_multigrid.h>
#include <d4est_solver_multigrid_mesh_analyzer.h>
#include <d4est_vtk.h>
#include <sc_reduce.h>
#include <time.h>
#include <stdio.h>


void
d4est_solver_multigrid_mesh_analyzer_destroy
(
 d4est_solver_multigrid_mesh_analyzer_t* mesh_analyzer
)
{
  /* free(mesh_analyzer->input_file); */
  P4EST_FREE(mesh_analyzer);
}

static
void d4est_solver_multigrid_mesh_analyzer_save_vtk
(
 p4est_t* p4est
)
{
  d4est_solver_multigrid_data_t* mg_data = p4est->user_pointer;
  d4est_solver_multigrid_mesh_analyzer_t* mesh_analyzer = mg_data->analyzer;
  d4est_operators_t* d4est_ops = mg_data->d4est_ops;
  /* char* output = "d4est_solver_multigrid_mesh_analyzer"; */
  d4est_vtk_save_aux
    (
     p4est,
     d4est_ops,
     mg_data->input_file,
     "d4est_vtk",
     (const char * []){NULL},
     (double* []){},
     (const char * []){NULL},
     (double* []){},
     NULL,
     NULL,
     "multigrid_mesh_analyzer",
     mesh_analyzer->stride
    );
  /* free(output); */
}

static
void d4est_solver_multigrid_mesh_analyzer_update
(
 p4est_t* p4est,
 int level,
 d4est_elliptic_data_t* vecs
)
{
  d4est_solver_multigrid_data_t* mg_data = p4est->user_pointer;
  d4est_solver_multigrid_mesh_analyzer_t* mesh_analyzer = mg_data->analyzer;
  
  if(mg_data->mg_state == START){
    mesh_analyzer->stride = (P4EST_MAXLEVEL);
    /* printf("P4EST_MAXLEVEL = %d\n", P4EST_MAXLEVEL); */
    mesh_analyzer->levels = 1;
    mesh_analyzer->p4est_checksums[mesh_analyzer->stride]
      = p4est_checksum(p4est);
    mesh_analyzer->deg_checksums[mesh_analyzer->stride]
      = d4est_mesh_get_local_degree_sum(p4est);
    d4est_solver_multigrid_mesh_analyzer_save_vtk(p4est);
    mesh_analyzer->stride -= 1;
  }
  else if (mg_data->mg_state == PRE_V){
  }
  else if (mg_data->mg_state == DOWNV_PRE_SMOOTH){
  }
  else if (mg_data->mg_state == DOWNV_POST_SMOOTH){
  }
  else if (mg_data->mg_state == DOWNV_PRE_COARSEN){
  }
  else if (mg_data->mg_state == DOWNV_POST_COARSEN){
  }  
  else if (mg_data->mg_state == DOWNV_PRE_BALANCE){
  }
  else if (mg_data->mg_state == DOWNV_POST_BALANCE){
 
  }
  else if (mg_data->mg_state == DOWNV_PRE_RESTRICTION){
    mesh_analyzer->p4est_checksums[mesh_analyzer->stride]
      = p4est_checksum(p4est);
    mesh_analyzer->deg_checksums[mesh_analyzer->stride]
      = d4est_mesh_get_local_degree_sum(p4est);
    d4est_solver_multigrid_mesh_analyzer_save_vtk(p4est);
    mesh_analyzer->stride -= 1;
    mesh_analyzer->levels += 1;    
  }
  else if (mg_data->mg_state == DOWNV_POST_RESTRICTION){
  }    
  else if (mg_data->mg_state == COARSE_PRE_SOLVE){
  }
  else if (mg_data->mg_state == COARSE_POST_SOLVE){
    mesh_analyzer->stride += 1;
  }
  else if (mg_data->mg_state == UPV_PRE_REFINE){
  }
  else if (mg_data->mg_state == UPV_POST_REFINE){
    mesh_analyzer->stride += 1;
    int mesh_checksum = p4est_checksum(p4est);
    int deg_checksum = d4est_mesh_get_local_degree_sum(p4est);
    if (mesh_checksum != mesh_analyzer->p4est_checksums[mesh_analyzer->stride] ||
        deg_checksum != mesh_analyzer->deg_checksums[mesh_analyzer->stride]){
      printf("mesh_analyzer->levels = %d\n", mesh_analyzer->levels);
      printf("mesh_analyzer->stride = %d\n", mesh_analyzer->stride);
      printf("mesh_analyzer->p4est_checksums[mesh_analyzer->stride] = %d\n", (mesh_analyzer->p4est_checksums[mesh_analyzer->stride]));
      printf("mesh_analyzer->deg_checksums[mesh_analyzer->stride] = %d\n", (mesh_analyzer->deg_checksums[mesh_analyzer->stride]));
      printf("p4est_checksum = %d\n",mesh_checksum);
      printf("deg_checksum = %d\n",deg_checksum);
      for (int level = 0; level < mesh_analyzer->levels; level++){
        printf("level %d: p4est_checksum = %d\n", level, mesh_analyzer->p4est_checksums[(P4EST_MAXLEVEL) - level]);
        printf("level %d: deg_checksum = %d\n",level,mesh_analyzer->deg_checksums[(P4EST_MAXLEVEL) - level]);
      }
      D4EST_ABORT("checksum != mesh_analyzer->checksums[mesh_analyzer->stride]");
    }
  }
  else if (mg_data->mg_state == POST_V){

    double biggest [4];
    double biggest_global[4];
    zlog_category_t *c_default = zlog_get_category("d4est_solver_multigrid_mesh_analyzer");

    for (int lev = (mesh_analyzer->levels)-1; lev >= 0; lev--){

      biggest[0] = 1./(double)mg_data->elements_on_level_of_multigrid[lev];
      biggest[1] =(double) mg_data->elements_on_level_of_multigrid[lev];
      biggest[2] = 1./(double)mg_data->nodes_on_level_of_multigrid[lev];
      biggest[3] = (double)mg_data->nodes_on_level_of_multigrid[lev];

      double* biggest_ptr = &biggest[0];
      /* DEBUG_PRINT_ARR_DBL(biggest_ptr, 4); */
      
      sc_reduce
        (
         &biggest[0],
         &biggest_global[0],
         4,
         sc_MPI_DOUBLE,
         sc_MPI_MAX,
         0,
         p4est->mpicomm
        );

      if (p4est->mpirank == 0){
        zlog_debug(c_default,
                   "For each processor, we have Lev %d, min K = %.0f, max K = %.0f, min N = %.0f, max N = %.0f",
                   lev,
                   (1./biggest_global[0]),
                   (biggest_global[1]),
                   (1./biggest_global[2]),
                   (biggest_global[3])
                  );
      }
    }
  }
  else if (mg_data->mg_state == END){
  }  
  else {
    return;
  }
}

d4est_solver_multigrid_mesh_analyzer_t*
d4est_solver_multigrid_mesh_analyzer_init
(
 /* const char* input_file */
)
{
  d4est_solver_multigrid_mesh_analyzer_t* mesh_analyzer = P4EST_ALLOC(d4est_solver_multigrid_mesh_analyzer_t, 1);
  /* mesh_analyzer->input_file = malloc(strlen(input_file) + 1); */
  /* printf("input_file = %s\n", input_file); */
  /* strcpy(mesh_analyzer->input_file, input_file); */
  /* printf("input_file_stored = %s\n", mesh_analyzer->input_file); */
  mesh_analyzer->update = d4est_solver_multigrid_mesh_analyzer_update;
  return mesh_analyzer;
}
