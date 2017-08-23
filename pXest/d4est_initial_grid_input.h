#ifndef D4EST_INITIAL_GRID_INPUT_H
#define D4EST_INITIAL_GRID_INPUT_H 

#include <ini.h>
#include <pXest.h>
#include <d4est_util.h>
#include <d4est_element_data.h>

typedef struct {

  int min_quadrants;
  int min_level;
  int fill_uniform;
  int print_elements_per_proc;
  int deg;
  int deg_quad;
  
} d4est_initial_grid_t;

static
int d4est_initial_grid_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_initial_grid_t* pconfig = (d4est_initial_grid_t*) user;
  if (d4est_util_match_couple(section,"initial_grid",name,"min_quadrants")) {
    D4EST_ASSERT(pconfig->min_quadrants == -2);
    pconfig->min_quadrants = atoi(value);
  }
  else if (d4est_util_match_couple(section,"initial_grid",name, "min_level")) {
    D4EST_ASSERT(pconfig->min_level == -2);
    pconfig->min_level = atoi(value);
  }
  else if (d4est_util_match_couple(section,"initial_grid",name,"fill_uniform")) {
    D4EST_ASSERT(pconfig->fill_uniform == -2);
    pconfig->fill_uniform = atoi(value);
  }
  else if (d4est_util_match_couple(section,"initial_grid",name,"print_elements_per_proc")) {
    D4EST_ASSERT(pconfig->print_elements_per_proc == 0);
    pconfig->print_elements_per_proc = atoi(value);
  }
  else if (d4est_util_match_couple(section,"initial_grid",name,"deg")) {
    D4EST_ASSERT(pconfig->deg == -1);
    pconfig->deg = atoi(value);
  }
  else if (d4est_util_match_couple(section,"initial_grid",name,"deg_quad")) {
    D4EST_ASSERT(pconfig->deg_quad == -1);
    pconfig->deg_quad = atoi(value);
  }   
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


void
d4est_initial_grid_degree_init
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  d4est_initial_grid_t* input = user_ctx;
  elem_data->deg = input->deg;
  elem_data->deg_vol_quad = input->deg_quad;
}


d4est_initial_grid_t
d4est_initial_grid_parse(const char* input_file)
{

  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;  
  /* SC_CHECK_MPI(mpiret); */

  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  
  d4est_initial_grid_t d4est_initial_grid;
  d4est_initial_grid.min_quadrants = -2;
  d4est_initial_grid.min_level = -2;
  d4est_initial_grid.fill_uniform = -2;
  d4est_initial_grid.print_elements_per_proc = 0;
  d4est_initial_grid.deg = -1;
  d4est_initial_grid.deg_quad = -1;
  
  if (ini_parse(input_file, d4est_initial_grid_handler, &d4est_initial_grid) < 0) {
    printf("[D4EST_ERROR]: pXest input_file = %s\n", input_file);
    D4EST_ABORT("[D4EST_ERROR]: Can't load pXest input file");
  }

  D4EST_CHECK_INPUT("initial_grid", d4est_initial_grid.min_quadrants, -2);
  D4EST_CHECK_INPUT("initial_grid", d4est_initial_grid.min_level, -2);
  D4EST_CHECK_INPUT("initial_grid", d4est_initial_grid.fill_uniform, -2);
  D4EST_CHECK_INPUT("initial_grid", d4est_initial_grid.deg, -1);
  D4EST_CHECK_INPUT("initial_grid", d4est_initial_grid.deg_quad, -1);
  
  if(
     proc_size != 1
     &&
     d4est_util_int_pow_int((P4EST_CHILDREN), d4est_initial_grid.min_level) < proc_size
     &&
     d4est_initial_grid.min_quadrants < proc_size
  ){
    if (proc_rank == 0){
      printf("[D4EST_ERROR]: proc_size = %d\n", proc_size);
      printf("[D4EST_ERROR]: min_level = %d\n", d4est_initial_grid.min_level);
      D4EST_ABORT("[D4EST_ERROR]: Starting p4est with elements < processes\n");
    }
  }

  return d4est_initial_grid;
}

#endif


