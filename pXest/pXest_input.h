#ifndef PXEST_INPUT_H
#define PXEST_INPUT_H 

#include <ini.h>
#include <pXest.h>
#include <d4est_util.h>

typedef struct {

  int min_quadrants;
  int min_level;
  int fill_uniform;
  int print_elements_per_proc;
  
} pXest_input_t;

static
int pXest_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  pXest_input_t* pconfig = (pXest_input_t*)user;
  if (d4est_util_match_couple(section,"initial_grid",name,"min_quadrants")) {
    pconfig->min_quadrants = atoi(value);
  } else if (d4est_util_match_couple(section,"initial_grid",name, "min_level")) {
    pconfig->min_level = atoi(value);
  } else if (d4est_util_match_couple(section,"initial_grid",name,"fill_uniform")) {
    pconfig->fill_uniform = atoi(value);
  } else if (d4est_util_match_couple(section,"initial_grid",name,"print_elements_per_proc")) {
    pconfig->print_elements_per_proc = atoi(value);
  } else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

pXest_input_t
pXest_input_parse(const char* input_file){


  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;  
  /* SC_CHECK_MPI(mpiret); */

  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  
  pXest_input_t pXest_input;
  pXest_input.min_quadrants = -2;
  pXest_input.min_level = -2;
  pXest_input.fill_uniform = -2;
  pXest_input.print_elements_per_proc = 0;
  
  if (ini_parse(input_file, pXest_input_handler, &pXest_input) < 0) {
    printf("[D4EST_ERROR]: pXest input_file = %s\n", input_file);
    D4EST_ABORT("[D4EST_ERROR]: Can't load pXest input file");
  }

  D4EST_CHECK_INPUT("initial_grid", pXest_input.min_quadrants, -2);
  D4EST_CHECK_INPUT("initial_grid", pXest_input.min_level, -2);
  D4EST_CHECK_INPUT("initial_grid", pXest_input.fill_uniform, -2);
  
  if(
     proc_size != 1
     &&
     d4est_util_int_pow_int((P4EST_CHILDREN), pXest_input.min_level) < proc_size
     &&
     pXest_input.min_quadrants < proc_size
  ){
    if (proc_rank == 0){
      printf("[D4EST_ERROR]: proc_size = %d\n", proc_size);
      printf("[D4EST_ERROR]: min_level = %d\n", pXest_input.min_level);
      D4EST_ABORT("[D4EST_ERROR]: Starting p4est with elements < processes\n");
    }
  }

  return pXest_input;
}

#endif
