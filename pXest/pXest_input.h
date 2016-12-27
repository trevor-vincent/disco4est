#ifndef PXEST_INPUT_H
#define PXEST_INPUT_H 

#include "../InputParser/ini.h"
#include "./pXest.h"
#include "../Utilities/util.h"

typedef struct {

  int min_quadrants;
  int min_level;
  int fill_uniform;

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
  if (util_match(section,"p4est",name,"min_quadrants")) {
    pconfig->min_quadrants = atoi(value);
  } else if (util_match(section,"p4est",name, "min_level")) {
    pconfig->min_level = atoi(value);
  } else if (util_match(section,"p4est",name,"fill_uniform")) {
    pconfig->fill_uniform = atoi(value);
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
  pXest_input.min_quadrants = -1;
  pXest_input.min_level = 0;
  pXest_input.fill_uniform = 1;
  
  if (ini_parse("options.input", pXest_input_handler, &pXest_input) < 0) {
    mpi_abort("[D4EST_ERROR]: Can't load pXest input file");
  }

  if(
     proc_size != 1
     &&
     util_int_pow_int((P4EST_CHILDREN), pXest_input.min_level) < proc_size
     &&
     pXest_input.min_quadrants < proc_size
  ){
    if (proc_rank == 0){
      printf("[D4EST_WARNING]: Starting p4est with elements < processes\n");
    }
  }

  return pXest_input;
}

#endif
