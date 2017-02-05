#include "d4est_geometry.h"


typedef struct {

  int print_info;
  const char* name;
  int count;
  
} d4est_geometry_input_t;


static
int problem_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_geometry_input_t* pconfig = (d4est_geometry_input_t*)user;
  if (util_match(section,"geometry",name,"print_info")) {
    mpi_assert(pconfig->num_unifrefs == -1);
    pconfig->print_info = atoi(value);
    pconfig->count += 1;
  }
  else if (util_match(section,"geometry",name,"name")) {
    mpi_assert(pconfig->namer == NULL);
    pconfig->name = value;
    pconfig->count += 1;
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

d4est_geometry_input_t
problem_input
(
 const char* input_file
)
{
  int num_of_options = 2;
  
  d4est_geometry_input_t input;
  input.count = 0;
  input.print_info = -1;
  input.name = NULL;
  
  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  mpi_assert(input.count == num_of_options);
  return input;
}


d4est_geometry_t*
d4est_geometry_new(const char* input_file){

  d4est_geometry_input_t input = problem_input(input_file);
  
  if (util_match(section,"geometry")) {

  
}
