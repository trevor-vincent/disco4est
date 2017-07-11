#include <d4est_geometry.h>
#include <d4est_geometry_cubed_sphere_outer_shell_block.h>
#include <d4est_util.h>
#include <ini.h>


typedef struct {

  double R0;
  double R1;
  int compactify
  
} d4est_geometry_cubed_sphere_outer_shell_block_input_t;

static
int d4est_geometry_cubed_sphere_outer_shell_block_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_geometry_cubed_sphere_outer_shell_block_input_t* pconfig = (d4est_geometry_cubed_sphere_outer_shell_block_input_t*)user;
  if (d4est_util_match_couple(section,"geometry",name,"R0")) {
    D4EST_ASSERT(pconfig->R0 == -1);
    pconfig->R0 = atof(value);
  }
  else if (d4est_util_match_couple(section,"geometry",name,"R1")) {
    D4EST_ASSERT(pconfig->R1 == -1);
    pconfig->R1 = atof(value);
  }
  else if (d4est_util_match_couple(section,"geometry",name,"compactify")) {
    D4EST_ASSERT(pconfig->compactify == -1);
    pconfig->compactify = atoi(value);
    D4EST_ASSERT(pconfig->compactify == 0 || pconfig->compactify == 1);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

d4est_geometry_cubed_sphere_outer_shell_block_input_t*
d4est_geometry_cubed_sphere_outer_shell_block_input
(
 const char* input_file
)
{
  
  d4est_geometry_cubed_sphere_outer_shell_block_input_t* input
    = P4EST_ALLOC(d4est_geometry_cubed_sphere_outer_shell_block_input_t, 1);
  
  input->R0 = -1;
  input->R1 = -1;
  
  if (ini_parse(input_file, d4est_geometry_cubed_sphere_outer_shell_block_input_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("geometry", input->R0, -1);
  D4EST_CHECK_INPUT("geometry", input->R1, -1);
  
  return input;
}




static void
d4est_geometry_cubed_sphere_outer_shell_block_destroy
(
 p4est_geometry_t* geom
)
{
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}


void
d4est_geometry_cubed_sphere_outer_shell_block_new
(
 const char* input_file,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_cubed_sphere_outer_shell_block_input_t input = d4est_geometry_cubed_sphere_outer_shell_block_input(input_file);
  p4est_geometry_t *cubed_sphere_outer_shell_block_geometry = P4EST_ALLOC(p4est_geometry_t,1);
  p4est_connectivity_t* conn = p8est_connectivity_new_unitcube();

  double* radii = P4EST_ALLOC(double, 2);
  radii[0] = input.R0;
  radii[1] = input.R1;
  
  cubed_sphere_outer_shell_block_geometry->user = (void*)radii;
  cubed_sphere_outer_shell_block_geometry->destroy = d4est_geometry_cubed_sphere_outer_shell_block_destroy;
  cubed_sphere_outer_shell_block_geometry->X = d4est_geometry_cubed_sphere_outer_shell_block_X;

  d4est_geom->p4est_geom = cubed_sphere_outer_shell_block_geometry;
  d4est_geom->p4est_conn = conn;
  
  printf("[GEOMETRY_INFO]: NAME = cubed_sphere_outer_shell_block\n");
  printf("[GEOMETRY_INFO]: R0 = %.25f\n", radii[0]);
  printf("[GEOMETRY_INFO]: R1 = %.25f\n", radii[1]);
}
