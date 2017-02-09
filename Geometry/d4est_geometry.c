#include <util.h>
#include <ini.h>
#include <d4est_geometry.h>
#include <d4est_geometry_sphere.h>
#include <d4est_geometry_trap.h>
#include <d4est_geometry_disk.h>
#include <d4est_geometry_pizza_half.h>
#include <d4est_geometry_2pac.h>

typedef struct {

  char* name;
  int count;
  
} d4est_geometry_input_t;


static
int d4est_geometry_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_geometry_input_t* pconfig = (d4est_geometry_input_t*)user;
  if (util_match_couple(section,"geometry",name,"name")) {
    mpi_assert(pconfig->name == NULL);
    D4EST_ASPRINTF(pconfig->name,"%s",value);
    pconfig->count += 1;
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

d4est_geometry_input_t
d4est_geometry_input
(
 const char* input_file
)
{
  int num_of_options = 1;
  
  d4est_geometry_input_t input;
  input.count = 0;
  input.name = NULL;
  
  if (ini_parse(input_file, d4est_geometry_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  mpi_assert(input.count == num_of_options);
  return input;
}


d4est_geometry_t*
d4est_geometry_new(const char* input_file){

  d4est_geometry_input_t input = d4est_geometry_input(input_file);
  d4est_geometry_t* d4est_geom = P4EST_ALLOC(d4est_geometry_t, 1);
  
  if (util_match(input.name,"sphere")) {
    d4est_geometry_sphere_new(input_file, d4est_geom);
  }
  else if (util_match(input.name,"compact_sphere")) {
    d4est_geometry_compactified_sphere_new(input_file, d4est_geom);
  }
  else if (util_match(input.name,"disk")) {
    d4est_geometry_disk_new(input_file, d4est_geom);
  }
  /* These are for debugging purposes */
  else if (util_match(input.name,"2pac_aligned")) {
    d4est_geometry_2pac_aligned_new(input_file, d4est_geom);
  }
  else if (util_match(input.name,"2pac_nonaligned")) {
    d4est_geometry_2pac_nonaligned_new(input_file, d4est_geom);
  }
  else if (util_match(input.name,"2pac_aligned_traps")) {
    d4est_geometry_2pac_aligned_traps_new(input_file, d4est_geom);
  }    
  else if (util_match(input.name,"2pac_aligned_squares")) {
    d4est_geometry_2pac_aligned_squares_new(input_file, d4est_geom);
  }
  else if (util_match(input.name,"trapezoid")) {
    d4est_geometry_trap_new(input_file, d4est_geom);
  }
  else if (util_match(input.name,"pizza_half")) {
    d4est_geometry_pizza_half_new(input_file, d4est_geom);
  }
  else {
    printf("[ERROR]: You tried to use %s geometry\n", input.name);
    mpi_abort("[ERROR]: this geometry is currently not supported");
  }

  free(input.name);
  return d4est_geom;
}

void
d4est_geometry_destroy(d4est_geometry_t* d4est_geom){
  p4est_connectivity_destroy(d4est_geom->p4est_conn);
  p4est_geometry_destroy (d4est_geom->p4est_geom);
  P4EST_FREE(d4est_geom);
}

static void
d4est_geometry_octree_to_vertex (p8est_geometry_t * geom,
                               p4est_topidx_t which_tree,
                               const double abc[3], double xyz[3])
{
  d4est_geometry_sphere_attr_t *sphere_att = (d4est_geometry_sphere_attr_t *) geom->user;
  p4est_connectivity_t *connectivity = (p4est_connectivity_t *) sphere_att->conn;
  const p4est_topidx_t *tree_to_vertex = connectivity->tree_to_vertex;
  const double       *v = connectivity->vertices;
  double              eta_x, eta_y, eta_z = 0.;
  int                 j, k;
  p4est_topidx_t      vt[P4EST_CHILDREN];

  /* retrieve corners of the tree */
  for (k = 0; k < P4EST_CHILDREN; ++k) {
    vt[k] = tree_to_vertex[which_tree * P4EST_CHILDREN + k];
  }

  /* these are reference coordinates in [0, 1]**d */
  eta_x = abc[0];
  eta_y = abc[1];
  eta_z = abc[2];

  /* bi/trilinear transformation */
  for (j = 0; j < 3; ++j) {
    /* *INDENT-OFF* */
    xyz[j] =
           ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                  eta_x  * v[3 * vt[1] + j]) +
                                  eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                                  eta_x  * v[3 * vt[3] + j]))
#ifdef P4_TO_P8
            +     eta_z  * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[4] + j] +
                                                  eta_x  * v[3 * vt[5] + j]) +
                                  eta_y  * ((1. - eta_x) * v[3 * vt[6] + j] +
                                                  eta_x  * v[3 * vt[7] + j]))
#endif
           );
    /* *INDENT-ON* */
  }
}
/* p4est_connectivity_t* */
/* problem_build_conn() */
/* {   */
/*   p4est_connectivity_t* conn; */
/* #if (P4EST_DIM)==2    */
/* #if defined DISK */
/*   conn = p4est_connectivity_new_disk(); */
/* #elif defined SQUARE */
/*   conn = p4est_connectivity_new_unitsquare(); */
/* #elif defined TUPAC_ALIGNED */
/*   conn = p4est_connectivity_new_2pac_aligned(); */
/* #elif defined TUPAC_ALIGNED_CUBE */
/*   conn = p4est_connectivity_new_2pac_aligned_CUBE(); */
/* #elif defined TUPAC_NONALIGNED */
/*   conn = p4est_connectivity_new_2pac_nonaligned(); */
/* #elif defined TUPAC_NONALIGNED_CUBE */
/*   conn = p4est_connectivity_new_2pac_nonaligned_CUBE(); */
/* #elif defined TRAP */
/*   conn = p4est_connectivity_new_trap(); */
/* #elif defined RECTANGLE */
/*   conn = p4est_connectivity_new_unitsquare(); */
/* #elif defined PIZZA_HALF */
/*   conn = p4est_connectivity_new_unitsquare();   */
/* #else */
/*   mpi_abort("Pick a 2d connectivity and geometry"); */
/* #endif */
/* #endif */

/* #if (P4EST_DIM)==3 */
/* #if defined SPHERE */
/*   conn = p8est_connectivity_new_sphere(); */
/* #elif defined CUBE */
/*   conn = p8est_connectivity_new_unitcube(); */
/* #else */
/*   mpi_abort("Pick a 3d connectivity and geometry"); */
/* #endif */
/* #endif */
/*   return conn; */
/* } */

/* p4est_geometry_t* */
/* problem_build_geom */
/* ( */
/*  p4est_connectivity_t* conn */
/* ) */
/* { */
/*   p4est_geometry_t* geom; */
/*   /\* double x0 = -1.; *\/ */
/*   /\* double x1 = 1.; *\/ */
/*   /\* double y0 = 0.; *\/ */
/*   /\* double y1 = 2.; *\/ */

  
/* #if (P4EST_DIM)==2    */
/* #if defined DISK */
/*   printf("hi from disk!"); */
/*   geom = d4est_geometry_new_disk(conn, R*.5, R); */
/* #elif defined SQUARE */
/*   geom = p4est_geometry_new_connectivity(conn); */
/* #elif defined RECTANGLE */
/*   geom = d4est_geometry_new_rectangle(conn, x0, x1, y0, y1); */
/* #elif defined TUPAC_ALIGNED */
/*   /\* geom = p4est_geometry_new_connectivity(conn); *\/ */
/*   /\* geom = p4est_geometry_new_connectivity(conn); *\/ */
/*   geom = d4est_geometry_new_2pac_aligned(conn, R*.5, R); */
/*   /\* geom = d4est_geometry_new_2pac_aligned_traps(conn, R*.5, R); *\/ */
/*   /\* geom = d4est_geometry_new_2pac_aligned_squares(conn, R*.5, R); *\/ */
/* #elif defined TUPAC_ALIGNED_CUBE */
/*   geom = p4est_geometry_new_connectivity(conn); */
/* #elif defined TUPAC_NONALIGNED */
/*   geom = p4est_geometry_new_connectivity(conn); */
/* #elif defined TUPAC_NONALIGNED_CUBE */
/*   geom = p4est_geometry_new_connectivity(conn); */
/* #elif defined TRAP */
/*   geom = p4est_geometry_new_connectivity(conn); */
/* #elif defined PIZZA_HALF */
/*   printf("hi from pizza half\n"); */
/*   geom = d4est_geometry_new_pizza_half(conn, R*.5, R); */
/* #else */
/*   mpi_abort("Pick a 2d connectivity and geometry"); */
/* #endif */
/* #endif */
