#define _GNU_SOURCE
#include <pXest.h>
#include <d4est_geometry.h>
#include <d4est_geometry_general_wedge.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_geometry_hole_in_a_box.h>
#include <d4est_connectivity_cubed_sphere.h>
#include <ini.h>

typedef struct {

  double zmin;
  double zmax;
  double inner_radius;
  double box_length;
  char* input_section;
  
} d4est_geometry_hole_in_a_box_attr_t;


static
int d4est_geometry_hole_in_a_box_get_number_of_regions
(
 d4est_geometry_t* d4est_geom
){
  return 1;
}

static
int d4est_geometry_hole_in_a_box_get_region
(
 d4est_geometry_t* d4est_geom,
 p4est_qcoord_t q [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int tree
){
  return 0;
}

static
int d4est_geometry_hole_in_a_box_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_geometry_hole_in_a_box_attr_t* pconfig = user;
  if (d4est_util_match_couple(section,pconfig->input_section,name,"inner_radius")) {
    D4EST_ASSERT(pconfig->inner_radius == -1);
    pconfig->inner_radius = atof(value);
    D4EST_ASSERT(pconfig->inner_radius > 0);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"box_length")) {
    D4EST_ASSERT(pconfig->box_length == -1);
    pconfig->box_length = atof(value);
    D4EST_ASSERT(pconfig->box_length > 0);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

static
d4est_geometry_hole_in_a_box_attr_t*
d4est_geometry_hole_in_a_box_input
(
 const char* input_file,
 const char* input_section
)
{
  d4est_geometry_hole_in_a_box_attr_t* hiab_attrs = P4EST_ALLOC(d4est_geometry_hole_in_a_box_attr_t, 1);
  asprintf(&hiab_attrs->input_section,"%s",input_section);
  hiab_attrs->inner_radius = -1;
  hiab_attrs->box_length = -1;

  if (ini_parse(input_file, d4est_geometry_hole_in_a_box_input_handler, hiab_attrs) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, hiab_attrs->inner_radius, -1);
  D4EST_CHECK_INPUT(input_section, hiab_attrs->box_length, -1);

  return hiab_attrs;
}

static void
d4est_geometry_hole_in_a_box_DX
(
 d4est_geometry_t* d4est_geom,
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 const double rst[(P4EST_DIM)], 
 double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]
)
{
  d4est_geometry_hole_in_a_box_attr_t* hiab_attrs = d4est_geom->user;
  double dxyz_drst_top[(P4EST_DIM)][(P4EST_DIM)];

  d4est_geometry_general_wedge_noncompactified_DX
    (
     d4est_geom,
     which_tree,
     q0,
     dq,
     rst, /* [-1,1]^3 */
     dxyz_drst_top,
     1,
     0,
     hiab_attrs->zmin,
     hiab_attrs->zmax,
     0,
     FULL_WEDGE
    );


  d4est_geometry_cubed_sphere_DX_aux_rotate
    (
     which_tree,
     dxyz_drst_top,
     dxyz_drst
    );
}

static void
d4est_geometry_hole_in_a_box_X
(
 d4est_geometry_t * geom,
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [3],
 p4est_qcoord_t dq,
 const double coords[3],
 coords_type_t coords_type,
 double xyz[3]
)
{

  d4est_geometry_hole_in_a_box_attr_t* hiab_attrs = geom->user;
  double xyz_top[(P4EST_DIM)];

  d4est_geometry_general_wedge_X
    (
     geom,
     which_tree,
     q0,
     dq,
     coords,
     coords_type,
     xyz_top,
     1,
     0,
     hiab_attrs->zmin,
     hiab_attrs->zmax,
     0,
     FULL_WEDGE
    );

  d4est_geometry_cubed_sphere_X_aux_rotate
    (
     which_tree,
     xyz_top,
     xyz
    );
}

static void
d4est_geometry_hole_in_a_box_destroy
(
 d4est_geometry_t* geom
)
{
  free(((d4est_geometry_hole_in_a_box_attr_t*)geom->user)->input_section);
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}

void
d4est_geometry_hole_in_a_box_new
(
 int mpirank,
 const char* input_file,
 const char* input_section,
 const char* printf_prefix,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_hole_in_a_box_attr_t* hiab_attrs = d4est_geometry_hole_in_a_box_input(input_file, input_section);
  p4est_connectivity_t* conn = d4est_connectivity_new_sphere_with_hole();

  hiab_attrs->zmin = hiab_attrs->inner_radius/sqrt(3);
  hiab_attrs->zmax = hiab_attrs->box_length/2;
  
  d4est_geom->p4est_conn = conn; 
  d4est_geom->user = hiab_attrs;
  d4est_geom->X = d4est_geometry_hole_in_a_box_X;
  d4est_geom->DX = d4est_geometry_hole_in_a_box_DX; 
  d4est_geom->JAC = NULL;
  d4est_geom->destroy = d4est_geometry_hole_in_a_box_destroy;
  d4est_geom->get_number_of_regions = d4est_geometry_hole_in_a_box_get_number_of_regions;
  d4est_geom->get_region = d4est_geometry_hole_in_a_box_get_region;
  
  if (mpirank == 0){
    printf("%s: NAME = hole in a box\n", printf_prefix );
    printf("%s: inner radius = %.25f\n", printf_prefix , hiab_attrs->inner_radius);
    printf("%s: box length = %.25f\n", printf_prefix , hiab_attrs->box_length);
    printf("%s: zmin = %.25f\n", printf_prefix , hiab_attrs->zmin);
    printf("%s: zmax = %.25f\n", printf_prefix , hiab_attrs->zmax);
  }
}
