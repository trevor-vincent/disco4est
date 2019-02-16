#include <pXest.h>
#include <d4est_geometry.h>
#include <ini.h>
#include <p8est_connectivity.h>
#include <d4est_util.h>
#include <d4est_geometry_cubed_sphere_vtk.h>
#include <d4est_connectivity_cubed_sphere.h>
#include <zlog.h>

static
int d4est_geometry_cubed_sphere_vtk_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_geometry_cubed_sphere_vtk_attr_t* pconfig = user;
  if (d4est_util_match_couple(section,pconfig->input_section,name,"R0")) {
    D4EST_ASSERT(pconfig->R0 == -1);
    pconfig->R0 = atof(value);
    D4EST_ASSERT(pconfig->R0 > 0);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"outer_angle_multiplier")) {
    D4EST_ASSERT(pconfig->outer_angle_multiplier == -1);
    pconfig->outer_angle_multiplier = atof(value);
    D4EST_ASSERT(pconfig->outer_angle_multiplier > 0);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"inner_angle_multiplier")) {
    D4EST_ASSERT(pconfig->inner_angle_multiplier == -1);
    pconfig->inner_angle_multiplier = atof(value);
    D4EST_ASSERT(pconfig->inner_angle_multiplier > 0);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"inner_nonz_multiplier")) {
    D4EST_ASSERT(pconfig->inner_nonz_multiplier == -1);
    pconfig->inner_nonz_multiplier = atof(value);
    D4EST_ASSERT(pconfig->inner_nonz_multiplier > 0);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"R2_a")) {
    D4EST_ASSERT(pconfig->R2_a == -1);
    pconfig->R2_a = atof(value);
    D4EST_ASSERT(pconfig->R2_a > 0);
  }
    else if (d4est_util_match_couple(section,pconfig->input_section,name,"R2_b")) {
    D4EST_ASSERT(pconfig->R2_b == -1);
    pconfig->R2_b = atof(value);
    D4EST_ASSERT(pconfig->R2_b > 0);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"R1_a")) {
    D4EST_ASSERT(pconfig->R1_a == -1);
    pconfig->R1_a = atof(value);
    D4EST_ASSERT(pconfig->R1_a > 0);
  }
    else if (d4est_util_match_couple(section,pconfig->input_section,name,"R1_b")) {
    D4EST_ASSERT(pconfig->R1_b == -1);
    pconfig->R1_b = atof(value);
    D4EST_ASSERT(pconfig->R1_b > 0);
  }  
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


static
d4est_geometry_cubed_sphere_vtk_attr_t*
d4est_geometry_cubed_sphere_vtk_input
(
 const char* input_file,
 const char* input_section
)
{
  d4est_geometry_cubed_sphere_vtk_attr_t* sphere_attrs = P4EST_ALLOC(d4est_geometry_cubed_sphere_vtk_attr_t, 1);
  snprintf (sphere_attrs->input_section, sizeof(sphere_attrs->input_section), "%s", input_section);
  sphere_attrs->R0 = -1;
  sphere_attrs->R1_a = -1;
  sphere_attrs->R1_b = -1;
  sphere_attrs->R2_a = -1;
  sphere_attrs->R2_b = -1;
  sphere_attrs->outer_angle_multiplier = -1;
  sphere_attrs->inner_angle_multiplier = -1;
  sphere_attrs->inner_nonz_multiplier = -1;

  if (ini_parse(input_file, d4est_geometry_cubed_sphere_vtk_input_handler, sphere_attrs) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, sphere_attrs->R0, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->outer_angle_multiplier, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->inner_angle_multiplier, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->inner_nonz_multiplier, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->R1_a, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->R1_b, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->R2_a, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->R2_b, -1);

  /* variables useful for the center cube */
  sphere_attrs->Clength = sphere_attrs->R0 / sqrt (3.);

  return sphere_attrs;
}

static
int d4est_geometry_cubed_sphere_vtk_get_number_of_regions
(
 d4est_geometry_t* d4est_geom
){
  if (d4est_geom->geom_type == GEOM_CUBED_SPHERE_VTK){
    return 3;
  }
  else {
    D4EST_ABORT("Geometry not supported yet");
    return -1;
  }
}

static
int d4est_geometry_cubed_sphere_vtk_get_region
(
 d4est_geometry_t* d4est_geom,
 p4est_qcoord_t q [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int tree
){
  if (d4est_geom->geom_type == GEOM_CUBED_SPHERE_VTK){
    if (tree < 6) {         /* outer shell */
      return 0;
    }
    else if (tree < 12) {   /* inner shell */
      return 1;
    }
    else {                  /* center cube */
      return 2;
    }
  }
  else {
    D4EST_ABORT("Geometry not supported yet");
  }
}

static void
d4est_geometry_cubed_sphere_vtk_X(
                              d4est_geometry_t * geom,
                              p4est_topidx_t which_tree,
                              p4est_qcoord_t q0 [3],
                              p4est_qcoord_t dq,
                              const double coords[3],
                              coords_type_t coords_type,
                              double xyz[3]
)
{
  d4est_geometry_cubed_sphere_vtk_attr_t* sphere = geom->user;

  double inner_angle_multiplier = sphere->inner_angle_multiplier;
  double inner_nonz_multiplier = sphere->inner_nonz_multiplier;
  double outer_angle_multiplier = sphere->outer_angle_multiplier;
  double              x, y, R, q;
  double              abc[3];

  double tcoords [3];
  d4est_geometry_get_tree_coords_in_range_0_to_1(q0, dq, coords, coords_type, tcoords);
  
  /* transform from the reference cube into vertex space */
  d4est_geometry_octree_to_vertex (geom->p4est_conn, which_tree, tcoords, abc);

  /* assert that input points are in the expected range */
  P4EST_ASSERT (0 <= which_tree && which_tree < 13);
  if (which_tree < 6) {         /* outer shell */
    x = tan (abc[0] * M_PI_4 * outer_angle_multiplier);
    y = tan (abc[1] * M_PI_4 * outer_angle_multiplier);
    R = sphere->R2_a*(2. - abc[2]) + sphere->R2_b*(abc[2] - 1.);
    q = R / sqrt (x * x + y * y + 1.);
  }
  else if (which_tree < 12) {   /* inner shell */
    double              p, tanx, tany;
    p = 2. - abc[2];
    tanx = tan (abc[0] * M_PI_4 * inner_angle_multiplier);
    tany = tan (abc[1] * M_PI_4 * inner_angle_multiplier);
    x = p * abc[0] + (1. - p) * tanx;
    y = p * abc[1] + (1. - p) * tany;

    /* if (abc[2] == 0){ */
    x*= inner_nonz_multiplier;
    y*= inner_nonz_multiplier;
    /* } */
    
    R = sphere->R1_a*(2. - abc[2]) + sphere->R1_b*(abc[2] - 1);
    q = R / sqrt (1. + (1. - p) * (tanx * tanx + tany * tany) + 2. * p);
  }
  else {                        /* center cube */
    xyz[0] = abc[0] * sphere->Clength;
    xyz[1] = abc[1] * sphere->Clength;
    xyz[2] = abc[2] * sphere->Clength;

    return;
  }
  switch (which_tree % 6) {
  case 0:                      /* front */
    xyz[0] = +q * x;
    xyz[1] = -q;
    xyz[2] = +q * y;
    break;
  case 1:                      /* top */
    xyz[0] = +q * x;
    xyz[1] = +q * y;
    xyz[2] = +q;
    break;
  case 2:                      /* back */
    xyz[0] = +q * x;
    xyz[1] = +q;
    xyz[2] = -q * y;
    break;
  case 3:                      /* right */
    xyz[0] = +q;
    xyz[1] = -q * x;
    xyz[2] = -q * y;
    break;
  case 4:                      /* bottom */
    xyz[0] = -q * y;
    xyz[1] = -q * x;
    xyz[2] = -q;
    break;
  case 5:                      /* left */
    xyz[0] = -q;
    xyz[1] = -q * x;
    xyz[2] = +q * y;
    break;
  default:
    SC_ABORT_NOT_REACHED();
  }
}

static void
d4est_geometry_cubed_sphere_vtk_destroy
(
 d4est_geometry_t* geom
)
{
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}

void
d4est_geometry_cubed_sphere_vtk_new
(
 int mpirank,
 const char* input_file,
 const char* input_section,
 zlog_category_t *c_default,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_cubed_sphere_vtk_attr_t* sphere_attrs = d4est_geometry_cubed_sphere_vtk_input(input_file, input_section);
  p4est_connectivity_t* conn = p8est_connectivity_new_sphere();
  
  d4est_geom->user = sphere_attrs;
  d4est_geom->X = d4est_geometry_cubed_sphere_vtk_X;
  d4est_geom->DX = NULL;
  d4est_geom->D2X = NULL;
  d4est_geom->JAC = NULL;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_vtk_destroy;
  d4est_geom->p4est_conn = conn;
  d4est_geom->get_number_of_regions = d4est_geometry_cubed_sphere_vtk_get_number_of_regions;
  d4est_geom->get_region = d4est_geometry_cubed_sphere_vtk_get_region;
  
  if (mpirank == 0){
    zlog_debug(c_default, "NAME = cubed sphere");
    zlog_debug(c_default, "R0 = %.25f", sphere_attrs->R0);
    zlog_debug(c_default, "outer_angle_multiplier = %.25f", sphere_attrs->outer_angle_multiplier);
    zlog_debug(c_default, "inner_angle_multiplier = %.25f", sphere_attrs->inner_angle_multiplier);
    zlog_debug(c_default, "inner_nonz_multiplier = %.25f", sphere_attrs->inner_nonz_multiplier);
    zlog_debug(c_default, "R1_a = %.25f", sphere_attrs->R1_a);
    zlog_debug(c_default, "R1_b = %.25f", sphere_attrs->R1_b);
    zlog_debug(c_default, "R2_a = %.25f", sphere_attrs->R2_a);
    zlog_debug(c_default, "R2_b = %.25f", sphere_attrs->R2_b);
  }
}
