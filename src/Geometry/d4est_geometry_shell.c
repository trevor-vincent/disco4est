#include <d4est_geometry.h>
#include <d4est_geometry_shell.h>
#include <ini.h>
#include <p8est_connectivity.h>
#include <p8est_geometry.h>
#include <d4est_util.h>

typedef struct {

  double R1;
  double R2;
  int count;
  
} d4est_geometry_shell_input_t;

static
int d4est_geometry_shell_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_geometry_shell_input_t* pconfig = (d4est_geometry_shell_input_t*)user;
  if (d4est_util_match_couple(section,"geometry",name,"R1")) {
    D4EST_ASSERT(pconfig->R1 == -1);
    pconfig->R1 = atof(value);
    pconfig->count += 1;
  }
  else if (d4est_util_match_couple(section,"geometry",name,"R2")) {
    D4EST_ASSERT(pconfig->R2 == -1);
    pconfig->R2 = atof(value);
    pconfig->count += 1;
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

static
d4est_geometry_shell_input_t
d4est_geometry_shell_input
(
 const char* input_file
)
{
  int num_of_options = 2;
  
  d4est_geometry_shell_input_t input;
  input.count = 0;
  input.R1 = -1;
  input.R2 = -1;
  
  if (ini_parse(input_file, d4est_geometry_shell_input_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_ASSERT(input.count == num_of_options);
  return input;
}


typedef struct
{
  double              R2, R1;
  double              R2byR1, R1sqrbyR2, Rlog;

  p8est_connectivity_t* conn;
  
} d4est_geometry_shell_attr_t;

static void
d4est_geometry_shell_destroy
(
 p8est_geometry_t* geom
)
{
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}


static void
d4est_geometry_compactified_shell_X
(
 p8est_geometry_t * geom,
 p4est_topidx_t which_tree,
 const double rst[3], double xyz[3])
{
  d4est_geometry_shell_attr_t* shell = geom->user;
  double              x, y, R, q;
  double              abc[3];

  /* transform from the reference cube into vertex space */
  d4est_geometry_octree_to_vertex (shell->conn, which_tree, rst, abc);

  /* assert that input points are in the expected range */
  P4EST_ASSERT (0 <= which_tree && which_tree < 24);
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[2] < 2.0 + SC_1000_EPS && abc[2] > 1.0 - SC_1000_EPS);

  /* transform abc[0] and y in-place for nicer grading */
  x = tan (abc[0] * M_PI_4);
  y = tan (abc[1] * M_PI_4);

  /* compute transformation ingredients */
  double m = (2. - 1.)/((1./shell->R2) - (1./shell->R1));
  double t = (1.*shell->R1 - 2.*shell->R2)/(shell->R1 - shell->R2);
  R = m/(abc[2] - t);
  q = R / sqrt (x * x + y * y + 1.);

  /* assign correct coordinates based on patch id */
  switch (which_tree / 4) {
  case 3:                      /* top */
    xyz[0] = +q * y;
    xyz[1] = -q * x;
    xyz[2] = +q;
    break;
  case 2:                      /* left */
    xyz[0] = -q;
    xyz[1] = -q * x;
    xyz[2] = +q * y;
    break;
  case 1:                      /* bottom */
    xyz[0] = -q * y;
    xyz[1] = -q * x;
    xyz[2] = -q;
    break;
  case 0:                      /* right */
    xyz[0] = +q;
    xyz[1] = -q * x;
    xyz[2] = -q * y;
    break;
  case 4:                      /* back */
    xyz[0] = -q * x;
    xyz[1] = +q;
    xyz[2] = +q * y;
    break;
  case 5:                      /* front */
    xyz[0] = +q * x;
    xyz[1] = -q;
    xyz[2] = +q * y;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}


static void
d4est_geometry_shell_X (p8est_geometry_t * geom,
                        p4est_topidx_t which_tree,
                        const double rst[3], double xyz[3])
{
  d4est_geometry_shell_attr_t* shell = geom->user;
  double              x, y, R, q;
  double              abc[3];

  /* transform from the reference cube into vertex space */
  d4est_geometry_octree_to_vertex (shell->conn, which_tree, rst, abc);

  /* assert that input points are in the expected range */
  P4EST_ASSERT (0 <= which_tree && which_tree < 24);
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[2] < 2.0 + SC_1000_EPS && abc[2] > 1.0 - SC_1000_EPS);

  /* transform abc[0] and y in-place for nicer grading */
  x = tan (abc[0] * M_PI_4);
  y = tan (abc[1] * M_PI_4);

  /* compute transformation ingredients */
  R = shell->R1*(2. - abc[2]) + shell->R2*(abc[2] - 1.);  
  q = R / sqrt (x * x + y * y + 1.);

  /* assign correct coordinates based on patch id */
  switch (which_tree / 4) {
  case 3:                      /* top */
    xyz[0] = +q * y;
    xyz[1] = -q * x;
    xyz[2] = +q;
    break;
  case 2:                      /* left */
    xyz[0] = -q;
    xyz[1] = -q * x;
    xyz[2] = +q * y;
    break;
  case 1:                      /* bottom */
    xyz[0] = -q * y;
    xyz[1] = -q * x;
    xyz[2] = -q;
    break;
  case 0:                      /* right */
    xyz[0] = +q;
    xyz[1] = -q * x;
    xyz[2] = -q * y;
    break;
  case 4:                      /* back */
    xyz[0] = -q * x;
    xyz[1] = +q;
    xyz[2] = +q * y;
    break;
  case 5:                      /* front */
    xyz[0] = +q * x;
    xyz[1] = -q;
    xyz[2] = +q * y;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

void
d4est_geometry_shell_new
(
 const char* input_file,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_shell_input_t input = d4est_geometry_shell_input(input_file);
  p8est_geometry_t* shell_geom = P4EST_ALLOC(p8est_geometry_t, 1);
  p8est_connectivity_t* conn = p8est_connectivity_new_shell();
  
  d4est_geometry_shell_attr_t* shell_attrs = P4EST_ALLOC(d4est_geometry_shell_attr_t, 1);

  shell_attrs->conn = conn;
  shell_attrs->R2 =input.R2;
  shell_attrs->R1 =input.R1;
  shell_attrs->R2byR1 =input.R2 /input.R1;
  shell_attrs->R1sqrbyR2 =input.R1 *input.R1 /input.R2;
  shell_attrs->Rlog = log (input.R2 /input.R1);

  shell_geom->name = "shell";
  shell_geom->user = shell_attrs;
  shell_geom->X = d4est_geometry_shell_X;
  shell_geom->destroy = d4est_geometry_shell_destroy;
  
  d4est_geom->p4est_geom = shell_geom;
  d4est_geom->p4est_conn = conn;
  
  printf("[GEOMETRY_INFO]: NAME = shell\n");
  printf("[GEOMETRY_INFO]: R1 = %.25f\n", shell_attrs->R1);
  printf("[GEOMETRY_INFO]: R2 = %.25f\n", shell_attrs->R2);
}

void
d4est_geometry_compactified_shell_new
(
 const char* input_file,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_shell_input_t input = d4est_geometry_shell_input(input_file);
  p8est_geometry_t* shell_geom = P4EST_ALLOC(p8est_geometry_t, 1);
  p4est_connectivity_t* conn = p8est_connectivity_new_shell();
  
  d4est_geometry_shell_attr_t* shell_attrs = P4EST_ALLOC(d4est_geometry_shell_attr_t, 1);

  shell_attrs->conn = conn;
  shell_attrs->R2 =input.R2;
  shell_attrs->R1 =input.R1;
  shell_attrs->R2byR1 =input.R2 /input.R1;
  shell_attrs->R1sqrbyR2 =input.R1 *input.R1 /input.R2;
  shell_attrs->Rlog = log (input.R2 /input.R1);

  shell_geom->name = "compactified_shell";
  shell_geom->user = shell_attrs;
  shell_geom->X = d4est_geometry_compactified_shell_X;
  shell_geom->destroy = d4est_geometry_shell_destroy;
  
  d4est_geom->p4est_geom = shell_geom;
  d4est_geom->p4est_conn = conn;
  
  printf("[GEOMETRY_INFO]: NAME = compactified_shell\n");
  printf("[GEOMETRY_INFO]: R1 = %.25f\n", shell_attrs->R1);
  printf("[GEOMETRY_INFO]: R2 = %.25f\n", shell_attrs->R2);
}
