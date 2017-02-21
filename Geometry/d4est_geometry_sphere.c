#include <pXest.h>
#include <d4est_geometry.h>
#include <d4est_geometry_sphere.h>
#include <ini.h>
#include <p8est_connectivity.h>
#include <util.h>
#include <petscsnes.h>

typedef struct {

  double R0;
  double R1;
  double R2;
  int count;
  
} d4est_geometry_sphere_input_t;

static
int d4est_geometry_sphere_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_geometry_sphere_input_t* pconfig = (d4est_geometry_sphere_input_t*)user;
  if (util_match_couple(section,"geometry",name,"R0")) {
    mpi_assert(pconfig->R0 == -1);
    pconfig->R0 = atof(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"geometry",name,"R1")) {
    mpi_assert(pconfig->R1 == -1);
    pconfig->R1 = atof(value);
    pconfig->count += 1;
  }
  else if (util_match_couple(section,"geometry",name,"R2")) {
    mpi_assert(pconfig->R2 == -1);
    pconfig->R2 = atof(value);
    pconfig->count += 1;
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

d4est_geometry_sphere_input_t
d4est_geometry_sphere_input
(
 const char* input_file
)
{
  int num_of_options = 3;
  
  d4est_geometry_sphere_input_t input;
  input.count = 0;
  input.R0 = -1;
  input.R1 = -1;
  input.R2 = -1;
  
  if (ini_parse(input_file, d4est_geometry_sphere_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  mpi_assert(input.count == num_of_options);
  return input;
}



static void
d4est_geometry_sphere_destroy
(
 p8est_geometry_t* geom
)
{
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}


static void
d4est_geometry_sphere_X(p8est_geometry_t * geom,
                       p4est_topidx_t which_tree,
                       const double rst[3], double xyz[3])
{
  d4est_geometry_sphere_attr_t* sphere = geom->user;

  double              x, y, R, q;
  double              abc[3];

  /* transform from the reference cube into vertex space */
  d4est_geometry_octree_to_vertex (sphere->conn, which_tree, rst, abc);

  /* assert that input points are in the expected range */
  P4EST_ASSERT (0 <= which_tree && which_tree < 13);
  /* if (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS){ */
  /*   printf("which_tree = %d\n", which_tree); */
  /*   printf("abc[0] = %.20f\n", abc[0]); */
  /*   printf("abc[1] = %.20f\n", abc[1]); */
  /*   SC_ABORT_NOT_REACHED(); */
  /* } */
  /* if (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS){ */
  /*   printf("which_tree = %d\n", which_tree); */
  /*   printf("abc[0] = %.20f\n", abc[0]); */
  /*   printf("abc[1] = %.20f\n", abc[1]); */
  /*   SC_ABORT_NOT_REACHED(); */
  /* } */
 
  /* P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS); */
  /* P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS); */

  if (which_tree < 6) {         /* outer shell */
    x = tan (abc[0] * M_PI_4);
    y = tan (abc[1] * M_PI_4);
    /* R = sphere->R1sqrbyR2 * pow (sphere->R2byR1, abc[2]); */
    R = sphere->R1*(2. - abc[2]) + sphere->R2*(abc[2] - 1.);
    q = R / sqrt (x * x + y * y + 1.);
  }
  else if (which_tree < 12) {   /* inner shell */
    double              p, tanx, tany;

    p = 2. - abc[2];
    tanx = tan (abc[0] * M_PI_4);
    tany = tan (abc[1] * M_PI_4);
    x = p * abc[0] + (1. - p) * tanx;
    y = p * abc[1] + (1. - p) * tany;
    /* R = sphere->R0sqrbyR1 * pow (sphere->R1byR0, abc[2]); */
    R = sphere->R0*(2. - abc[2]) + sphere->R1*(abc[2] - 1.);
    
    q = R / sqrt (1. + (1. - p) * (tanx * tanx + tany * tany) + 2. * p);
  }
  else {                        /* center cube */
    xyz[0] = abc[0] * sphere->Clength;
    xyz[1] = abc[1] * sphere->Clength;
    xyz[2] = abc[2] * sphere->Clength;

    return;
  }
  
  /* assign correct coordinates based on direction */
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
d4est_geometry_compactified_sphere_X(p8est_geometry_t * geom,
                                p4est_topidx_t which_tree,
                                const double rst[3], double xyz[3])
{
  d4est_geometry_sphere_attr_t* sphere = geom->user;

  double              x, y, R, q;
  double              abc[3];

  /* transform from the reference cube into vertex space */
  d4est_geometry_octree_to_vertex (sphere->conn, which_tree, rst, abc);

  /* assert that input points are in the expected range */
  P4EST_ASSERT (0 <= which_tree && which_tree < 13);
  /* if (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS){ */
  /*   printf("which_tree = %d\n", which_tree); */
  /*   printf("abc[0] = %.20f\n", abc[0]); */
  /*   printf("abc[1] = %.20f\n", abc[1]); */
  /*   SC_ABORT_NOT_REACHED(); */
  /* } */
  /* if (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS){ */
  /*   printf("which_tree = %d\n", which_tree); */
  /*   printf("abc[0] = %.20f\n", abc[0]); */
  /*   printf("abc[1] = %.20f\n", abc[1]); */
  /*   SC_ABORT_NOT_REACHED(); */
  /* } */
 
  /* P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS); */
  /* P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS); */

  if (which_tree < 6) {         /* outer shell */
    x = tan (abc[0] * M_PI_4);
    y = tan (abc[1] * M_PI_4);
    /* R = sphere->R1sqrbyR2 * pow (sphere->R2byR1, abc[2]); */
    /* R = sphere->R1*(2. - abc[2]) + sphere->R2*(abc[2] - 1.); */


    /* m = (x2 - x1)/((1/r2) - (1/r1)) */
    /* t = (x1*r1 - x2*r2)/(r1 - r2) */
    /* r[x_] := m/(x - t)` */
    
    double m = (2. - 1.)/((1./sphere->R2) - (1./sphere->R1));
    double t = (1.*sphere->R1 - 2.*sphere->R2)/(sphere->R1 - sphere->R2);
    R = m/(abc[2] - t);
    
    q = R / sqrt (x * x + y * y + 1.);
  }
  else if (which_tree < 12) {   /* inner shell */
    double              p, tanx, tany;

    p = 2. - abc[2];
    tanx = tan (abc[0] * M_PI_4);
    tany = tan (abc[1] * M_PI_4);
    x = p * abc[0] + (1. - p) * tanx;
    y = p * abc[1] + (1. - p) * tany;
    /* R = sphere->R0sqrbyR1 * pow (sphere->R1byR0, abc[2]); */
    R = sphere->R0*(2. - abc[2]) + sphere->R1*(abc[2] - 1.);
    
    q = R / sqrt (1. + (1. - p) * (tanx * tanx + tany * tany) + 2. * p);
  }
  else {                        /* center cube */
    xyz[0] = abc[0] * sphere->Clength;
    xyz[1] = abc[1] * sphere->Clength;
    xyz[2] = abc[2] * sphere->Clength;

    return;
  }
  
  /* assign correct coordinates based on direction */
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


void
d4est_geometry_sphere_new
(
 const char* input_file,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_sphere_input_t input = d4est_geometry_sphere_input(input_file);
  p8est_geometry_t* sphere_geom = P4EST_ALLOC(p8est_geometry_t, 1);
  p4est_connectivity_t* conn = p8est_connectivity_new_sphere();
  
  d4est_geometry_sphere_attr_t* sphere_attrs = P4EST_ALLOC(d4est_geometry_sphere_attr_t, 1);
  
  sphere_attrs->R2 = input.R2;
  sphere_attrs->R1 = input.R1;
  sphere_attrs->R0 = input.R0;
  sphere_attrs->conn = conn;
  
  /* variables useful for the outer shell */
  sphere_attrs->R2byR1 = input.R2 / input.R1;
  sphere_attrs->R1sqrbyR2 = input.R1 * input.R1 / input.R2;
  sphere_attrs->R1log = log ( input.R2 / input.R1);

  /* variables useful for the inner shell */
  sphere_attrs->R1byR0 = input.R1 / input.R0;
  sphere_attrs->R0sqrbyR1 = input.R0 * input.R0 / input.R1;
  sphere_attrs->R0log = log ( input.R1 / input.R0);

  /* variables useful for the center cube */
  sphere_attrs->Clength = input.R0 / sqrt (3.);
  sphere_attrs->CdetJ = pow ( input.R0 / sqrt (3.), 3.);

  sphere_geom->name = "sphere";
  sphere_geom->user = sphere_attrs;
  sphere_geom->X = d4est_geometry_sphere_X;
  sphere_geom->destroy = d4est_geometry_sphere_destroy;
  
  d4est_geom->p4est_geom = sphere_geom;
  d4est_geom->p4est_conn = conn;

  sc_MPI_Comm mpicomm = PETSC_COMM_WORLD;  
  int proc_rank;
  MPI_Comm_rank(mpicomm, &proc_rank);

  if (proc_rank == 0){
    printf("[GEOMETRY_INFO]: NAME = sphere\n");
    printf("[GEOMETRY_INFO]: R0 = %.25f\n", sphere_attrs->R0);
    printf("[GEOMETRY_INFO]: R1 = %.25f\n", sphere_attrs->R1);
    printf("[GEOMETRY_INFO]: R2 = %.25f\n", sphere_attrs->R2);
  }
}

void
d4est_geometry_compactified_sphere_new
(
 const char* input_file,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_sphere_input_t input = d4est_geometry_sphere_input(input_file);
  p8est_geometry_t* sphere_geom = P4EST_ALLOC(p8est_geometry_t, 1);
  p4est_connectivity_t* conn = p8est_connectivity_new_sphere();
  d4est_geometry_sphere_attr_t* sphere_attrs = P4EST_ALLOC(d4est_geometry_sphere_attr_t, 1);
  
  sphere_attrs->R2 = input.R2;
  sphere_attrs->R1 = input.R1;
  sphere_attrs->R0 = input.R0;
  sphere_attrs->conn = conn;
  
  /* variables useful for the outer shell */
  sphere_attrs->R2byR1 = input.R2 / input.R1;
  sphere_attrs->R1sqrbyR2 = input.R1 * input.R1 / input.R2;
  sphere_attrs->R1log = log ( input.R2 / input.R1);

  /* variables useful for the inner shell */
  sphere_attrs->R1byR0 = input.R1 / input.R0;
  sphere_attrs->R0sqrbyR1 = input.R0 * input.R0 / input.R1;
  sphere_attrs->R0log = log ( input.R1 / input.R0);

  /* variables useful for the center cube */
  sphere_attrs->Clength = input.R0 / sqrt (3.);
  sphere_attrs->CdetJ = pow ( input.R0 / sqrt (3.), 3.);

  sphere_geom->name = "compact_sphere";
  sphere_geom->user = sphere_attrs;
  sphere_geom->X = d4est_geometry_compactified_sphere_X;
  sphere_geom->destroy = d4est_geometry_sphere_destroy;
  
  d4est_geom->p4est_geom = sphere_geom;
  d4est_geom->p4est_conn = conn;

  sc_MPI_Comm mpicomm = PETSC_COMM_WORLD;  
  int proc_rank;
  MPI_Comm_rank(mpicomm, &proc_rank);

  if (proc_rank == 0){
    printf("[GEOMETRY_INFO]: NAME = compact_sphere\n");
    printf("[GEOMETRY_INFO]: R0 = %.25f\n", sphere_attrs->R0);
    printf("[GEOMETRY_INFO]: R1 = %.25f\n", sphere_attrs->R1);
    printf("[GEOMETRY_INFO]: R2 = %.25f\n", sphere_attrs->R2);
  }
}
/* typically only used for vtk output */
p8est_geometry_t*
d4est_geometry_compactified_sphere_from_param
(
 double R0,
 double R1,
 double R2,
 p8est_connectivity_t* conn
)
{
  p8est_geometry_t* sphere_geom = P4EST_ALLOC(p8est_geometry_t, 1);
  d4est_geometry_sphere_attr_t* sphere_attrs = P4EST_ALLOC(d4est_geometry_sphere_attr_t, 1);
  
  sphere_attrs->R2 = R2;
  sphere_attrs->R1 = R1;
  sphere_attrs->R0 = R0;
  sphere_attrs->conn = conn;
  
  /* variables useful for the outer shell */
  sphere_attrs->R2byR1 = R2 / R1;
  sphere_attrs->R1sqrbyR2 = R1 * R1 / R2;
  sphere_attrs->R1log = log ( R2 / R1);

  /* variables useful for the inner shell */
  sphere_attrs->R1byR0 = R1 / R0;
  sphere_attrs->R0sqrbyR1 = R0 * R0 / R1;
   sphere_attrs->R0log = log ( R1 / R0);

  /* variables useful for the center cube */
  sphere_attrs->Clength = R0 / sqrt (3.);
  sphere_attrs->CdetJ = pow ( R0 / sqrt (3.), 3.);

  sphere_geom->name = "compact_sphere";
  sphere_geom->user = sphere_attrs;
  sphere_geom->X = d4est_geometry_compactified_sphere_X;
  sphere_geom->destroy = d4est_geometry_sphere_destroy;
  
  return sphere_geom;
}
