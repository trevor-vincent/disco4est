#include <d4est_geometry.h>
#include <d4est_geometry_pizza_half.h>
#include <util.h>
#include <ini.h>

#undef P4_TO_P8

typedef struct {

  double R0;
  double R1;
  int count;
  
} d4est_geometry_pizza_half_input_t;

static
int d4est_geometry_pizza_half_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_geometry_pizza_half_input_t* pconfig = (d4est_geometry_pizza_half_input_t*)user;
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
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

d4est_geometry_pizza_half_input_t
d4est_geometry_pizza_half_input
(
 const char* input_file
)
{
  int num_of_options = 2;
  
  d4est_geometry_pizza_half_input_t input;
  input.count = 0;
  input.R0 = -1;
  input.R1 = -1;
  
  if (ini_parse(input_file, d4est_geometry_pizza_half_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  mpi_assert(input.count == num_of_options);
  return input;
}



static
void d4est_geometry_pizza_half_linear_map
(
 double xref,
 double refmin,
 double refmax,
 double emin,
 double emax,
 double *x
){
  *x = emin + (emax-emin)*(xref - refmin)/(refmax - refmin);
}

static
void
d4est_geometry_pizza_half_map_cube_to_slab(
                    double xref,
                    double yref,
                    double cmin,
                    double cmax,
                    double emin,
                    double emax,
                    double* x,
                    double* y
){
  double xbar, ybar;
    /* map x from [0,1] to [emin, emax] */
  d4est_geometry_pizza_half_linear_map(xref, 0., 1., emin, emax, &xbar);
    /* map y from [0,1] to [-1,1] */
  d4est_geometry_pizza_half_linear_map(yref, 0., 1., -1., 1., &ybar);
  double xmin = (1. - cmin)*emin + emin*cmin/sqrt(1. + ybar*ybar);
  double xmax = (1. - cmax)*emax + emax*cmax/sqrt(1. + ybar*ybar);
  *x = xmin + (xmax - xmin)*(xbar - emin)/(emax - emin);
  *y = (*x)*ybar;
}


void 
d4est_geometry_pizza_half_X(p4est_geometry_t * geom,
                       p4est_topidx_t which_tree,
                       const double rst[3], double xyz[3])
{
  double* radii = (double*)geom->user;
  
  double R0 = radii[0];
  double R1 = radii[1];

  double xref = rst[0];
  double yref = rst[1];
  double x,y;
 
  d4est_geometry_pizza_half_map_cube_to_slab(xref, yref, 0, 1, R0/sqrt(2), R1, &x, &y);

  xyz[0] = x - R0/sqrt(2);
  xyz[1] = y;
  xyz[2] = 0.;
}


static void
d4est_geometry_pizza_half_destroy
(
 p4est_geometry_t* geom
)
{
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}


void
d4est_geometry_pizza_half_new
(
 const char* input_file,
 d4est_geometry_t* d4est_geom
)
{
  mpi_assert((P4EST_DIM)==2);
  d4est_geometry_pizza_half_input_t input = d4est_geometry_pizza_half_input(input_file);
  p4est_geometry_t *pizza_half_geometry = P4EST_ALLOC(p4est_geometry_t,1);
  p4est_connectivity_t* conn = NULL;
#if (P4EST_DIM)==2
  conn = p4est_connectivity_new_unitsquare();
#endif
  double* radii = P4EST_ALLOC(double, 2);
  radii[0] = input.R0;
  radii[1] = input.R1;
  
  pizza_half_geometry->user = (void*)radii;
  pizza_half_geometry->destroy = d4est_geometry_pizza_half_destroy;
  pizza_half_geometry->X = d4est_geometry_pizza_half_X;

  d4est_geom->p4est_geom = pizza_half_geometry;
  d4est_geom->p4est_conn = conn;
  
  printf("[GEOMETRY_INFO]: NAME = pizza_half\n");
  printf("[GEOMETRY_INFO]: R0 = %.25f\n", radii[0]);
  printf("[GEOMETRY_INFO]: R1 = %.25f\n", radii[1]);
}
