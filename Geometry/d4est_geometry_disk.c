#include <ini.h>
#include <util.h>
#include <pXest.h>
#include <d4est_geometry.h>

#if (P4EST_DIM)==2
#include <d4est_geometry_disk.h>
#include <p4est_connectivity.h>



typedef struct {

  double R0;
  double R1;
  
} d4est_geometry_disk_input_t;

static
int d4est_geometry_disk_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_geometry_disk_input_t* pconfig = (d4est_geometry_disk_input_t*)user;
  if (util_match_couple(section,"geometry",name,"R0")) {
    mpi_assert(pconfig->R0 == -1);
    pconfig->R0 = atof(value);
  }
  else if (util_match_couple(section,"geometry",name,"R1")) {
    mpi_assert(pconfig->R1 == -1);
    pconfig->R1 = atof(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

d4est_geometry_disk_input_t*
d4est_geometry_disk_input
(
 const char* input_file
)
{
  d4est_geometry_disk_input_t* input = P4EST_ALLOC(d4est_geometry_disk_input_t, 1);
  input->R0 = -1;
  input->R1 = -1;
  
  if (ini_parse(input_file, d4est_geometry_disk_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  D4EST_CHECK_INPUT("geometry", input->R0, -1);
  D4EST_CHECK_INPUT("geometry", input->R1, -1);
  
  return input;
}


static
void d4est_geometry_disk_linear_map
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
d4est_geometry_disk_map_cube_to_slab(
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
  d4est_geometry_disk_linear_map(xref, 0., 1., emin, emax, &xbar);
    /* map y from [0,1] to [-1,1] */
  d4est_geometry_disk_linear_map(yref, 0., 1., -1., 1., &ybar);
  double xmin = (1. - cmin)*emin + emin*cmin/sqrt(1. + ybar*ybar);
  double xmax = (1. - cmax)*emax + emax*cmax/sqrt(1. + ybar*ybar);
  *x = xmin + (xmax - xmin)*(xbar - emin)/(emax - emin);
  *y = (*x)*ybar;
}


static void 
d4est_geometry_disk_X(d4est_geometry_t * geom,
                      p4est_topidx_t which_tree,
                      p4est_qcoord_t q0 [2],
                      p4est_qcoord_t dq,
                      const double coords[2],
                      coords_type_t coords_type,
                      double xyz[2]
                      )
{
  double* radii = (double*)geom->user;
  
  double R0 = radii[0];
  double R1 = radii[1];

  double tcoords[2];

  d4est_geometry_get_tree_coords_in_range_0_to_1(q0, dq, coords, coords_type, tcoords);
  
  double xref = tcoords[0];
  double yref = tcoords[1];
  double x,y;
  
  if (which_tree == 0){
    /* bottom */
    d4est_geometry_disk_map_cube_to_slab(yref, xref, 1., 0., -R1, -R0/sqrt(2), &y, &x);
    x *= -1.;
  }
  else if (which_tree == 1){
    /* left */
    d4est_geometry_disk_map_cube_to_slab(xref, yref, 1., 0., -R1, -R0/sqrt(2), &x, &y);
    y *= -1;
  }
  else if (which_tree == 2){
    /* center */
    d4est_geometry_disk_linear_map(xref, 0., 1., -R0/sqrt(2), R0/sqrt(2), &x);
    d4est_geometry_disk_linear_map(yref, 0., 1., -R0/sqrt(2), R0/sqrt(2), &y);
  }
  else if (which_tree == 3){
    /* right */
    d4est_geometry_disk_map_cube_to_slab(xref, yref, 0, 1, R0/sqrt(2), R1, &x, &y);
  }
  else if (which_tree == 4){
    /* top */
    d4est_geometry_disk_map_cube_to_slab(yref, xref, 0, 1, R0/sqrt(2), R1, &y, &x);
  }
  else{
    SC_ABORT_NOT_REACHED();
  }

  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = 0.;
}

static void
d4est_geometry_disk_destroy
(
 d4est_geometry_t* geom
)
{
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}


void
d4est_geometry_disk_new
(
 int mpirank,
 const char* input_file,
 d4est_geometry_t* d4est_geom
)
{
  mpi_assert((P4EST_DIM)==2);
  d4est_geometry_disk_input_t* input = d4est_geometry_disk_input(input_file);
  p4est_connectivity_t* conn = p4est_connectivity_new_disk();
  
  d4est_geom->user = input;
  d4est_geom->destroy = d4est_geometry_disk_destroy;
  d4est_geom->X = d4est_geometry_disk_X;
  d4est_geom->DX = NULL;

  if (mpirank == 0){
    printf("[GEOMETRY_INFO]: NAME = disk\n");
    printf("[GEOMETRY_INFO]: R0 = %.25f\n", input->R0);
    printf("[GEOMETRY_INFO]: R1 = %.25f\n", input->R1);
  }
}

#endif
