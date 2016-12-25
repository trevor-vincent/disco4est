
#include <d4est_geometry_compact_disk.h>

static
void d4est_geometry_compact_disk_linear_map
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
d4est_geometry_compact_disk_map_cube_to_slab(
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
  d4est_geometry_compact_disk_linear_map(xref, 0., 1., emin, emax, &xbar);
    /* map y from [0,1] to [-1,1] */
  d4est_geometry_compact_disk_linear_map(yref, 0., 1., -1., 1., &ybar);
  double xmin = (1. - cmin)*emin + emin*cmin/sqrt(1. + ybar*ybar);
  double xmax = (1. - cmax)*emax + emax*cmax/sqrt(1. + ybar*ybar);
  *x = xmin + (xmax - xmin)*(xbar - emin)/(emax - emin);
  *y = (*x)*ybar;
}


static void 
d4est_geometry_compact_disk_X(p4est_geometry_t * geom,
                       p4est_topidx_t which_tree,
                       const double rst[3], double xyz[3])
{
  double* radii = (double*)geom->user;
  
  double R0 = radii[0];
  double R1 = radii[1];

  double xref = rst[0];
  double yref = rst[1];
  double x,y;
  
  if (which_tree == 0){
    /* bottom */
    d4est_geometry_compact_disk_map_cube_to_slab(yref, xref, 1., 0., -R1, -R0/sqrt(2), &y, &x);
    x *= -1.;
  }
  else if (which_tree == 1){
    /* left */
    d4est_geometry_compact_disk_map_cube_to_slab(xref, yref, 1., 0., -R1, -R0/sqrt(2), &x, &y);
    y *= -1;
  }
  else if (which_tree == 2){
    /* center */
    d4est_geometry_compact_disk_linear_map(xref, 0., 1., -R0/sqrt(2), R0/sqrt(2), &x);
    d4est_geometry_compact_disk_linear_map(yref, 0., 1., -R0/sqrt(2), R0/sqrt(2), &y);
  }
  else if (which_tree == 3){
    /* right */
    d4est_geometry_compact_disk_map_cube_to_slab(xref, yref, 0, 1, R0/sqrt(2), R1, &x, &y);
  }
  else if (which_tree == 4){
    /* top */
    d4est_geometry_compact_disk_map_cube_to_slab(yref, xref, 0, 1, R0/sqrt(2), R1, &y, &x);
  }
  else{
    SC_ABORT_NOT_REACHED();
  }
z
  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = 0.;

  /* compactify */
  if (which_tree != 2){
    double xnc = xyz[0];
    double ync = xyz[1];
    double R = sqrt(xnc*xnc + ync*ync);
    double r = (R1 + (w - R1)*(R - R1)/(R2 - R1))/(1 - (R - R1)*(1 - w/Rinf)/(R2 - R1));
    xyz[0] = r*xnc/R;
    xyz[1] = r*ync/R;
  }
}

static void
d4est_geometry_compact_disk_destroy
(
 p4est_geometry_t* geom
)
{
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}


p4est_geometry_t*
d4est_geometry_new_compact_disk
(
 p4est_connectivity_t * conn,
 double R0,
 double R1,
 double w,
 double Rinf,
)
{
  p4est_geometry_t *compact_disk_geometry = P4EST_ALLOC(p4est_geometry_t,1);

  double* radii = P4EST_ALLOC(double, 4);
  radii[0] = R0;
  radii[1] = R1;
  radii[2] = w;
  radii[3] = Rinf;
  
  compact_disk_geometry->user = (void*)radii;
  compact_disk_geometry->destroy = d4est_geometry_compact_disk_destroy;
  compact_disk_geometry->X = d4est_geometry_compact_disk_X;

  return compact_disk_geometry;
}
