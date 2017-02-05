#include <d4est_geometry_pizza_half.h>

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



void 
d4est_geometry_pizza_half_dxdr2(p4est_geometry_t * geom,
                                p4est_topidx_t which_tree,
                                const double rst[(P4EST_DIM)], // \in [-1,1] such as GL or GLL points
                                double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)],
                                p4est_qcoord_t q0 [(P4EST_DIM)],
                                p4est_qcoord_t dq)
{
  
  double* radii = (double*)geom->user;
  
  double R0 = radii[0];
  double R1 = radii[1];
  
  double r = rst[0];
  double s = rst[1];
  /* double x,y; */

  double amin = q0[0]; amin /= (double)(P4EST_ROOT_LEN);
  double amax = q0[0] + dq; amax /= (double)(P4EST_ROOT_LEN);
  double bmin = q0[1]; bmin /= (double)(P4EST_ROOT_LEN);
  double bmax = q0[1] + dq; bmax /= (double)(P4EST_ROOT_LEN);
  
  double cmin, cmax, emin, emax;  
  cmin = 0.;
  cmax = 1.;
  emin = R0/sqrt(2);
  emax = R1;
  
  dxyz_drst[0][0] = ((amax - amin)*(-1.*(-1. + cmax)*emax + (-1. + cmin)*emin + 
       (cmax*emax)/sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2)) - 
                                    (cmin*emin)/sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2))))/2.;

  dxyz_drst[0][1] = -((-bmax + bmin)*(2*cmin*emin - amin*(cmax*emax - cmin*emin)*(-1 + r) + 
        amax*(cmax*emax - cmin*emin)*(1 + r))*(1 + bmin*(-1 + s) - bmax*(1 + s)))/
                    (2.*pow(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2),1.5));

  dxyz_drst[1][0] = ((amax - amin)*(-1 + bmax + bmin + bmax*s - bmin*s)*
     (-1.*(-1. + cmax)*emax + (-1. + cmin)*emin + 
       (cmax*emax)/sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2)) - 
      (cmin*emin)/sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2))))/2.;
  
  dxyz_drst[1][1] = ((bmax - bmin)*(-(((2*cmin*emin - amin*(cmax*emax - cmin*emin)*(-1 + r) + 
              amax*(cmax*emax - cmin*emin)*(1 + r))*
            pow(-1 + bmax + bmin + bmax*s - bmin*s,2))/
          (sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2))*
            (2. + 2*bmin*(-1 + s) + pow(bmin,2)*pow(-1 + s,2) + 
              pow(bmax,2)*pow(1 + s,2) - 2*bmax*(1 + s + bmin*(-1 + pow(s,2)))))) + 
       2*(-1.*(-1. + cmin)*emin + (cmin*emin)/
           sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2)) + 
          ((amax + amin + amax*r - amin*r)*
             (-1.*(-1. + cmax)*emax + (-1. + cmin)*emin + 
               (cmax*emax)/sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2)) - 
               (cmin*emin)/sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2))))/2.)))/
                    2.;

}

void 
d4est_geometry_pizza_half_dxdr_face(
                                    p4est_geometry_t * geom,
                                    p4est_topidx_t which_tree,
                                    const double rs[(P4EST_DIM)-1], // \in [-1,1] such as GL or GLL points
                                    double dxyz_drs[(P4EST_DIM)][(P4EST_DIM)-1],
                                    p4est_qcoord_t q0 [(P4EST_DIM)],
                                    p4est_qcoord_t dqa [(P4EST_DIM)-1][(P4EST_DIM)]
                                   )
{
  
  double* radii = (double*)geom->user;
  
  double R0 = radii[0];
  double R1 = radii[1];
  
  double t = rs[0];

  double amin = q0[0]; amin /= (double)(P4EST_ROOT_LEN);
  double amax = q0[0] + dqa[0][0]; amax /= (double)(P4EST_ROOT_LEN);
  double bmin = q0[1]; bmin /= (double)(P4EST_ROOT_LEN);
  double bmax = q0[1] + dqa[0][1]; bmax /= (double)(P4EST_ROOT_LEN);
  
  double cmin, cmax, emin, emax;  
  cmin = 0.;
  cmax = 1.;
  emin = R0/sqrt(2);
  emax = R1;

  dxyz_drs[0][0] = ((-2*(bmax - bmin)*cmin*emin*(-1 + bmax + bmin + bmax*t - bmin*t))/pow(1 + pow(-1 + bmax + bmin + bmax*t - bmin*t,2),1.5) - 
     ((bmax - bmin)*(cmax*emax - cmin*emin)*(amax + amin + amax*t - amin*t)*(-1 + bmax + bmin + bmax*t - bmin*t))/pow(1 + pow(-1 + bmax + bmin + bmax*t - bmin*t,2),1.5) + 
     (amax - amin)*(-1.*(-1. + cmax)*emax + (-1. + cmin)*emin + (cmax*emax)/sqrt(1 + pow(-1 + bmax + bmin + bmax*t - bmin*t,2)) - (cmin*emin)/sqrt(1 + pow(-1 + bmax + bmin + bmax*t - bmin*t,2))))
                   /2.;

  dxyz_drs[1][0] = ((-1 + bmax + bmin + bmax*t - bmin*t)*((-2*(bmax - bmin)*cmin*emin*(-1 + bmax + bmin + bmax*t - bmin*t))/pow(1 + pow(-1 + bmax + bmin + bmax*t - bmin*t,2),1.5) - 
        ((bmax - bmin)*(cmax*emax - cmin*emin)*(amax + amin + amax*t - amin*t)*(-1 + bmax + bmin + bmax*t - bmin*t))/pow(1 + pow(-1 + bmax + bmin + bmax*t - bmin*t,2),1.5) + 
        (amax - amin)*(-1.*(-1. + cmax)*emax + (-1. + cmin)*emin + (cmax*emax)/sqrt(1 + pow(-1 + bmax + bmin + bmax*t - bmin*t,2)) - 
           (cmin*emin)/sqrt(1 + pow(-1 + bmax + bmin + bmax*t - bmin*t,2)))))/2. + 
   (bmax - bmin)*(-1.*(-1. + cmin)*emin + (cmin*emin)/sqrt(1 + pow(-1 + bmax + bmin + bmax*t - bmin*t,2)) + 
      ((amax + amin + amax*t - amin*t)*(-1.*(-1. + cmax)*emax + (-1. + cmin)*emin + (cmax*emax)/sqrt(1 + pow(-1 + bmax + bmin + bmax*t - bmin*t,2)) - 
                                        (cmin*emin)/sqrt(1 + pow(-1 + bmax + bmin + bmax*t - bmin*t,2))))/2.);
  
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


p4est_geometry_t*
d4est_geometry_new_pizza_half
(
 p4est_connectivity_t * conn,
 double R0,
 double R1
)
{
  p4est_geometry_t *pizza_half_geometry = P4EST_ALLOC(p4est_geometry_t,1);

  double* radii = P4EST_ALLOC(double, 2);
  radii[0] = R0;
  radii[1] = R1;
  
  pizza_half_geometry->user = (void*)radii;
  pizza_half_geometry->destroy = d4est_geometry_pizza_half_destroy;
  pizza_half_geometry->X = d4est_geometry_pizza_half_X;

  return pizza_half_geometry;
}
