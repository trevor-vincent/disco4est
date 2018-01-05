#include <pXest.h>
#include <d4est_geometry.h>
#include <d4est_geometry_general_wedge.h>

void
d4est_geometry_general_wedge_X
(
 d4est_geometry_t * geom,
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [3],
 p4est_qcoord_t dq,
 const double coords[3],
 coords_type_t coords_type,
 double xyz[3],
 double curvature_at_zmin,
 double curvature_at_zmax,
 double zmin,
 double zmax,
 int compactified,
 d4est_wedge_part_t which_part
)
{  
  if (compactified == 1){
    D4EST_ABORT("Compactification not supported yet");
  }

  double tcoords [3];
  /* transform coords of element to [0,1]^3*/
  d4est_geometry_get_tree_coords_in_range_0_to_1(q0, dq, coords, coords_type, tcoords);

  double a = tcoords[0];
  double b = tcoords[1];
  double c = tcoords[2];

  if (which_part == LEFT_WEDGE){
    a -= 1;
    b -= 1;
  }
  if (which_part == FULL_WEDGE){
    a = 2*a - 1;
    b = 2*b - 1;
  }
  
  double R = c*(zmax - zmin) + zmin;
  double x = tan(.25*M_PI*a);
  double y = tan(.25*M_PI*b);
  double p = 1./sqrt(1 + x*x + y*y);
  double fmin = zmin*(1+curvature_at_zmin*(p-1));
  double fmax = zmax*(1+curvature_at_zmax*(p-1));
  double q = fmin + (fmax - fmin)*(R-zmin)/(zmax - zmin);

  xyz[0] = +q * x;
  xyz[1] = +q * y;
  xyz[2] = +q; 
}


void
d4est_geometry_general_wedge_noncompactified_DX
(
 d4est_geometry_t* d4est_geom,
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 const double rst[(P4EST_DIM)], /* [-1,1]^3 */
 double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)],
 double curvature_at_zmin,
 double curvature_at_zmax,
 double zmin,
 double zmax,
 int compactified,
 d4est_wedge_part_t which_part
)
{  
  if (compactified == 1){
    D4EST_ABORT("Compactification not supported yet");
  }
  
  double r = rst[0];
  double s = rst[1];
  double t = rst[2];

  /* topological coordinates of element corners */
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;
  double cmin = q0[2];
  double cmax = q0[2] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;
  cmin /= (double)P4EST_ROOT_LEN;
  cmax /= (double)P4EST_ROOT_LEN;

  /* transform element corners to [-1,1]^2 x [0,1] topological space used in the cubed-sphere mapping*/
  if (which_part == LEFT_WEDGE){
    amin -= 1;
    amax -= 1;
    bmin -= 1;
    bmax -= 1;
  }
  if (which_part == FULL_WEDGE){
    amin = 2*amin - 1;
    amax = 2*amax - 1;
    bmin = 2*bmin - 1;
    bmax = 2*bmax - 1;
  }

  dxyz_drst[0][0] = ((amax - amin)*M_PI*pow(d4est_util_secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*(((curvature_at_zmin*zmin*(-2 + cmax + cmin + cmax*t - cmin*t) - curvature_at_zmax*zmax*(cmax + cmin + cmax*t - cmin*t))*pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))/pow(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5) + 2*(zmin*(1 + curvature_at_zmin*(-1 + 1/sqrt(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)))) + ((cmax + cmin + cmax*t - cmin*t)*(zmin*(-1 + curvature_at_zmin - curvature_at_zmin/sqrt(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))) + zmax*(1 + curvature_at_zmax*(-1 + 1/sqrt(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))))))/2.)))/16.;

  dxyz_drst[0][1] = ((bmax - bmin)*M_PI*(curvature_at_zmin*zmin*(-2 + cmax + cmin + cmax*t - cmin*t) - curvature_at_zmax*zmax*(cmax + cmin + cmax*t - cmin*t))*pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(16.*pow(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));

  dxyz_drst[0][2] = ((cmax - cmin)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*(zmin*(-1 + curvature_at_zmin - curvature_at_zmin/sqrt(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))) + zmax*(1 + curvature_at_zmax*(-1 + 1/sqrt(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))))))/2.;

  dxyz_drst[1][0] = ((amax - amin)*M_PI*(curvature_at_zmin*zmin*(-2 + cmax + cmin + cmax*t - cmin*t) - curvature_at_zmax*zmax*(cmax + cmin + cmax*t - cmin*t))*pow(d4est_util_secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(16.*pow(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));

  dxyz_drst[1][1] = ((bmax - bmin)*M_PI*pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*(2*(zmin*(1 + curvature_at_zmin*(-1 + 1/sqrt(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)))) + ((cmax + cmin + cmax*t - cmin*t)*(zmin*(-1 + curvature_at_zmin - curvature_at_zmin/sqrt(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))) + zmax*(1 + curvature_at_zmax*(-1 + 1/sqrt(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))))))/2.) + ((curvature_at_zmin*zmin*(-2 + cmax + cmin + cmax*t - cmin*t) - curvature_at_zmax*zmax*(cmax + cmin + cmax*t - cmin*t))*pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/pow(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5)))/16.;

  dxyz_drst[1][2] = ((cmax - cmin)*(zmin*(-1 + curvature_at_zmin - curvature_at_zmin/sqrt(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))) + zmax*(1 + curvature_at_zmax*(-1 + 1/sqrt(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)))))*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/2.;

  dxyz_drst[2][0] = ((amax - amin)*M_PI*(curvature_at_zmin*zmin*(-2 + cmax + cmin + cmax*t - cmin*t) - curvature_at_zmax*zmax*(cmax + cmin + cmax*t - cmin*t))*pow(d4est_util_secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(16.*pow(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));

  dxyz_drst[2][1] = ((bmax - bmin)*M_PI*(curvature_at_zmin*zmin*(-2 + cmax + cmin + cmax*t - cmin*t) - curvature_at_zmax*zmax*(cmax + cmin + cmax*t - cmin*t))*pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(16.*pow(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));

  dxyz_drst[2][2] = ((cmax - cmin)*(zmin*(-1 + curvature_at_zmin - curvature_at_zmin/sqrt(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))) + zmax*(1 + curvature_at_zmax*(-1 + 1/sqrt(pow(d4est_util_secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))))))/2.;

}
