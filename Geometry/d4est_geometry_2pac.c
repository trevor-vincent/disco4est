#include <d4est_geometry_2pac.h>
#include <d4est_geometry_pizza_half.h>


typedef struct{

  double radius;
  double base_height_half;

} tupac_attr_t;


static
void d4est_geometry_2pac_aligned_linear_map
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
d4est_geometry_2pac_aligned_map_cube_to_slab(
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
  d4est_geometry_2pac_aligned_linear_map(xref, 0., 1., emin, emax, &xbar);
    /* map y from [0,1] to [-1,1] */
  d4est_geometry_2pac_aligned_linear_map(yref, 0., 1., -1., 1., &ybar);
  double xmin = (1. - cmin)*emin + emin*cmin/sqrt(1. + ybar*ybar);
  double xmax = (1. - cmax)*emax + emax*cmax/sqrt(1. + ybar*ybar);
  *x = xmin + (xmax - xmin)*(xbar - emin)/(emax - emin);
  *y = (*x)*ybar;
}


static void 
d4est_geometry_2pac_aligned_X(p4est_geometry_t * geom,
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
    /* left */
    d4est_geometry_2pac_aligned_map_cube_to_slab(xref, yref, 0, 1, R0/sqrt(2), R1, &x, &y);
    /* double x0 = x; */
    x -= R0/sqrt(2) - R1;
    x -= (R1 - R0/sqrt(2));
    /* printf("0,xref,yref,x,y = %f,%f,%f,%f\n",xref,yref,x,y); */
    /* printf("xbefore,xafter = %f,%f\n",x0,x); */
    /* printf("x,y = %f, %f\n", x,y); */
  }
  else if (which_tree == 1){
    /* right */
    d4est_geometry_2pac_aligned_map_cube_to_slab(xref, yref, 1, 1, R1, 2*R1, &x, &y);
    /* printf("1,xref,yref,x,y = %f,%f,%f,%f\n",xref,yref,x,y); */
    /* x -= R0/sqrt(2); */
    /* printf("x,y = %f, %f\n", x,y); */
  }
  else{
    SC_ABORT_NOT_REACHED();
  }

  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = 0.;
}

void 
d4est_geometry_2pac_aligned_dxdr_old(p4est_geometry_t * geom,
                               p4est_topidx_t which_tree,
                               const double rst[3], // \in [-1,1] such as GL or GLL points
                                     const p4est_qcoord_t q [(P4EST_DIM)],
                                     const p4est_qcoord_t dq,
                               double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)])
{
  double* radii = (double*)geom->user;
  
  double R0 = radii[0];
  double R1 = radii[1];
  
  double r = rst[0];
  double s = rst[1];
  /* double x,y; */

  double cmin, cmax, emin, emax;

  double amin = ((double)q[0])/(double)P4EST_ROOT_LEN;
  double amax = ((double)q[0] + (double)dq)/(double)P4EST_ROOT_LEN;
  double bmin = ((double)q[1])/(double)P4EST_ROOT_LEN;
  double bmax = ((double)q[1] + (double)dq)/(double)P4EST_ROOT_LEN;
  
  if (which_tree == 0){
    /* left */
    cmin = 0.;
    cmax = 1.;
    emin = R0/sqrt(2);
    emax = R1;
  }
  else if (which_tree == 1){
    cmin = 1.;
    cmax = 1.;
    emin = R1;
    emax = 2*R1;
  }
  else{
    SC_ABORT_NOT_REACHED();
  }

  dxyz_drst[0][0] = ((amax - amin)*(-1.*(-1. + cmax)*emax + (-1. + cmin)*emin + (cmax*emax)/sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2)) - (cmin*emin)/sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2))))/2.;

  dxyz_drst[0][1] = -((-bmax + bmin)*(2*cmin*emin - amin*(cmax*emax - cmin*emin)*(-1 + r) + amax*(cmax*emax - cmin*emin)*(1 + r))*(1 + bmin*(-1 + s) - bmax*(1 + s)))/
                    (2.*pow(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2),1.5));

  dxyz_drst[1][0] = ((amax - amin)*(-1 + bmax + bmin + bmax*s - bmin*s)*(-1.*(-1. + cmax)*emax + (-1. + cmin)*emin + (cmax*emax)/sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2)) - 
                                                                         (cmin*emin)/sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2))))/2.;

  dxyz_drst[1][1] = ((bmax - bmin)*(-(((2*cmin*emin - amin*(cmax*emax - cmin*emin)*(-1 + r) + amax*(cmax*emax - cmin*emin)*(1 + r))*pow(-1 + bmax + bmin + bmax*s - bmin*s,2))/
          (sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2))*(2. + 2*bmin*(-1 + s) + pow(bmin,2)*pow(-1 + s,2) + pow(bmax,2)*pow(1 + s,2) - 2*bmax*(1 + s + bmin*(-1 + pow(s,2)))))) + 
       2*(-1.*(-1. + cmin)*emin + (cmin*emin)/sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2)) + 
          ((amax + amin + amax*r - amin*r)*(-1.*(-1. + cmax)*emax + (-1. + cmin)*emin + (cmax*emax)/sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2)) - 
                                            (cmin*emin)/sqrt(1. + pow(-1 + bmax + bmin + bmax*s - bmin*s,2))))/2.)))/2.;
 
}


void 
d4est_geometry_2pac_aligned_dxdr(p4est_geometry_t * geom,
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

  double cmin, cmax, emin, emax;

  double amin = q0[0]; amin /= (double)(P4EST_ROOT_LEN);
  double amax = q0[0] + dq; amax /= (double)(P4EST_ROOT_LEN);
  double bmin = q0[1]; bmin /= (double)(P4EST_ROOT_LEN);
  double bmax = q0[1] + dq; bmax /= (double)(P4EST_ROOT_LEN);
  
  if (which_tree == 0){
    /* left */
    cmin = 0.;
    cmax = 1.;
    emin = R0/sqrt(2);
    emax = R1;
  }
  else if (which_tree == 1){
    cmin = 1.;
    cmax = 1.;
    emin = R1;
    emax = 2*R1;
  }
  else{
    SC_ABORT_NOT_REACHED();
  }

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
d4est_geometry_2pac_aligned_dxdr_face(
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
  if (which_tree == 0){
    /* left */
    cmin = 0.;
    cmax = 1.;
    emin = R0/sqrt(2);
    emax = R1;
  }
  else if (which_tree == 1){
    cmin = 1.;
    cmax = 1.;
    emin = R1;
    emax = 2*R1;
  }
  else{
    SC_ABORT_NOT_REACHED();
  }

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
d4est_geometry_2pac_aligned_X_traps(p4est_geometry_t * geom,
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
    /* left */
    d4est_geometry_2pac_aligned_map_cube_to_slab(xref, yref, 0, 0, R0/sqrt(2), R1, &x, &y);

    /* double x0 = x; */
    
    x -= R0/sqrt(2) - R1;
    x -= (R1 - R0/sqrt(2));
    /* printf("x,y = %f, %f\n", x,y); */
  }
  else if (which_tree == 1){
    /* right */
    d4est_geometry_2pac_aligned_map_cube_to_slab(xref, yref, 0, 0, R1, 2*R1, &x, &y);
    /* x -= R0/sqrt(2); */
    /* printf("x,y = %f, %f\n", x,y); */
  }
  else{
    SC_ABORT_NOT_REACHED();
  }


  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = 0.;
}

static void 
d4est_geometry_2pac_aligned_X_squares(p4est_geometry_t * geom,
                       p4est_topidx_t which_tree,
                       const double rst[3], double xyz[3])
{
  /* double* radii = (double*)geom->user; */
  
  /* double R0 = radii[0]; */
  /* double R1 = radii[1]; */
  
  double xref = rst[0];
  double yref = rst[1];
  double x,y;

  if (which_tree == 0){
    d4est_geometry_2pac_aligned_linear_map(yref, 0., 1., 0., 1.0, &y);
    d4est_geometry_2pac_aligned_linear_map(xref, 0., 1., -1., 0, &x);    
  }
  else if (which_tree == 1){
    d4est_geometry_2pac_aligned_linear_map(yref, 0., 1., 0., 1., &y);
    d4est_geometry_2pac_aligned_linear_map(xref, 0., 1., 0., 1., &x);    
  }
  else{
    SC_ABORT_NOT_REACHED();
  }

  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = 0.;
}

/* static void  */
/* d4est_geometry_2pac_aligned_X */
/* ( */
/*  p4est_geometry_t * geom, */
/*  p4est_topidx_t which_tree, */
/*  const double rst[3], */
/*  double xyz[3] */
/* ) */
/* { */
/*   double* radii = (double*)geom->user; */
  
/*   double R0 = radii[0]; */
/*   double R1 = radii[1]; */
  
/*   double xref = rst[0]; */
/*   double yref = rst[1]; */
/*   double x,y; */

/*   if (which_tree == 0){ */
/*     /\* left *\/ */
/*     d4est_geometry_2pac_aligned_map_cube_to_slab(xref, yref, 1., 0., -R1, -R0/sqrt(2), &x, &y); */
/*     y *= -1; */
/*     x += R0/sqrt(2); */
/*     /\* printf("x,y = %f, %f\n", x,y); *\/ */
/*   } */
/*   else if (which_tree == 1){ */
/*     /\* right *\/ */
/*     d4est_geometry_2pac_aligned_map_cube_to_slab(xref, yref, 0, 1, R0/sqrt(2), R1, &x, &y); */
/*     x -= R0/sqrt(2); */
/*     /\* printf("x,y = %f, %f\n", x,y); *\/ */
/*   } */
/*   else{ */
/*     SC_ABORT_NOT_REACHED(); */
/*   } */
/*   xyz[0] = x; */
/*   xyz[1] = y; */
/*   xyz[2] = 0.; */
/* } */



static void
d4est_geometry_2pac_aligned_destroy
(
 p4est_geometry_t* geom
)
{
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}


p4est_geometry_t*
d4est_geometry_new_2pac_aligned
(
 p4est_connectivity_t * conn,
 double R0,
 double R1
)
{
  p4est_geometry_t* tupac_aligned_geometry = P4EST_ALLOC(p4est_geometry_t,1);

  double* radii = P4EST_ALLOC(double, 2);
  radii[0] = R0;
  radii[1] = R1;
    
  tupac_aligned_geometry->user = (void*)radii;
  tupac_aligned_geometry->destroy = d4est_geometry_2pac_aligned_destroy;
  tupac_aligned_geometry->X = d4est_geometry_2pac_aligned_X;

  return tupac_aligned_geometry;
}

p4est_geometry_t*
d4est_geometry_new_2pac_aligned_traps
(
 p4est_connectivity_t * conn,
 double R0,
 double R1
)
{
  p4est_geometry_t* tupac_aligned_geometry = P4EST_ALLOC(p4est_geometry_t,1);

  double* radii = P4EST_ALLOC(double, 2);
  radii[0] = R0;
  radii[1] = R1;
    
  tupac_aligned_geometry->user = (void*)radii;
  tupac_aligned_geometry->destroy = d4est_geometry_2pac_aligned_destroy;
  tupac_aligned_geometry->X = d4est_geometry_2pac_aligned_X_traps;

  return tupac_aligned_geometry;
}

p4est_geometry_t*
d4est_geometry_new_2pac_aligned_squares
(
 p4est_connectivity_t * conn,
 double R0,
 double R1
)
{
  p4est_geometry_t* tupac_aligned_geometry = P4EST_ALLOC(p4est_geometry_t,1);

  double* radii = P4EST_ALLOC(double, 2);
  radii[0] = R0;
  radii[1] = R1;
    
  tupac_aligned_geometry->user = (void*)radii;
  tupac_aligned_geometry->destroy = d4est_geometry_2pac_aligned_destroy;
  tupac_aligned_geometry->X = d4est_geometry_2pac_aligned_X_squares;

  return tupac_aligned_geometry;
}

//

static
void d4est_geometry_2pac_nonaligned_linear_map
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
d4est_geometry_2pac_nonaligned_map_cube_to_slab(
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
  d4est_geometry_2pac_nonaligned_linear_map(xref, 0., 1., emin, emax, &xbar);
    /* map y from [0,1] to [-1,1] */
  d4est_geometry_2pac_nonaligned_linear_map(yref, 0., 1., -1., 1., &ybar);
  double xmin = (1. - cmin)*emin + emin*cmin/sqrt(1. + ybar*ybar);
  double xmax = (1. - cmax)*emax + emax*cmax/sqrt(1. + ybar*ybar);
  *x = xmin + (xmax - xmin)*(xbar - emin)/(emax - emin);
  *y = (*x)*ybar;
}


static void 
d4est_geometry_2pac_nonaligned_X(p4est_geometry_t * geom,
                       p4est_topidx_t which_tree,
                       const double rst[3], double xyz[3])
{
  double R0 = *((double*)geom->user);
  
  double xref = rst[0];
  double yref = rst[1];
  double x,y;
  

  if (which_tree == 0){
    /* left */
    d4est_geometry_2pac_nonaligned_map_cube_to_slab(xref, yref, 1., 0., -R0, 0., &x, &y);
    y *= -1;
  }
  else if (which_tree == 1){
    /* right */
    d4est_geometry_2pac_nonaligned_map_cube_to_slab(yref, xref, 0, 1, 0, R0, &x, &y);
    y *= -1.;
  }
  else{
    SC_ABORT_NOT_REACHED();
  }

  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = 0.;
}

static void
d4est_geometry_2pac_nonaligned_destroy
(
 p4est_geometry_t* geom
)
{
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}


p4est_geometry_t*
d4est_geometry_new_2pac_nonaligned
(
 p4est_connectivity_t * conn,
 double R0
)
{
  p4est_geometry_t * tupac_nonaligned_geometry = P4EST_ALLOC(p4est_geometry_t,1);

  double* R0p = P4EST_ALLOC(double, 1);
  *R0p = R0;
                
  tupac_nonaligned_geometry->user = (void*)R0p;
  tupac_nonaligned_geometry->destroy = d4est_geometry_2pac_nonaligned_destroy;
  tupac_nonaligned_geometry->X = d4est_geometry_2pac_nonaligned_X;

  return tupac_nonaligned_geometry;
}
