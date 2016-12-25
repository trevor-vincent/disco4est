#include <util.h>
#include <d4est_geometry_compact_sphere.h>

typedef struct
{
  double R2, R1, R0, w, Rinf;
  double R2byR1, R1sqrbyR2, R1log;
  double R1byR0, R0sqrbyR1, R0log;
  double Clength, CdetJ;

  p4est_connectivity_t* conn;
  
} d4est_geometry_compact_sphere_attr_t;

static void
d4est_geometry_compact_sphere_destroy
(
 p8est_geometry_t* geom
)
{
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}

static void
d4est_geometry_octree_to_vertex (p8est_geometry_t * geom,
                               p4est_topidx_t which_tree,
                               const double abc[3], double xyz[3])
{
  d4est_geometry_compact_sphere_attr_t *compact_sphere_att = (d4est_geometry_compact_sphere_attr_t *) geom->user;
  p4est_connectivity_t *connectivity = (p4est_connectivity_t *) compact_sphere_att->conn;
  const p4est_topidx_t *tree_to_vertex = connectivity->tree_to_vertex;
  const double       *v = connectivity->vertices;
  double              eta_x, eta_y, eta_z = 0.;
  int                 j, k;
  p4est_topidx_t      vt[P4EST_CHILDREN];

  /* retrieve corners of the tree */
  for (k = 0; k < P4EST_CHILDREN; ++k) {
    vt[k] = tree_to_vertex[which_tree * P4EST_CHILDREN + k];
  }

  /* these are reference coordinates in [0, 1]**d */
  eta_x = abc[0];
  eta_y = abc[1];
  eta_z = abc[2];

  /* bi/trilinear transformation */
  for (j = 0; j < 3; ++j) {
    /* *INDENT-OFF* */
    xyz[j] =
           ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                  eta_x  * v[3 * vt[1] + j]) +
                                  eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                                  eta_x  * v[3 * vt[3] + j]))
#ifdef P4_TO_P8
            +     eta_z  * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[4] + j] +
                                                  eta_x  * v[3 * vt[5] + j]) +
                                  eta_y  * ((1. - eta_x) * v[3 * vt[6] + j] +
                                                  eta_x  * v[3 * vt[7] + j]))
#endif
           );
    /* *INDENT-ON* */
  }
}


static void
d4est_geometry_compact_sphere_X(p8est_geometry_t * geom,
                       p4est_topidx_t which_tree,
                       const double rst[3], double xyz[3])
{
  d4est_geometry_compact_sphere_attr_t* compact_sphere = geom->user;

  double              x, y, R, q;
  double              abc[3];

  /* double R0 = compact_sphere->R0; */
  double R1 = compact_sphere->R1;
  double R2 = compact_sphere->R2;
  double Rinf = compact_sphere->Rinf;
  double w = compact_sphere->w;

  /* DEBUG_PRINT_DBL(R1); */
  /* DEBUG_PRINT_DBL(R2); */
  /* DEBUG_PRINT_DBL(Rinf); */
  /* DEBUG_PRINT_DBL(w); */
  
  /* transform from the reference cube into vertex space */
  d4est_geometry_octree_to_vertex (geom, which_tree, rst, abc);

  /* assert that input points are in the expected range */
  P4EST_ASSERT (0 <= which_tree && which_tree < 13);

  if (which_tree < 6) {         /* outer shell */
    x = tan (abc[0] * M_PI_4);
    y = tan (abc[1] * M_PI_4);
    /* R = compact_sphere->R1sqrbyR2 * pow (compact_sphere->R2byR1, abc[2]); */
    R = compact_sphere->R1*(2. - abc[2]) + compact_sphere->R2*(abc[2] - 1.);
    q = R / sqrt (x * x + y * y + 1.);
  }
  else if (which_tree < 12) {   /* inner shell */
    double              p, tanx, tany;

    p = 2. - abc[2];
    tanx = tan (abc[0] * M_PI_4);
    tany = tan (abc[1] * M_PI_4);
    x = p * abc[0] + (1. - p) * tanx;
    y = p * abc[1] + (1. - p) * tany;
    /* R = compact_sphere->R0sqrbyR1 * pow (compact_sphere->R1byR0, abc[2]); */
    R = compact_sphere->R0*(2. - abc[2]) + compact_sphere->R1*(abc[2] - 1.);
    
    q = R / sqrt (1. + (1. - p) * (tanx * tanx + tany * tany) + 2. * p);
  }
  else {                        /* center cube */
    xyz[0] = abc[0] * compact_sphere->Clength;
    xyz[1] = abc[1] * compact_sphere->Clength;
    xyz[2] = abc[2] * compact_sphere->Clength;

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

  /* compactify */
  if (which_tree < 6){
    double xnc = xyz[0];
    double ync = xyz[1];
    double znc = xyz[2];
    double R = sqrt(xnc*xnc + ync*ync + znc*znc);
    double r = (R1 + (w - R1)*(R - R1)/(R2 - R1))/(1 - (R - R1)*(1 - w/Rinf)/(R2 - R1));
    xyz[0] = r*xnc/R;
    xyz[1] = r*ync/R;
    xyz[2] = r*znc/R;
  }
}

p8est_geometry_t*
d4est_geometry_new_compact_sphere
(
 p4est_connectivity_t* conn,
 double R2,
 double R1,
 double R0,
 double w,
 double Rinf
)
{
  p8est_geometry_t* compact_sphere_geom = P4EST_ALLOC(p8est_geometry_t, 1);

  d4est_geometry_compact_sphere_attr_t* compact_sphere_attrs = P4EST_ALLOC(d4est_geometry_compact_sphere_attr_t, 1);
  
  compact_sphere_attrs->R2 = R2;
  compact_sphere_attrs->R1 = R1;
  compact_sphere_attrs->R0 = R0;
  compact_sphere_attrs->w = w;
  compact_sphere_attrs->Rinf = Rinf;
  compact_sphere_attrs->conn = conn;
  
  /* variables useful for the outer shell */
  compact_sphere_attrs->R2byR1 = R2 / R1;
  compact_sphere_attrs->R1sqrbyR2 = R1 * R1 / R2;
  compact_sphere_attrs->R1log = log (R2 / R1);

  /* variables useful for the inner shell */
  compact_sphere_attrs->R1byR0 = R1 / R0;
  compact_sphere_attrs->R0sqrbyR1 = R0 * R0 / R1;
  compact_sphere_attrs->R0log = log (R1 / R0);

  /* variables useful for the center cube */
  compact_sphere_attrs->Clength = R0 / sqrt (3.);
  compact_sphere_attrs->CdetJ = pow (R0 / sqrt (3.), 3.);

  compact_sphere_geom->name = "d4est_geometry_compact_sphere";
  compact_sphere_geom->user = compact_sphere_attrs;
  compact_sphere_geom->X = d4est_geometry_compact_sphere_X;
  compact_sphere_geom->destroy = d4est_geometry_compact_sphere_destroy;
  
  return compact_sphere_geom;
}
