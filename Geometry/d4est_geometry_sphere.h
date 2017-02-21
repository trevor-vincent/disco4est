#ifndef D4EST_GEOMETRY_SPHERE_H
#define D4EST_GEOMETRY_SPHERE_H

#include <pXest.h>
#include <p8est_geometry.h>


typedef struct
{
  double R2, R1, R0;
  double R2byR1, R1sqrbyR2, R1log;
  double R1byR0, R0sqrbyR1, R0log;
  double Clength, CdetJ;

  p4est_connectivity_t* conn;
  
} d4est_geometry_sphere_attr_t;

void d4est_geometry_compactified_sphere_new(const char *input_file,d4est_geometry_t *d4est_geom);
void d4est_geometry_sphere_new(const char *input_file,d4est_geometry_t *d4est_geom);

p8est_geometry_t*
d4est_geometry_compactified_sphere_from_param
(
 double R0,
 double R1,
 double R2,
 p8est_connectivity_t* conn
);

#endif
