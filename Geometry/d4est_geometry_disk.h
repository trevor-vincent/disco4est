#ifndef D4EST_GEOMETRY_5TREEDISK_H
#define D4EST_GEOMETRY_5TREEDISK_H

typedef struct {

  //R0 and R1 for 5 tree disk, inner square and wedge respectively
  double R0;
  double R1; 
  double R2; //added radii for 9 tree disk
  int compactify_outer_wedge;
  const char* input_section;
  
} d4est_geometry_disk_attr_t;

#include <d4est_geometry.h>
void d4est_geometry_disk_outer_wedge_new(int mpirank,const char *input_file,const char *input_section,const char *printf_prefix,d4est_geometry_t *d4est_geom);
void d4est_geometry_5treedisk_new(int mpirank,const char *input_file,const char *input_section,const char *printf_prefix,d4est_geometry_t *d4est_geom);
#endif
