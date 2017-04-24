#ifndef D4EST_GEOMETRY_5TREEDISK_H
#define D4EST_GEOMETRY_5TREEDISK_H

#include <d4est_geometry.h>
void d4est_geometry_disk_outer_wedge_new(int mpirank,const char *input_file,const char *input_section,const char *printf_prefix,d4est_geometry_t *d4est_geom);
void d4est_geometry_5treedisk_new(int mpirank,const char *input_file,const char *input_section,const char *printf_prefix,d4est_geometry_t *d4est_geom);
#endif
