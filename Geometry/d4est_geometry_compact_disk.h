#ifndef D4EST_GEOMETRY_COMPACT_DISK_H
#define D4EST_GEOMETRY_COMPACT_DISK_H 


#include <p4est_connectivity.h>
#include <p4est_geometry.h>

p4est_geometry_t *d4est_geometry_new_compact_disk(p4est_connectivity_t *conn, double R0,
                                                  double R1,double w, double Rinf);



#endif
