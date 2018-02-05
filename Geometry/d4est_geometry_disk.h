#ifndef D4EST_GEOMETRY_5TREEDISK_H
#define D4EST_GEOMETRY_5TREEDISK_H

#include <zlog.h>


typedef struct {

  //R0 and R1 for 5 tree disk, inner square and wedge respectively
  double R0;
  double R1;
  double R2; //added radii for 9 tree disk
  int compactify_outer_wedge;
  const char* input_section;
  
} d4est_geometry_disk_attr_t;

#include <d4est_geometry.h>
void d4est_geometry_disk_outer_wedge_new(int mpirank,const char *input_file,const char *input_section,zlog_category_t *c_default,d4est_geometry_t *d4est_geom);
void d4est_geometry_disk_outer_wedge_new_aux(d4est_geometry_t *d4est_geom,d4est_geometry_disk_attr_t *disk_attrs);
void d4est_geometry_disk_outer_wedge_sj_div_jac_analytic(d4est_geometry_t *d4est_geom,p4est_topidx_t which_tree,p4est_qcoord_t q0[2],p4est_qcoord_t dq,const double rst[2],double *sj_div_jac);
void d4est_geometry_disk_outer_wedge_jac_analytic(d4est_geometry_t *d4est_geom,p4est_topidx_t which_tree,p4est_qcoord_t q0[2],p4est_qcoord_t dq,const double rst[2],double *j);
void d4est_geometry_disk_outer_wedge_sj_analytic(d4est_geometry_t *d4est_geom,p4est_topidx_t which_tree,p4est_qcoord_t q0[2],p4est_qcoord_t dq,const double rst[2],double *sj);
void d4est_geometry_5treedisk_new(int mpirank,const char *input_file,const char *input_section,zlog_category_t *c_default,d4est_geometry_t *d4est_geom);
#endif




