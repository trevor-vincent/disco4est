#ifndef D4EST_GEOMETRY_GENERAL_CUBED_SPHERE_H
#define D4EST_GEOMETRY_GENERAL_CUBED_SPHERE_H 

typedef enum {LEFT_WEDGE, FULL_WEDGE, RIGHT_WEDGE} d4est_wedge_part_t;

/* This file was automatically generated.  Do not edit! */
void d4est_geometry_general_wedge_3D_DX(d4est_geometry_t *d4est_geom,p4est_topidx_t which_tree,p4est_qcoord_t q0[(P4EST_DIM)],p4est_qcoord_t dq,const double rst[(P4EST_DIM)],double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)],double curvature_at_zmin,double curvature_at_zmax,double zmin,double zmax,int compactified,d4est_wedge_part_t which_part);
void d4est_geometry_general_wedge_3D_X(d4est_geometry_t *geom,p4est_topidx_t which_tree,p4est_qcoord_t q0[3],p4est_qcoord_t dq,const double coords[3],coords_type_t coords_type,double xyz[3],double curvature_at_zmin,double curvature_at_zmax,double zmin,double zmax,int compactified,d4est_wedge_part_t which_part);


#endif
