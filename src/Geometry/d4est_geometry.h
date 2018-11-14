#ifndef D4EST_GEOMETRY_H
#define D4EST_GEOMETRY_H

#include <p8est_connectivity.h>
#include <pXest.h>
#include <d4est_operators.h>
#include <zlog.h>

typedef struct d4est_geometry d4est_geometry_t;

typedef enum {COMPUTE_NORMAL_USING_JACOBIAN, COMPUTE_NORMAL_USING_CROSS_PRODUCT} normal_compute_method_t;
typedef enum {GEOM_COMPUTE_NUMERICAL, GEOM_COMPUTE_ANALYTIC, GEOM_COMPUTE_NOT_SET} geometric_quantity_compute_method_t;

//         +----------------------+
//         |			  |
//         |			  |	     b
//         |			  |	     |
//         |			  |	     |
//         |	   face f	  |	     |
//         |	       	       	  |  	     |
//         |			  |	     /---------	a
//         |			  |	    /
//         |			  |	   /
//         |			  |	  c
//         |----------------------+

typedef struct {
  int a;      /* "x" coord on this face (z-ordering) */
  int b;      /* "y" coord on this face (z-ordering) */
  int c;      /* "normal" coord on this face (z-ordering) */
  double sgn; /* is the normal in the - or +  c-direction */
} d4est_geometry_face_info_t;


typedef enum {GEOM_CUBED_SPHERE_13TREE,
              GEOM_CUBED_SPHERE_7TREE,
              GEOM_CUBED_SPHERE_OUTER_WEDGE,
              GEOM_CUBED_SPHERE_INNER_WEDGE,
              GEOM_CUBED_SPHERE_INNEROUTER_WEDGE,
              GEOM_CUBED_SPHERE_WITH_CUBE_HOLE,
              GEOM_CUBED_SPHERE_WITH_SPHERE_HOLE,
              GEOM_DISK_5TREE,
              GEOM_DISK_OUTER_WEDGE,
              GEOM_HOLE_IN_A_BOX,
              GEOM_BRICK,
              GEOM_NONE} d4est_geometry_type_t;

/**
* There are three types of coords used to describe a point in an element.
* The first are the integration coords for gaussian quadrature. These are
* specified on the reference cube [-1,1]^DIM and represent a local coordinate
* system for the element. These are called COORDS_INTEG_RST. A point in an
* element can also have coordinates in the tree it is in, here d4est commonly
* used two different coords. The first is COORDS_P4EST_INT which is the built-in
* coords used by p4est, which represent points in the tree in the cube [0, P4EST_ROOT_LEN]^DIM
* where P4EST_ROOT_LEN is some big integer. These coords are stored as integers.
* The second is COORDS_TREE_UNITCUBE which are stored as dbls in the range [0,1]^DIM
*
*/

typedef enum {COORDS_INTEG_RST, COORDS_P4EST_INT, COORDS_TREE_UNITCUBE} coords_type_t;

typedef void        (*d4est_geometry_VEC_t) (d4est_geometry_t*,
                                           p4est_topidx_t ,
                                           p4est_qcoord_t [(P4EST_DIM)],
                                           p4est_qcoord_t,
                                           const double [(P4EST_DIM)],
                                           coords_type_t,
                                           double [(P4EST_DIM)]
                                          );

typedef void        (*d4est_geometry_MAT_t) (d4est_geometry_t*,
                                            p4est_topidx_t ,
                                            p4est_qcoord_t [(P4EST_DIM)],
                                            p4est_qcoord_t,
                                            const double [(P4EST_DIM)],
                                            double [(P4EST_DIM)][(P4EST_DIM)]
                                            );

typedef void        (*d4est_geometry_3DMAT_t) (d4est_geometry_t*,
                                            p4est_topidx_t ,
                                            p4est_qcoord_t [(P4EST_DIM)],
                                            p4est_qcoord_t,
                                            const double [(P4EST_DIM)],
                                            double [(P4EST_DIM)][(P4EST_DIM)][(P4EST_DIM)]
                                            );


typedef void        (*d4est_geometry_SCA_t) (d4est_geometry_t*,
                                            p4est_topidx_t ,
                                            p4est_qcoord_t [(P4EST_DIM)],
                                            p4est_qcoord_t,
                                            const double [(P4EST_DIM)],
                                            double*
                                            );


typedef void        (*d4est_geometry_destroy_t) (d4est_geometry_t * geom);

struct d4est_geometry {

  p4est_connectivity_t* p4est_conn;
  d4est_geometry_type_t geom_type;
  
  /* Mapping from [-1,1]^3 to grid coordinates*/
  /* geometric_quantity_compute_method_t X_compute_method; /\* only analytic atm *\/ */
  geometric_quantity_compute_method_t DX_compute_method; /* analytic and numerical available */
  geometric_quantity_compute_method_t JAC_compute_method; /* analytic and numerical available usually */
  
  /* Analytic derivatives of the mapping */
  d4est_geometry_SCA_t JAC;
  d4est_geometry_VEC_t X;
  d4est_geometry_MAT_t DX;
  d4est_geometry_3DMAT_t D2X;

  int(*get_number_of_regions)(d4est_geometry_t*);
  int(*get_region)(d4est_geometry_t*, p4est_qcoord_t [(P4EST_DIM)], p4est_qcoord_t, int);
  
  d4est_geometry_destroy_t destroy;
  void* user;
  
};
/* This file was automatically generated.  Do not edit! */
double d4est_geometry_compute_bounds(double *xyz[(P4EST_DIM)],int deg,double xi[(P4EST_DIM)],double xf[(P4EST_DIM)]);
double d4est_geometry_compute_lebesgue_measure(d4est_operators_t *d4est_ops,int deg_GL,double *jac_GL);
int d4est_geometry_is_face_on_boundary(p4est_t *p4est,p4est_quadrant_t *q,int which_tree,int face);
int d4est_geometry_does_element_touch_boundary(p4est_t *p4est,p4est_quadrant_t *q,int which_tree);
void d4est_geometry_get_tree_coords_in_range_0_to_1(p4est_qcoord_t q0[3],p4est_qcoord_t dq,const double coords[3],coords_type_t coords_type,double tcoords[3]);
void d4est_geometry_compute_xyz(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_rst_t rst_points,int which_tree,int deg,p4est_qcoord_t q[(P4EST_DIM)],p4est_qcoord_t dq,double *xyz[(P4EST_DIM)]);
void d4est_geometry_compute_drst_dxyz_times_jacobian(double *dxyz_drst[(P4EST_DIM)][(P4EST_DIM)],double *drst_dxyz_times_jac[(P4EST_DIM)][(P4EST_DIM)],int volume_nodes);
void d4est_geometry_compute_drst_dxyz(double *dxyz_drst[(P4EST_DIM)][(P4EST_DIM)],double *jac,double *drst_dxyz[(P4EST_DIM)][(P4EST_DIM)],int volume_nodes);
void d4est_geometry_compute_jacobian(double *dxyz_drst[(P4EST_DIM)][(P4EST_DIM)],double *jac,int volume_nodes);
void d4est_geometry_compute_xyz_face_analytic(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_rst_t rst_points,p4est_qcoord_t q0[(P4EST_DIM)],p4est_qcoord_t dq,int which_tree,int face,int deg,double *xyz[(P4EST_DIM)]);
void d4est_geometry_get_face_info(int f,d4est_geometry_face_info_t *face_info);
void d4est_geometry_compute_dxyz_drst_face_analytic(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_rst_t rst_points,p4est_qcoord_t q0[(P4EST_DIM)],p4est_qcoord_t dq,int which_tree,int face,int deg,double *dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]);
void d4est_geometry_compute_dxyz_drst_numerically(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_rst_t lobatto_rst_points,p4est_topidx_t which_tree,p4est_qcoord_t q0[(P4EST_DIM)],p4est_qcoord_t dq,int deg,double *dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]);
void d4est_geometry_compute_dxyz_drst(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_rst_t rst_points,p4est_topidx_t which_tree,p4est_qcoord_t q0[(P4EST_DIM)],p4est_qcoord_t dq,int deg,double *dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]);
void d4est_geometry_compute_dxyz_drst_analytic(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_rst_t rst_points,p4est_topidx_t which_tree,p4est_qcoord_t q0[(P4EST_DIM)],p4est_qcoord_t dq,int deg,double *dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]);
void d4est_geometry_quadtree_to_vertex(p4est_connectivity_t *connectivity,p4est_topidx_t which_tree,const double abc[3],double xyz[3]);
void d4est_geometry_octree_to_vertex(p8est_connectivity_t *connectivity,p4est_topidx_t which_tree,const double abc[3],double xyz[3]);
void d4est_geometry_destroy(d4est_geometry_t *d4est_geom);
d4est_geometry_t *d4est_geometry_new(int mpirank,const char *input_file,const char *input_section,zlog_category_t *c_default);

#endif

