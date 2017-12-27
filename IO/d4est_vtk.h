#ifndef D4EST_VTK_H
#define D4EST_VTK_H

#include <d4est_operators.h>
#include <pXest.h>
#include <d4est_geometry.h>
/* SC_EXTERN_C_BEGIN; */



/** Opaque context type for writing VTK output with multiple function calls.
 */
typedef struct d4est_vtk_context d4est_vtk_context_t;
typedef enum {D4EST_VTK_DG_GRID, D4EST_VTK_CORNER_GRID} d4est_vtk_grid_type_t;
typedef enum {D4EST_VTK_ASCII, D4EST_VTK_BINARY} d4est_vtk_output_type_t;

typedef void
(*d4est_vtk_user_fcn_t)
(
 d4est_vtk_context_t*,
 void*
);

/* This file was automatically generated.  Do not edit! */
d4est_vtk_context_t *d4est_vtk_write_element_data(d4est_vtk_context_t *cont,int write_tree,int write_level,int write_rank,int wrap_rank,int write_deg,int num_element_fields,const char **names,double **values);
d4est_vtk_context_t *d4est_vtk_write_dg_data(d4est_vtk_context_t *cont,int num_point_scalars,const char **names,double **values,d4est_vtk_grid_type_t grid_type);
void d4est_vtk_save(p4est_t *p4est,d4est_operators_t *d4est_ops,const char *input_file,const char *input_section,const char *save_as_filename,const char **dg_field_names,double **dg_fields,const char **element_field_names,double **element_fields,d4est_vtk_grid_type_t grid_type,d4est_vtk_output_type_t output_type);
d4est_vtk_context_t *d4est_vtk_write_point_scalar(d4est_vtk_context_t *cont,const char *scalar_name,double *values,d4est_vtk_grid_type_t grid_type);
double *d4est_vtk_convert_nodal_to_vtk(p4est_t *p4est,d4est_vtk_context_t *cont,d4est_operators_t *d4est_ops,double *values,d4est_vtk_grid_type_t grid_type);
int d4est_vtk_write_footer(d4est_vtk_context_t *cont);
d4est_vtk_context_t *d4est_vtk_write_dg_cell_scalar(d4est_vtk_context_t *cont,const char *scalar_name,double *values);
d4est_vtk_context_t *d4est_vtk_write_dg_header(d4est_vtk_context_t *cont,d4est_operators_t *d4est_ops);
void d4est_vtk_context_destroy(d4est_vtk_context_t *context);
void d4est_vtk_context_set_continuous(d4est_vtk_context_t *cont,int continuous);
void d4est_vtk_context_set_scale(d4est_vtk_context_t *cont,double scale);
void d4est_vtk_context_set_deg_array(d4est_vtk_context_t *cont,int *deg_array);
void d4est_vtk_context_set_geom(d4est_vtk_context_t *cont,d4est_geometry_t *geom);
d4est_vtk_context_t *d4est_vtk_dg_context_new(p4est_t *p4est,d4est_operators_t *d4est_ops,const char *filename);

#endif /* !D4EST_VTK_H */
