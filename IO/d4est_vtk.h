#ifndef D4EST_VTK_H
#define D4EST_VTK_H

#include <d4est_operators.h>
#include <pXest.h>
#include <d4est_geometry.h>
/* SC_EXTERN_C_BEGIN; */



/** Opaque context type for writing VTK output with multiple function calls.
 */
typedef struct d4est_vtk_context d4est_vtk_context_t;

typedef void
(*d4est_vtk_user_fcn_t)
(
 d4est_vtk_context_t*,
 void*
);


/* This file was automatically generated.  Do not edit! */
void d4est_vtk_save(p4est_t *p4est,d4est_operators_t *d4est_ops,const char *input_file,const char *input_section,const char *print_prefix,const char **dg_field_names,double **dg_fields,const char **element_field_names,double **element_fields,int folder_number);

#endif /* !D4EST_VTK_H */
