#ifndef D4EST_VTK_H
#define D4EST_VTK_H

#include <d4est_operators.h>
#include <pXest.h>
#include <d4est_geometry.h>

/** Opaque context type for writing VTK output with multiple function calls.
 */
typedef struct d4est_vtk_context d4est_vtk_context_t;

typedef void
(*d4est_vtk_user_fcn_t)
(
 d4est_vtk_context_t*,
 void*
);


typedef struct {

  p4est_t* p4est;
  d4est_operators_t* d4est_ops;

  int name_size;
  int max_num_fields;
  int local_nodes;
  int local_elements;

  double** nodal_dbl_fields;
  char** nodal_dbl_fields_names;
  int num_nodal_dbl_fields;
  int* nodal_dbl_did_we_alloc;
  
  double** cell_dbl_fields;
  char** cell_dbl_fields_names;
  int num_cell_dbl_fields;
  int* cell_dbl_did_we_alloc;
  
  int** cell_int_fields;
  char** cell_int_fields_names;
  int num_cell_int_fields;
  int* cell_int_did_we_alloc;
  
} d4est_vtk_helper_array_t;


/* This file was automatically generated.  Do not edit! */
void d4est_vtk_helper_array_destroy(d4est_vtk_helper_array_t *array);
void d4est_vtk_save_helper_array(d4est_vtk_helper_array_t *array,char *input_file);
int *d4est_vtk_helper_array_alloc_and_add_cell_int_field(d4est_vtk_helper_array_t *array,const char *prefix,int suffix_id);
double *d4est_vtk_helper_array_alloc_and_add_cell_dbl_field(d4est_vtk_helper_array_t *array,const char *prefix,int suffix_id);
double *d4est_vtk_helper_array_alloc_and_add_nodal_dbl_field(d4est_vtk_helper_array_t *array,const char *prefix,int suffix_id);
void d4est_vtk_helper_array_add_nodal_dbl_field(d4est_vtk_helper_array_t *array,const char *prefix,int suffix_id,double *field);
d4est_vtk_helper_array_t *d4est_vtk_helper_array_init(p4est_t *p4est,d4est_operators_t *d4est_ops,int local_nodes,int name_size,int max_num_fields);
void d4est_vtk_save(p4est_t *p4est,d4est_operators_t *d4est_ops,const char *input_file,const char *input_section,const char **dg_field_names,double **dg_fields,const char **element_float_field_names,double **element_float_fields,const char **element_int_field_names,int **element_int_fields,int sub_folder_number);
void d4est_vtk_save_aux(p4est_t *p4est,d4est_operators_t *d4est_ops,const char *input_file,const char *input_section,const char **dg_field_names,double **dg_fields,const char **element_float_field_names,double **element_float_fields,const char **element_int_field_names,int **element_int_fields,const char *folder,int sub_folder_number);

#endif
