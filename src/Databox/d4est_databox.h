#ifndef D4EST_DATABOX_H
#define D4EST_DATABOX_H

#include <d4est_dictionary.h>
#include <d4est_field.h>
#include <d4est_element_data.h>
/* #include <d4est_ghost_data.h> */

typedef double
(*d4est_databox_init_fcn_t)
(
 double,
 double,
#if (P4EST_DIM)==3
 double,
#endif 
 void*
);


typedef struct {
  
  int local_nodes;       /**< size of local nodal space */
  
} d4est_databox_sizes_t;


typedef struct {
  int mpirank; /**< mpirank or id of subdomain */
  d4est_dictionary_t fields; /**< dictionary for fields defined on local grid */
  int field_sizes [D4EST_FIELD_TYPES];
 
} d4est_databox_t;

void d4est_databox_destroy(d4est_databox_t *s);
d4est_databox_t *d4est_databox_init(int mpirank,int *loc_sizes);
void d4est_databox_get_vector_field(d4est_databox_t *s,const char *prefix,double *vfield[(P4EST_DIM)]);
void d4est_databox_delete_field(d4est_databox_t *s,const char *name);
int d4est_databox_add_field(d4est_databox_t *s,const char *name,d4est_field_type_t field_type);
int d4est_databox_get_field_type(d4est_databox_t *s,const char *name,d4est_field_type_t *type);
double *d4est_databox_get_field(d4est_databox_t *s,const char *name);
void d4est_databox_init_field(d4est_databox_t *s,const char *name,const char *x,const char *y,const char *z,d4est_databox_init_fcn_t init_fcn,void *arg);

#endif
