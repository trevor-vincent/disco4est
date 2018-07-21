#include <d4est_databox.h>
#include <d4est_element_data.h>
#include <d4est_field.h>
#include <d4est_util.h>
#define _GNU_SOURCE
#include <stdio.h>

static
int d4est_databox_free_fields_callback(const char *key, void *val, void *arg)
{
  d4est_field_t* field = (d4est_field_t*)val;
  field->field_type = NO_TYPE;
  D4EST_FREE(field->field_data);
  D4EST_FREE(field);
  return 1;
}

static
int d4est_databox_realloc_fields_callback(const char *key, void *val, void *arg)
{
  d4est_field_t* field = val;
  int* sizes = arg;
  D4EST_FIELD_CHECK_TYPE(field->field_type);
  int size = sizes[field->field_type];
  field->field_data = D4EST_REALLOC(field->field_data, double, size);
  return 1;
}

void
d4est_databox_init_field
(
 d4est_databox_t *s,
 const char *name,
 const char* x,
 const char* y,
 const char* z,
 d4est_databox_init_fcn_t init_fcn,
 void *arg
)
{
  double* restrict field = d4est_databox_get_field(s, name); 
  double* restrict xf = d4est_databox_get_field(s, x);
  double* restrict yf = d4est_databox_get_field(s, y);
#if (P4EST_DIM)==3
  double* restrict zf = d4est_databox_get_field(s, z);
#endif

  if (field == NULL ||
#if (P4EST_DIM)==3
      zf == NULL ||
#endif
      yf == NULL ||
      xf == NULL
     ){D4EST_ABORT("Init field error, can't find field, x, y, or z arrays");}

  d4est_field_type_t typef, typex, typey, typez;
  int rval = d4est_databox_get_field_type(s, name, &typef); D4EST_ASSERT(rval == 1);
  rval = d4est_databox_get_field_type(s, x, &typex); D4EST_ASSERT(rval == 1);
  rval = d4est_databox_get_field_type(s, y, &typey); D4EST_ASSERT(rval == 1);
  
#if (P4EST_DIM)==3
  rval = d4est_databox_get_field_type(s, z, &typez); D4EST_ASSERT(rval == 1);
  D4EST_ASSERT(typey == typez);
#endif
  D4EST_ASSERT(typef == typex && typex == typey);

  D4EST_FIELD_CHECK_TYPE(typef);
  int size = s->field_sizes[typef];
  for (int i = 0; i < size; i++){
    field[i] = init_fcn
               (
                xf[i],
                yf[i],
#if (P4EST_DIM)==3
                zf[i],
#endif
                arg
               );
  }
}

int
d4est_databox_add_field
(
 d4est_databox_t *s,
 const char *name,
 d4est_field_type_t field_type
)
{
  int size = s->field_sizes[field_type];
  if(d4est_dictionary_get_value_ptr(&s->fields,name))
    return 1;

  d4est_field_t* field = D4EST_ALLOC(d4est_field_t, 1);
  double* field_data = D4EST_ALLOC(double, size);

  field->field_data = field_data;
  field->field_type = field_type;
  
  int rval = d4est_dictionary_insert_ptr(&s->fields, name, field);

  D4EST_ASSERT(rval != 1);
  if(rval == 0){
    D4EST_FREE(field_data);
    D4EST_FREE(field);
    D4EST_ABORT("Out of memory, can't store field in dictionary\n");
  }  
  return rval;
}

void
d4est_databox_delete_field
(
 d4est_databox_t *s,
 const char *name
)
{
  d4est_field_t* field = d4est_dictionary_get_value_ptr(&s->fields, name);
  if (field == NULL){
    D4EST_ABORT("Could not find field to delete");
    return;
  }
  field->field_type = NO_TYPE;
  D4EST_FREE(field->field_data);
  D4EST_FREE(field);
  if(d4est_dictionary_delete(&s->fields, name) == 0){
    D4EST_ABORT("Could not delete field from dictionary\n");
  }
}

double*
d4est_databox_get_field
(
 d4est_databox_t* s,
 const char *name
)
{
  d4est_field_t* field = d4est_dictionary_get_value_ptr(&s->fields, name);
  if (field == NULL)
    return NULL;
  else
    return field->field_data;
}

void
d4est_databox_get_vector_field
(
 d4est_databox_t* s,
 const char *prefix,
 double* vfield [(P4EST_DIM)]
)
{
  char* vx = NULL;
  char* vy = NULL;
  asprintf(&vx, "%sx", prefix);
  asprintf(&vy, "%sy", prefix);
  vfield[0] = d4est_databox_get_field(s, vx);
  vfield[1] = d4est_databox_get_field(s, vy);
  free(vx);
  free(vy);
  if ((P4EST_DIM) == 3){
    char* vz = NULL;
    asprintf(&vz, "%sz", prefix);
    vfield[2] = d4est_databox_get_field(s, vz);
    free(vz);
  }
}


int
d4est_databox_get_field_type
(
 d4est_databox_t* s,
 const char *name,
 d4est_field_type_t* type
)
{
  d4est_field_t* field = d4est_dictionary_get_value_ptr(&s->fields, name);
  if (field == NULL){
    type = NULL;
    return 0;
  }
  else {
    *type = field->field_type;
    return 1;
  }        
}



/* double* d4est_databox_get_field_on_ghost */
/* ( */
/*  d4est_element_data_t* ed, */
/*  int ghost_name_id, */
/*  const char* name, */
/*  d4est_ghost_data_t* dgd */
/* ){ */
/*   D4EST_ASSERT(dgd != NULL); */
/*   if (ghost_name_id == -1){ */
/*     D4EST_ASSERT(name != NULL); */
/*     for (int i = 0; i < dgd->num_names; i++){ */
/*       if(d4est_util_compare_strings(name, dgd->transfer_names[i])){ */
/*         ghost_name_id = i; */
/*         break; */
/*       } */
/*     } */
/*   } */
/*   D4EST_ASSERT(ghost_name_id >= 0 && ghost_name_id < dgd->num_names); */
/*   if (name != NULL){ */
/*     D4EST_ASSERT(d4est_util_compare_strings(dgd->transfer_names[ghost_name_id], name)); */
/*   } */
/*   int stride = dgd->receive_strides[ed->id][ghost_name_id]; */
/*   return &dgd->receive_data[stride]; */
/* } */


/* double* d4est_databox_get_field_on_element_or_ghost */
/* ( */
/*  d4est_element_data_t* ed, */
/*  const char* name, */
/*  d4est_field_type_t type, */
/*  double* field_on_local_mesh_if_available, */
/*  d4est_databox_data_t* dmd, */
/*  d4est_ghost_data_t* dgd, */
/*  int ghost_name_id //for faster lookup */
/* ) */
/* { */
/*   D4EST_ASSERT(dmd != NULL && dgd != NULL); */
/*   if (ed->mpi_rank == dmd->mpi_rank){ */
/*     d4est_databox_get_field_on_element(ed, name, type, field_on_local_mesh_if_available, dmd); */
/*   } */
/*   else { */
/*     d4est_databox_get_field_on_ghost(ed, ghost_name_id, name, dgd); */
/*   } */

/* } */



d4est_databox_t*
d4est_databox_init
(
 int mpi_rank,
 int* loc_sizes
)
{
  d4est_databox_t* s = D4EST_ALLOC(d4est_databox_t, 1);
  s->mpi_rank = mpi_rank;

  for (int i = 0; i < D4EST_FIELD_TYPES; i++){
    s->field_sizes[i] = loc_sizes[i];
  }
  
  d4est_dictionary_init(&s->fields);
  
  return s;
}

void
d4est_databox_destroy
(
 d4est_databox_t *s
)
{
  /* delete all fields */
  d4est_dictionary_allprefixed_ptr(&s->fields, "",
                                   &d4est_databox_free_fields_callback, NULL);

  d4est_dictionary_clear(&s->fields);
  s->mpi_rank = -1;
  for (int i = 0; i < D4EST_FIELD_TYPES; i++){
    s->field_sizes[i] = -1;
  }

  D4EST_FREE(s);
}
