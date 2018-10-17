#ifndef D4EST_FIELD_H
#define D4EST_FIELD_H 

#include <pXest.h>

#define D4EST_FIELD_TYPES 4

/* more types will be added in the future */

typedef enum {
              VOLUME, /* a field on each element volume */
              NODAL, /* a field on the LGL points of each element  */
              FACE, /* a field on each face of an element */
              MORTAR, /* a field on each mortar face touching the element, per element size = (P4EST_FACES)*(P4EST_HALF) */
              NO_TYPE
} d4est_field_type_t;

typedef struct {

  double* field_data;
  d4est_field_type_t field_type;
  
} d4est_field_t;

typedef enum {

              FIELD_ZEROED,
              FIELD_NOT_ZEROED
} d4est_field_init_type_t;


#define D4EST_FIELD_CHECK_TYPE(a) do {                  \
    D4EST_ASSERT(a >= 0 && a < D4EST_FIELD_TYPES);      \
  } while(0)


#endif
