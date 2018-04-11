#ifndef D4EST_CRITBIT_H
#define D4EST_CRITBIT_H 

/* Augmented version of a critbit tree in the
 * BFAM library (https://github.com/bfam/bfam)
 * @file   d4est_critbit.h
 * 
 * @brief  
 * 
 * 
 */

typedef struct
{
  void *root;
} d4est_critbit0_tree_t;

void
d4est_critbit0_clear(d4est_critbit0_tree_t*t);
int
d4est_critbit0_allprefixed(d4est_critbit0_tree_t *t, const char *prefix,
                          int(*handle)(const char*,void*), void *arg);
int
d4est_critbit0_delete(d4est_critbit0_tree_t*t,const char*u);
int
d4est_critbit0_insert(d4est_critbit0_tree_t *t, const char *u);

#endif
