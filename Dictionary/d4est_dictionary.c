#define D4EST_BUFSIZ 8192
#define D4EST_KEYVALUE_SPLIT '\255'
#define D4EST_PTR_STR_LEN D4EST_BUFSIZ

#include "d4est_dictionary.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void d4est_dictionary_init(d4est_dictionary_t *d)
{
  d->num_entries = 0;
  d->t.root = NULL;
}

void d4est_dictionary_clear(d4est_dictionary_t *d)
{
  d4est_critbit0_clear(&(d->t));
  d->num_entries = 0;
}

static int d4est_dictionary_contains_check(const char *value, void *arg)
{
  (*(int *)arg)++;
  return 1;
}

int d4est_dictionary_contains(d4est_dictionary_t *d, const char *key)
{
  const size_t keylen = strlen(key);

  char *u = malloc(sizeof(char) * (keylen + 2));

  memcpy(u, key, keylen);
  u[keylen] = D4EST_KEYVALUE_SPLIT;
  u[keylen + 1] = '\0';
  int found = 0;
  d4est_critbit0_allprefixed(&(d->t), u, &d4est_dictionary_contains_check,
                            &found);
  free(u);
  return found;
}

/** Inserting key and value pair into a dictionary
 *
 * It takes a dictionary, \a d, and possibly mutates it such that a \c NULL
 * terminated strings, \a key and \a value, is a member on exit.
 *
 * \param [in,out] d dictionary
 * \param [in] key possible key
 * \param [in] val possible value
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it key} {\rm was already a member} \cr
 *            2 &if {\it d} {\rm was mutated successfully}}$.
 */
int d4est_dictionary_insert(d4est_dictionary_t *d, const char *key,
                                  const char *val)
{
  const size_t keylen = strlen(key);
  const size_t vallen = strlen(val);

  char *keyval = malloc(sizeof(char) * (keylen + vallen + 2));

  memcpy(keyval, key, keylen);

  if (d4est_dictionary_contains(d, key))
  {
    free(keyval);
    return 1;
  }

  keyval[keylen] = D4EST_KEYVALUE_SPLIT;
  memcpy(&keyval[keylen + 1], val, vallen);
  keyval[keylen + vallen + 1] = '\0';

  int rval = d4est_critbit0_insert(&(d->t), keyval);

  if (rval == 2)
    ++d->num_entries;

  free(keyval);
  return rval;
}

int d4est_dictionary_insert_ptr(d4est_dictionary_t *d, const char *key,
                               const void *val_ptr)
{
  char val_str[D4EST_PTR_STR_LEN + 1];
  snprintf(val_str, D4EST_PTR_STR_LEN + 1, "%p", val_ptr);
  return d4est_dictionary_insert(d, key, val_str);
}

int d4est_dictionary_insert_int(d4est_dictionary_t *d, const char *key,
                               const int val)
{
  char val_str[D4EST_BUFSIZ];
  snprintf(val_str, D4EST_BUFSIZ, "%d", val);
  return d4est_dictionary_insert(d, key, val_str);
}

int d4est_dictionary_get_value_int(d4est_dictionary_t *d, const char *key,
                                         int *val)
{
  char *val_str = d4est_dictionary_get_value(d, key);
  if (val_str == NULL)
    return 0;
  int n = sscanf(val_str, "%d", val);
  return n;
}


static int d4est_dictionary_get_value_handle(const char *keyval, void *arg)
{
  char *key = (char *)((void **)arg)[0];
  char **val = (char **)((void **)arg)[1];
  const size_t keylen = strlen(key);

  *val = (char *)&keyval[keylen];

  return 1;
}

/** Return a value given a key
 *
 * It takes a dictionary, \a d, returns a pointer to the value associated with
 * a \c NULL terminated \a key.
 *
 * \param [in] d dictionary
 * \param [in] key possible key
 * \returns:
 *   $\cases{ \c NULL & if {\it key} {\rm is not a member} \cr
 *          {\rm pointer to value} & if {\it key} {\rm is a member}$
 */
char *d4est_dictionary_get_value(d4est_dictionary_t *d, const char *key)
{
  if (!d4est_dictionary_contains(d, key))
    return NULL;

  const size_t keylen = strlen(key);

  char *u = malloc(sizeof(char) * (keylen + 2));

  memcpy(u, key, keylen);
  u[keylen] = D4EST_KEYVALUE_SPLIT;
  u[keylen + 1] = '\0';

  char *value = NULL;
  void *arg[2];
  arg[0] = u;
  arg[1] = &value;

  d4est_critbit0_allprefixed(&(d->t), u, &d4est_dictionary_get_value_handle, arg);

  free(u);

  return value;
}

/** Delete key and value pair into a dictionary
 *
 * It takes a dictionary, \a d, and possibly mutates it such that a \c NULL
 * terminated string, \a key  is not member on exit.
 *
 * \param [in,out] d dictionary
 * \param [in] key possible key
 * \returns:
 *   $\cases{ 0 &if {\rm key not found} \cr
 *            1 &if {\it key deleted} }$.
 */
/* static int d4est_dictionary_delete(d4est_dictionary_t *d, const char *key) */
/* { */
/*   const size_t keylen = strlen(key); */
/*   char *val = d4est_dictionary_get_value(d, key); */
/*   if (val == NULL) */
/*     return 0; */
/*   const char *keyval = val - keylen - 1; */
/*   return d4est_critbit0_delete(&(d->t), keyval); */
/* } */

void *d4est_dictionary_get_value_ptr(d4est_dictionary_t *d, const char *key)
{
  char *val_str = d4est_dictionary_get_value(d, key);
  if (val_str == NULL)
    return NULL;
  void *val_ptr = NULL;
  sscanf(val_str, "%p", &val_ptr);
  return val_ptr;
}

typedef struct
{
  int (*handle)(const char *, const char *, void *);
  void *arg;
} d4est_dict_allprex;

static int d4est_dictionary_allprefixed_usercall(const char *keyval, void *arg)
{
  char *key = (char *)keyval;
  char *split = strchr(key, D4EST_KEYVALUE_SPLIT);
  *split = '\0';
  d4est_dict_allprex *s_arg = (d4est_dict_allprex *)arg;
  s_arg->handle(key, split + 1, s_arg->arg);
  *split = D4EST_KEYVALUE_SPLIT;
  return 1;
}

int d4est_dictionary_allprefixed(d4est_dictionary_t *d, const char *prefix,
                                int (*handle)(const char *, const char *,
                                              void *),
                                void *arg)
{
  d4est_dict_allprex args = {0, 0};
  args.handle = handle;
  args.arg = arg;
  return d4est_critbit0_allprefixed(
      &(d->t), prefix, &d4est_dictionary_allprefixed_usercall, &args);
}

typedef struct
{
  int (*handle)(const char *, void *, void *);
  void *arg;
} d4est_dict_allprex_ptr;

static int d4est_dictionary_allprefixed_usercall_ptr(const char *keyval,
                                                    void *arg)
{
  char *key = (char *)keyval;
  char *split = strchr(key, D4EST_KEYVALUE_SPLIT);
  *split = '\0';

  void *val_ptr = NULL;
  sscanf(split + 1, "%p", &val_ptr);

  d4est_dict_allprex_ptr *s_arg = (d4est_dict_allprex_ptr *)arg;
  s_arg->handle(key, val_ptr, s_arg->arg);

  *split = D4EST_KEYVALUE_SPLIT;
  return 1;
}

int d4est_dictionary_allprefixed_ptr(d4est_dictionary_t *d, const char *prefix,
                                    int (*handle)(const char *, void *, void *),
                                    void *arg)
{
  d4est_dict_allprex_ptr args = {0, 0};
  args.handle = handle;
  args.arg = arg;
  return d4est_critbit0_allprefixed(
      &(d->t), prefix, &d4est_dictionary_allprefixed_usercall_ptr, &args);
}

int d4est_dictionary_delete(d4est_dictionary_t *d, const char *key)
{
  const size_t keylen = strlen(key);
  char *val = d4est_dictionary_get_value(d, key);
  if (val == NULL)
    return 0;
  const char *keyval = val - keylen - 1;
  return d4est_critbit0_delete(&(d->t), keyval);
}

// }}}
