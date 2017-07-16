#ifndef D4EST_AMR_RANDOM_H
#define D4EST_AMR_RANDOM_H 

/* This file was automatically generated.  Do not edit! */
void d4est_amr_random_destroy(d4est_amr_scheme_t *scheme);
d4est_amr_scheme_t *d4est_amr_random_init(p4est_t *p4est,const char *input_file,d4est_amr_scheme_t *scheme,void *user);
void d4est_amr_random_mark_elements(p4est_iter_volume_info_t *info,void *user_data);


#endif
