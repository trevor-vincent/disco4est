#ifndef AMR_CURVED_UNIFORM_H
#define AMR_CURVED_UNIFORM_H 

amr_scheme_t *amr_curved_uniform_init();
void amr_curved_uniform_destroy(amr_scheme_t *scheme);
void amr_curved_uniform_set_refinement(p4est_iter_volume_info_t *info,void *user_data);

#endif
