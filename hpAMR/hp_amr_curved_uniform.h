#ifndef HP_AMR_CURVED_UNIFORM_H
#define HP_AMR_CURVED_UNIFORM_H 

hp_amr_scheme_t *hp_amr_curved_uniform_init();
void hp_amr_curved_uniform_destroy(hp_amr_scheme_t *scheme);
void hp_amr_curved_uniform_set_refinement(p4est_iter_volume_info_t *info,void *user_data);

#endif
