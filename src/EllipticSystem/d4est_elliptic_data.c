#include <d4est_elliptic_data.h>

void
d4est_elliptic_data_copy_ptrs
(
 d4est_elliptic_data_t* pd1,
 d4est_elliptic_data_t* pd2
)
{
  pd2->Au = pd1->Au;
  pd2->u = pd1->u;
  pd2->u0 = pd1->u0;
  pd2->rhs = pd1->rhs;
  pd2->local_nodes = pd1->local_nodes;
  pd2->mpirank = pd1->mpirank;
  pd2->user = pd1->user;
  pd2->field_types = pd1->field_types;
  pd2->num_of_fields = pd1->num_of_fields;
}
