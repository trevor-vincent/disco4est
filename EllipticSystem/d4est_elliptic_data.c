#include "problem_data.h"

void
problem_data_copy_ptrs
(
 d4est_elliptic_problem_data_t* pd1,
 d4est_elliptic_problem_data_t* pd2
)
{
  pd2->Au = pd1->Au;
  pd2->u = pd1->u;
  pd2->u0 = pd1->u0;
  pd2->rhs = pd1->rhs;
  pd2->local_nodes = pd1->local_nodes;
  pd2->mpi_rank = pd1->mpi_rank;
  pd2->user = pd1->user;
}
