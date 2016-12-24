#include "./pXest/pXest.h"
#include "./pXest/pXest_input.h"
#include "./Problems/problem.h"
 #include <petscsnes.h>

int main(int argc, char *argv[])
{  
  int mpiret = sc_MPI_Init(&argc, &argv);
  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;  
  SC_CHECK_MPI(mpiret);

  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);

  pXest_input_t pXest_input = pXest_input_parse("options.input");
  
#ifdef DEBUG
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DEBUG MODE ON\n");
  p4est_init(NULL, SC_LP_ALWAYS);
#else
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DEBUG MODE OFF\n");
  p4est_init(NULL, SC_LP_ERROR);
#endif 

#if (P4EST_DIM)==3
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DIM = 3\n");
  p4est_connectivity_t* conn = p8est_connectivity_new_unitcube();
#else
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DIM = 2\n");
  p4est_connectivity_t* conn = p4est_connectivity_new_unitsquare();
#endif

  problem_help();
  
  p4est_t* p4est = p4est_new_ext
    (
     mpicomm,
     conn,
     pXest_input.min_quadrants,
     pXest_input.min_level,
     pXest_input.fill_uniform,
     0,
     NULL,
     NULL
    );

  /* free pXest */
  p4est_destroy(p4est);
  p4est_connectivity_destroy(conn);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
