#include "./pXest/pXest.h"
#include "./pXest/pXest_input.h"
#include "./Problems/problem.h"
 #include <petscsnes.h>

int main(int argc, char *argv[])
{  
  sc_MPI_Comm mpicomm;
  int mpiret = sc_MPI_Init(&argc, &argv);
  SC_CHECK_MPI(mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  PETSC_COMM_WORLD = mpicomm;
  PetscInitialize(&argc,&argv,(char*)0,NULL);

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
#else
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DIM = 2\n");
#endif

  if (proc_rank == 0){
    printf("[D4EST_INFO]: mpisize = %d\n", proc_size);
    printf("[D4EST_INFO]: min_quadrants = %d\n", pXest_input.min_quadrants);
    printf("[D4EST_INFO]: min_level = %d\n", pXest_input.min_level);
    printf("[D4EST_INFO]: fill_uniform = %d\n", pXest_input.fill_uniform);
  }
  p4est_connectivity_t* conn = problem_build_conn();
  p4est_t* p4est = problem_build_p4est
                   (
                    mpicomm,
                    conn,
                    pXest_input.min_quadrants,
                    pXest_input.min_level,
                    pXest_input.fill_uniform
                   );

  
  p4est_geometry_t* p4est_geom = problem_build_geom(conn);

  
  /* start just-in-time dg-math */
  dgmath_jit_dbase_t* dgmath_jit_dbase = dgmath_jit_dbase_init();

  /* checkpoints not supported yet */
  int load_from_checkpoint = 0;
  
  /* Solve Problem */
  problem_init
    (
     argc,
     argv,
     p4est,
     p4est_geom,
     dgmath_jit_dbase,
     proc_size,
     mpicomm,
     load_from_checkpoint
    );
  
  dgmath_jit_dbase_destroy(dgmath_jit_dbase);
  
  /* free pXest */
  p4est_destroy(p4est);
  p4est_connectivity_destroy(conn);

  if(p4est_geom != NULL)
    p4est_geometry_destroy (p4est_geom);

  sc_finalize ();
  PetscFinalize();
  
  return 0;
}
