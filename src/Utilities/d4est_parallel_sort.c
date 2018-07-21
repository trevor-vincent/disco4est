#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <mpsort.h>

static void radix_double(const void * ptr, void * radix, void * arg) {
    *(double*)radix = *(const double*) ptr;
}

void
d4est_parallel_sort(sc_MpiComm mpi_comm, int mpirank, double* array,  int local_size){


      mpsort_mpi(array, mysize, sizeof(double),
            radix_double, sizeof(double),
            NULL, mpi_comm);

}
