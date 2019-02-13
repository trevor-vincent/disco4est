set(CMAKE_BUILD_TYPE Release)

set(ENABLE_BUNDLED_ZLOG ON)
set(ENABLE_BUNDLED_P4EST ON)
set(ENABLE_BUNDLED_ZLIB ON)

set(ENABLE_BUNDLED_HDF5 ON)

set(USE_TCMALLOC ON)
set(ENABLE_BUNDLED_TCMALLOC ON)

#for intel
set(BLA_VENDOR Intel10_64lp_seq)
set(ENABLE_BUNDLED_BLAS OFF)
set(PETSC_PREFIX "/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/MPI/intel2016.4/openmpi2.1/petsc/3.9.0/")
SET(CMAKE_CXX_COMPILER "icpc")
SET(CMAKE_C_COMPILER "icc")
