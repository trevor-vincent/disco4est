if(LAPTOP)
  set(OPENBLAS_DIR "/home/tvincent/Dropbox/Research/Codes/OpenBLAS/")
  set(OPENBLAS_ROOT "/home/tvincent/Dropbox/Research/Codes/OpenBLAS/")
  set(OPENBLAS_LIBRARIES
    "/home/tvincent/Dropbox/Research/Codes/OpenBLAS/libopenblas.a"
    )
  set(OPENBLAS_INCLUDE_DIR "/home/tvincent/Dropbox/Research/Codes/OpenBLAS/lapack-netlib/LAPACKE/include/")  
  set(USE_TCMALLOC ON)
  set(TCMALLOC_DIR "/home/tvincent/Dropbox/Research/Codes/gperftools/gperftools_install/")
  set(TCMALLOC_ROOT "/home/tvincent/Dropbox/Research/Codes/gperftools/gperftools_install/")
  set(TCMALLOC_LIBRARIES "/home/tvincent/Dropbox/Research/Codes/gperftools/gperftools_install/lib/libtcmalloc.a")
  set(TCMALLOC_INCLUDE_DIR "/home/tvincent/Dropbox/Research/Codes/gperftools/gperftools_install/include/gperftools/")

  set(HDF5_LIBRARIES "/home/tvincent/Dropbox/Research/Codes/hdf5-1.10.1/install/lib/libhdf5-static.a")
  set(HDF5_INCLUDE_DIR "/home/tvincent/Dropbox/Research/Codes/hdf5-1.10.1/install/include/")
endif()

if(LAPTOP_WITH_PREFIX)
  set(PETSC_PREFIX "/home/tvincent/Dropbox/Research/Codes/disco4est/src_new_build_debug/ThirdParty/petsc/install/")
  set(P4EST_PREFIX "/home/tvincent/Dropbox/Research/Codes/disco4est/src_new_build_debug/ThirdParty/p4est/install/")
  set(ZLIB_PREFIX
    "/home/tvincent/Dropbox/Research/Codes/disco4est/src_new_build_debug/ThirdParty/zlib")
  set(OPENBLAS_DIR "/home/tvincent/Dropbox/Research/Codes/OpenBLAS/")
  set(OPENBLAS_ROOT "/home/tvincent/Dropbox/Research/Codes/OpenBLAS/")
  set(OPENBLAS_LIBRARIES
    "/home/tvincent/Dropbox/Research/Codes/OpenBLAS/libopenblas.a"
    )
  set(OPENBLAS_INCLUDE_DIR "/home/tvincent/Dropbox/Research/Codes/OpenBLAS/lapack-netlib/LAPACKE/include/")  
  set(USE_TCMALLOC ON)
  set(TCMALLOC_DIR "/home/tvincent/Dropbox/Research/Codes/gperftools/gperftools_install/")
  set(TCMALLOC_ROOT "/home/tvincent/Dropbox/Research/Codes/gperftools/gperftools_install/")
  set(TCMALLOC_LIBRARIES "/home/tvincent/Dropbox/Research/Codes/gperftools/gperftools_install/lib/libtcmalloc.a")
  set(TCMALLOC_INCLUDE_DIR "/home/tvincent/Dropbox/Research/Codes/gperftools/gperftools_install/include/gperftools/")

  set(HDF5_LIBRARIES "/home/tvincent/Dropbox/Research/Codes/hdf5-1.10.1/install/lib/libhdf5-static.a")
  set(HDF5_INCLUDE_DIR "/home/tvincent/Dropbox/Research/Codes/hdf5-1.10.1/install/include/")
endif()


if(SCINET)
  set(OPENBLAS_DIR "/scinet/gpc/Libraries/OpenBLAS/singlethreaded/")
  set(OPENBLAS_ROOT "/scinet/gpc/Libraries/OpenBLAS/singlethreaded/")
  set(OPENBLAS_LIBRARIES
    "/scinet/gpc/Libraries/OpenBLAS/singlethreaded/lib/libopenblas.a"
    )
  set(OPENBLAS_INCLUDE_DIR "/scinet/gpc/Libraries/OpenBLAS/singlethreaded/include")
  set(HDF5_LIBRARIES "/scratch/p/pfeiffer/tvincent/hdf5-1.10.1/install/lib/libhdf5-static.a")
  set(HDF5_INCLUDE_DIR "/scratch/p/pfeiffer/tvincent/hdf5-1.10.1/install/include/") 
endif()

if(GRAHAM OR CEDAR)
  set(OPENBLAS_DIR "/scratch/tvincent/OpenBLAS/")
  set(OPENBLAS_ROOT "/scratch/tvincent/OpenBLAS/")
  set(OPENBLAS_LIBRARIES
    "/scratch/tvincent/OpenBLAS/libopenblas.a"
    )
  set(OPENBLAS_INCLUDE_DIR "/scratch/tvincent/OpenBLAS/lapack-netlib/LAPACKE/include/")
  set(HDF5_LIBRARIES "/scratch/tvincent/hdf5-1.10.1/install/lib/libhdf5-static.a")
  set(HDF5_INCLUDE_DIR "/scratch/tvincent/hdf5-1.10.1/install/include/")
endif()


if(CITA)
  set(OPENBLAS_DIR "/opt/openblas/0.2.19/")
  set(OPENBLAS_ROOT "/opt/openblas/0.2.19/")
  set(OPENBLAS_LIBRARIES
    "/opt/openblas/0.2.19/lib/libopenblas.a"
    )
  set(OPENBLAS_INCLUDE_DIR "/opt/openblas/0.2.19/include")
endif()
