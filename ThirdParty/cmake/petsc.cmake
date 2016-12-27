include(ProcessorCount)
ProcessorCount(N)
# if(NOT N EQUAL 0)
#   set(CTEST_BUILD_FLAGS -j${N})
#   set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
# endif()

#
# Bundled petsc paths.
#
set(PETSC_BUNDLED_PREFIX "${PROJECT_BINARY_DIR}/ThirdParty/petsc/install")
set(PETSC_BUNDLED_LIBRARIES
  ${PETSC_BUNDLED_PREFIX}/lib/libpetsc.a
  )



macro(petsc_use_bundled)
  set(PETSC_PREFIX "${PETSC_BUNDLED_PREFIX}")
  set(PETSC_INCLUDE_DIRS "${PETSC_BUNDLED_PREFIX}/include")
  set(PETSC_LIBRARIES "${PETSC_BUNDLED_LIBRARIES}")
  set(ENABLE_BUNDLED_PETSC True)
endmacro()

#
# Check if there is a usable petsc at the given prefix path.
#
macro (petsc_try_prefix)

  if(PETSC_INCLUDE_DIRS AND PETSC_LIB)
    # set(PETSC_INCLUDE_DIRS ${PETSC_PETSC_INCLUDE_DIR} ${PETSC_SC_INCLUDE_DIRS})
    set(PETSC_LIBRARIES ${PETSC_LIB})
    include_directories(${PETSC_INCLUDE_DIRS})
  else()
    message(FATAL_ERROR "Couldn't find petsc in '${PETSC_PREFIX}'")
  endif()
endmacro()

#
# petsc options.
#
option(ENABLE_BUNDLED_PETSC "Enable building of the bundled petsc" ON)
option(PETSC_PREFIX "Build with petsc at the given path" "")

if(PETSC_PREFIX AND ENABLE_BUNDLED_PETSC)
  message(FATAL_ERROR "Options PETSC_PREFIX and ENABLE_BUNDLED_PETSC "
    "are not compatible with each other.")
endif()

if (PETSC_PREFIX)
  petsc_try_prefix()
else()
  petsc_use_bundled()
endif()

include_directories(${PETSC_INCLUDE_DIRS})

message(STATUS "Use petsc includes: ${PETSC_INCLUDE_DIRS}")
message(STATUS "Use petsc library: ${PETSC_LIBRARIES}")

macro(petsc_build)
  if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    set(petsc_config_args "--with-debugging=1")
  else("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    set(petsc_config_args "--with-debugging=0")
  endif("${CMAKE_BUILD_TYPE}" MATCHES "Debug")


  set(blas_config_args "")
  if(USE_OPENBLAS)
    message("*****USING OPENBLAS*****")
    set(blas_config_args "--with-openblas=1 --with-openblas-dir=${OPENBLAS_ROOT} --with-openblas-include=${OPENBLAS_INCLUDE_DIR} --with-openblas-lib=${OPENBLAS_LIBRARIES} --with-blas-lib=${OPENBLAS_LIBRARIES}")
  endif()
  if(DOWNLOAD_BLAS)
    set(blas_config_args "--download-openblas=yes --download-openblas-make-options=USE_THREAD=0")
  endif()
  
  message("******blas_config_args =  ${blas_config_args}")
  ExternalProject_Add(petsc
    PREFIX    ${CMAKE_BINARY_DIR}/ThirdParty/petsc
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/ThirdParty/petsc/
    CONFIGURE_COMMAND cd ${CMAKE_SOURCE_DIR}/ThirdParty/petsc &&
    python configure
    # "CC=${MPI_C_COMPILER}"
    # "CXX=${MPI_CXX_COMPILER}"
    # "F77=${MPI_Fortran_COMPILER}"
    # "CPPFLAGS=-I${LUA_INCLUDE_DIR} ${zlib_include}"
    # "LIBS=${lua_lib} ${zlib_lib}"
    ${petsc_config_args}
    ${blas_config_args}
    # --enable-mpi
    # --with-fc=gfortran
    --with-x=0
    --with-ssl=0
    # --with-cc=icc
    # --download-fblaslapack
    --with-shared-libraries=0
    --prefix=${PETSC_BUNDLED_PREFIX}
    BUILD_COMMAND       cd ${CMAKE_SOURCE_DIR}/ThirdParty/petsc && make MAKE_NP=${N}
    INSTALL_COMMAND     cd ${CMAKE_SOURCE_DIR}/ThirdParty/petsc && make install
    )
    # add_dependencies(petsc petsc_bundled_libs)
add_dependencies(build_bundled_libs petsc)

endmacro()

if (ENABLE_BUNDLED_PETSC)
  petsc_build()
endif()
