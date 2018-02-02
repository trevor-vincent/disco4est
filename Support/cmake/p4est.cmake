include(ProcessorCount)
#
# Bundled p4est paths.
#
set(P4EST_BUNDLED_PREFIX "${PROJECT_BINARY_DIR}/ThirdParty/p4est/install")
set(P4EST_BUNDLED_LIBRARIES
  ${P4EST_BUNDLED_PREFIX}/lib/libp4est.a
  ${P4EST_BUNDLED_PREFIX}/lib/libsc.a
  )

macro(p4est_use_bundled)
  set(P4EST_PREFIX "${P4EST_BUNDLED_PREFIX}")
  set(P4EST_INCLUDE_DIRS "${P4EST_BUNDLED_PREFIX}/include")
  set(P4EST_LIBRARIES "${P4EST_BUNDLED_LIBRARIES}")
  set(ENABLE_BUNDLED_P4EST True)
endmacro()

macro(p4est_try_system)
  find_path(P4EST_SC_INCLUDE_DIRS sc.h PATH_SUFFIXES include)
  find_library(P4EST_SC_LIB NAMES sc PATH_SUFFIXES lib)

  find_path(P4EST_INCLUDE_DIRS p4est.h PATH_SUFFIXES include)
  find_library(P4EST_LIB NAMES p4est PATH_SUFFIXES lib)

  if(P4EST_INCLUDE_DIRS AND P4EST_LIB AND
     P4EST_SC_INCLUDE_DIRS AND P4EST_SC_LIB)
    message (STATUS "-- sc include: ${P4EST_SC_INCLUDE_DIR}, lib: ${P4EST_SC_LIB}")
    message (STATUS "-- p4est include: ${P4EST_INCLUDE_DIR}, lib: ${P4EST_LIB}")
    message (STATUS "Found a system-wide p4est.")
    set(P4EST_INCLUDE_DIRS ${P4EST_INCLUDE_DIR} ${P4EST_SC_INCLUDE_DIRS})
    set(P4EST_LIBRARIES ${P4EST_LIB} ${P4EST_SC_LIB})
  else()
    message (FATAL_ERROR "Not found a system p4est")
    #p4est_use_bundled()
  endif()
endmacro()

#
# Check if there is a usable p4est at the given prefix path.
#
macro (p4est_try_prefix)
  find_path(P4EST_SC_INCLUDE_DIRS sc.h ${P4EST_PREFIX}/include)
  find_library(P4EST_SC_LIB libsc.a ${P4EST_PREFIX}/lib)
  find_path(P4EST_INCLUDE_DIRS p4est.h ${P4EST_PREFIX}/include)
  find_library(P4EST_LIB libp4est.a ${P4EST_PREFIX}/lib)

  if(P4EST_INCLUDE_DIRS AND P4EST_LIB AND
     P4EST_SC_INCLUDE_DIRS AND P4EST_SC_LIB)
    set(P4EST_INCLUDE_DIRS ${P4EST_INCLUDE_DIRS} ${P4EST_SC_INCLUDE_DIRS})
    set(P4EST_LIBRARIES ${P4EST_LIB} ${P4EST_SC_LIB})
  else()
    message(FATAL_ERROR "Couldn't find p4est in '${P4EST_PREFIX}'")
  endif()
endmacro()

#
# p4est options.
#
# option(ENABLE_BUNDLED_P4EST "Enable building of the bundled p4est" ON)
# option(P4EST_PREFIX "Build with p4est at the given path" "")

if(P4EST_PREFIX AND ENABLE_BUNDLED_P4EST)
  message(FATAL_ERROR "Options P4EST_PREFIX and ENABLE_BUNDLED_P4EST "
    "are not compatible with each other.")
endif()

if (P4EST_PREFIX)
  # trying to build with specified p4est.
  p4est_try_prefix()
elseif (NOT ENABLE_BUNDLED_P4EST)
  # trying to build with system p4est, macro can turn on
  # building of p4est bundled with the server source.
  p4est_try_system()
else()
  p4est_use_bundled()
endif()

include_directories(${P4EST_INCLUDE_DIRS})

message(STATUS "Use p4est includes: ${P4EST_INCLUDE_DIRS}")
message(STATUS "Use p4est library: ${P4EST_LIBRARIES}")

set(num_of_build_threads "-j${N}")

if(GRAHAM)
  set(num_of_build_threads "-j1")
endif()
  

macro(p4est_build)  
  foreach(dir ${ZLIB_INCLUDE_DIRS})
    set(zlib_include "${zlib_include} -I${dir}")
  endforeach()

  foreach(lib ${ZLIB_LIBRARIES})
    set(zlib_lib "${zlib_lib} ${lib}")
  endforeach()

  if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    set(p4est_config_args "--enable-debug")
  else("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    set(p4est_config_args "")
  endif("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  ExternalProject_Add(p4est
    PREFIX    ${CMAKE_BINARY_DIR}/ThirdParty/p4est
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/ThirdParty/p4est/
    CONFIGURE_COMMAND cd ${CMAKE_SOURCE_DIR}/ThirdParty/p4est && git submodule init && git submodule update && ./bootstrap &&
    ./configure
    "CC=${MPI_C_COMPILER}"
    "CXX=${MPI_CXX_COMPILER}"
    # "F77=${MPI_Fortran_COMPILER}"
    # "CPPFLAGS=-fPIC"
    # "LDFLAGS=-fPIC"
    # "LIBS=${zlib_lib}"
    ${p4est_config_args}
    --enable-mpi
    --without-blas
    # --enable-mpi --disable-vtk-binary --without-blas
    # --enable-static
    # --disable-shared
    # --without-zlib
    --prefix=${P4EST_BUNDLED_PREFIX}
    BUILD_COMMAND       cd ${CMAKE_SOURCE_DIR}/ThirdParty/p4est && make ${num_of_build_threads}
    INSTALL_COMMAND     cd ${CMAKE_SOURCE_DIR}/ThirdParty/p4est && make install
  )
  add_dependencies(p4est p4est_bundled_libs)
  add_dependencies(build_bundled_libs p4est)
endmacro()

if (ENABLE_BUNDLED_P4EST)
  p4est_build()
endif()
