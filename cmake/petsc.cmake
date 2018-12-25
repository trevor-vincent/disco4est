include(ProcessorCount)

#
# Bundled petsc paths.
#
set(PETSC_BUNDLED_PREFIX "${PROJECT_BINARY_DIR}/third_party/petsc/install")
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
  find_path(PETSC_INCLUDE_DIRS petsc.h ${PETSC_PREFIX}/include)
  find_library(PETSC_LIB libpetsc.a libpetsc.so ${PETSC_PREFIX}/lib)

  if(PETSC_INCLUDE_DIRS AND PETSC_LIB)
    set(PETSC_LIBRARIES ${PETSC_LIB})
    include_directories(${PETSC_INCLUDE_DIRS})
  else()
    message(FATAL_ERROR "Couldn't find petsc in '${PETSC_PREFIX}'")
  endif()
endmacro()

#
# petsc options.
#
# option(ENABLE_BUNDLED_PETSC "Enable building of the bundled petsc" ON)
# option(PETSC_PREFIX "Build with petsc at the given path" "")

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
    MESSAGE ("*******PETSC DEBUGGING TURNED ON*******")
  else()
    set(petsc_config_args "--with-debugging=0")
    MESSAGE ("*******PETSC DEBUGGING TURNED OFF*******")
  endif("${CMAKE_BUILD_TYPE}" MATCHES "Debug")

  set(blas_config_args
    "--with-blas-lib=${BLAS_LIBRARIES}"
    )
  set(lapack_config_args
    "--with-lapack-lib=${BLAS_LIBRARIES}"
    )    
  ExternalProject_Add(petsc
    PREFIX    ${CMAKE_BINARY_DIR}/third_party/petsc
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/third_party/petsc/
    CONFIGURE_COMMAND cd ${CMAKE_SOURCE_DIR}/third_party/petsc &&
    python2 configure
    ${petsc_config_args}
    ${blas_config_args}
    ${lapack_config_args}
    --with-x=0
    --with-ssl=0
    --with-make-np=1
    --with-shared-libraries=0
    --with-silent-rules=1
    --prefix=${PETSC_BUNDLED_PREFIX}
    BUILD_COMMAND       cd ${CMAKE_SOURCE_DIR}/third_party/petsc && make -j1 --silent V=0
    INSTALL_COMMAND     cd ${CMAKE_SOURCE_DIR}/third_party/petsc && make install --silent
    )
  if(ENABLE_BUNDLED_BLAS)
    add_dependencies(petsc openblas)
  endif()
  add_dependencies(build_bundled_libs petsc)

endmacro()
if (ENABLE_BUNDLED_PETSC)
  petsc_build()
endif()

