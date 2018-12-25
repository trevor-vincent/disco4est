include (ExternalProject)

#
# Bundled openblas paths.
#
set(BLAS_BUNDLED_PREFIX "${PROJECT_BINARY_DIR}/third_party/OpenBLAS/install")
set(BLAS_BUNDLED_LIB    "${BLAS_BUNDLED_PREFIX}/lib/libopenblas.a")

macro(openblas_use_bundled)
  set(BLAS_PREFIX "${BLAS_BUNDLED_PREFIX}")
  set(BLAS_INCLUDE_DIRS "${BLAS_BUNDLED_PREFIX}/include")
  set(BLAS_LIBRARIES "${BLAS_BUNDLED_LIB}")
  set(ENABLE_BUNDLED_BLAS True)
  openblas_build()
endmacro()

macro(openblas_try_system)
  message(STATUS "******Trying to find blas on system******.")
  find_package(BLAS REQUIRED)
  message(STATUS "Found a system-wide openblas.")
  # message(STATUS "Found a system-wide openblas.")
endmacro()

macro(openblas_try_prefix)
  message(STATUS "*******Trying to use blas prefix*******.")
  find_path(BLAS_INCLUDE_DIRS lapacke.h ${BLAS_PREFIX}/include NO_DEFAULT_PATH)
  find_library(BLAS_LIB libopenblas.a ${BLAS_PREFIX}/lib NO_DEFAULT_PATH)

  if(BLAS_INCLUDE_DIRS AND BLAS_LIB)
    set(BLAS_LIBRARIES ${BLAS_LIB})
  else()
    message(FATAL_ERROR "Couldn't find openblas in '${BLAS_PREFIX}'")
  endif()
endmacro()


macro(openblas_build)
  ExternalProject_Add(openblas
    PREFIX              ${CMAKE_BINARY_DIR}/third_party/OpenBLAS
    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/third_party/OpenBLAS
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       cd ${CMAKE_SOURCE_DIR}/third_party/OpenBLAS && make PREFIX=${BLAS_BUNDLED_PREFIX}
    INSTALL_COMMAND     cd ${CMAKE_SOURCE_DIR}/third_party/OpenBLAS && make PREFIX=${BLAS_BUNDLED_PREFIX} install
    CMAKE_ARGS
    -DCMAKE_CXX_FLAGS:STRING=-Wunused-but-set-variable -Wunused-variable -Wincompatible-pointer-types -Wformat-nonliteral -Wstrict-overflow -Wconversion -Wunused-but-set-variable -Wunsafe-loop-optimizations -Wunused-parameter -Wlarger-than= -Wdiscarded-qualifiers -Wfloat-equal -Wmaybe-uninitialized -Wcomment
    -DCMAKE_C_FLAGS:STRING=-Wunused-but-set-variable -Wunused-variable -Wincompatible-pointer-types -Wformat-nonliteral -Wstrict-overflow -Wconversion -Wunused-but-set-variable -Wunsafe-loop-optimizations -Wunused-parameter -Wlarger-than= -Wdiscarded-qualifiers -Wfloat-equal -Wmaybe-uninitialized -Wcomment
    )
    add_dependencies(build_bundled_libs openblas)
endmacro()


if (BLAS_PREFIX)
  openblas_try_prefix()
elseif (NOT ENABLE_BUNDLED_BLAS)
  openblas_try_system()
else()
  openblas_use_bundled()
endif()

include_directories(${BLAS_INCLUDE_DIRS})
message(STATUS "Use openblas includes: ${BLAS_INCLUDE_DIRS}")
message(STATUS "Use openblas library: ${BLAS_LIBRARIES}")

# if(ENABLE_BUNDLED_BLAS)

# endif()
