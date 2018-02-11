include (ExternalProject)

#
# Bundled openblas paths.
#
set(OPENBLAS_BUNDLED_PREFIX "${PROJECT_BINARY_DIR}/ThirdParty/OpenBLAS/install")
set(OPENBLAS_BUNDLED_LIB    "${OPENBLAS_BUNDLED_PREFIX}/lib/libopenblas.a")

macro(openblas_use_bundled)
  set(OPENBLAS_PREFIX "${OPENBLAS_BUNDLED_PREFIX}")
  set(OPENBLAS_INCLUDE_DIRS "${OPENBLAS_BUNDLED_PREFIX}/include")
  set(OPENBLAS_LIBRARIES "${OPENBLAS_BUNDLED_LIB}")
  set(ENABLE_BUNDLED_OPENBLAS True)
endmacro()

macro(openblas_try_system)
  find_package(BLAS REQUIRED)
  message(STATUS "Found a system-wide openblas.")
endmacro()

macro(openblas_try_prefix)
  find_path(OPENBLAS_INCLUDE_DIRS openblas.h ${OPENBLAS_PREFIX}/include NO_DEFAULT_PATH)
  find_library(OPENBLAS_LIB libopenblas.a ${OPENBLAS_PREFIX}/lib NO_DEFAULT_PATH)

  if(OPENBLAS_INCLUDE_DIRS AND OPENBLAS_LIB)
    set(OPENBLAS_LIBRARIES ${OPENBLAS_LIB})
  else()
    message(FATAL_ERROR "Couldn't find openblas in '${OPENBLAS_PREFIX}'")
  endif()
endmacro()

if (OPENBLAS_PREFIX)
  openblas_try_prefix()
elseif (NOT ENABLE_BUNDLED_OPENBLAS)
  openblas_try_system()
else()
  openblas_use_bundled()
endif()

include_directories(${OPENBLAS_INCLUDE_DIRS})

message(STATUS "Use openblas includes: ${OPENBLAS_INCLUDE_DIRS}")
message(STATUS "Use openblas library: ${OPENBLAS_LIBRARIES}")

macro(openblas_build)
  ExternalProject_Add(openblas
    PREFIX              ${CMAKE_BINARY_DIR}/ThirdParty/OpenBLAS
    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/ThirdParty/OpenBLAS
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       cd ${CMAKE_SOURCE_DIR}/ThirdParty/OpenBLAS && make PREFIX=${OPENBLAS_BUNDLED_PREFIX}
    INSTALL_COMMAND     cd ${CMAKE_SOURCE_DIR}/ThirdParty/OpenBLAS && make PREFIX=${OPENBLAS_BUNDLED_PREFIX} install
    CMAKE_ARGS
    -DCMAKE_CXX_FLAGS:STRING=-Wunused-but-set-variable -Wunused-variable -Wincompatible-pointer-types -Wformat-nonliteral -Wstrict-overflow -Wconversion -Wunused-but-set-variable -Wunsafe-loop-optimizations -Wunused-parameter -Wlarger-than= -Wdiscarded-qualifiers -Wfloat-equal -Wmaybe-uninitialized -Wcomment
    -DCMAKE_C_FLAGS:STRING=-Wunused-but-set-variable -Wunused-variable -Wincompatible-pointer-types -Wformat-nonliteral -Wstrict-overflow -Wconversion -Wunused-but-set-variable -Wunsafe-loop-optimizations -Wunused-parameter -Wlarger-than= -Wdiscarded-qualifiers -Wfloat-equal -Wmaybe-uninitialized -Wcomment
    )
    add_dependencies(build_bundled_libs openblas)
endmacro()


if(ENABLE_BUNDLED_OPENBLAS)
  openblas_build()
endif()
