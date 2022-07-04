include (ExternalProject)

#
# Bundled tcmalloc paths.
#
set(TCMALLOC_BUNDLED_PREFIX "${PROJECT_BINARY_DIR}/third_party/gperftools/install")
set(TCMALLOC_BUNDLED_LIB    "${TCMALLOC_BUNDLED_PREFIX}/lib/libtcmalloc.a")

macro(tcmalloc_use_bundled)
  set(TCMALLOC_PREFIX "${TCMALLOC_BUNDLED_PREFIX}")
  set(TCMALLOC_INCLUDE_DIRS "${TCMALLOC_BUNDLED_PREFIX}/include")
  set(TCMALLOC_LIBRARIES "${TCMALLOC_BUNDLED_LIB}")
  set(ENABLE_BUNDLED_TCMALLOC True)
endmacro()

macro(tcmalloc_try_system)
  find_package(TCMALLOC REQUIRED)
  message(STATUS "Found a system-wide tcmalloc.")
endmacro()

macro(tcmalloc_try_prefix)
  find_path(TCMALLOC_INCLUDE_DIRS tcmalloc.h ${TCMALLOC_PREFIX}/include NO_DEFAULT_PATH)
  find_library(TCMALLOC_LIB libtcmalloc.a ${TCMALLOC_PREFIX}/lib NO_DEFAULT_PATH)

  if(TCMALLOC_INCLUDE_DIRS AND TCMALLOC_LIB)
    set(TCMALLOC_LIBRARIES ${TCMALLOC_LIB})
  else()
    message(FATAL_ERROR "Couldn't find tcmalloc in '${TCMALLOC_PREFIX}'")
  endif()
endmacro()

if (TCMALLOC_PREFIX)
  tcmalloc_try_prefix()
elseif (NOT ENABLE_BUNDLED_TCMALLOC)
  tcmalloc_try_system()
else()
  tcmalloc_use_bundled()
endif()

include_directories(${TCMALLOC_INCLUDE_DIRS})

message(STATUS "Use tcmalloc includes: ${TCMALLOC_INCLUDE_DIRS}")
message(STATUS "Use tcmalloc library: ${TCMALLOC_LIBRARIES}")

macro(tcmalloc_build)
  ExternalProject_Add(tcmalloc
    PREFIX              ${CMAKE_BINARY_DIR}/third_party/gperftools
    URL          ${CMAKE_SOURCE_DIR}/third_party/gperftools-gperftools-2.10.tar.gz
    CONFIGURE_COMMAND   cd ${CMAKE_SOURCE_DIR}/third_party/gperftools && libtoolize --force && aclocal && autoheader && automake --add-missing && autoconf && ./configure --prefix=${TCMALLOC_BUNDLED_PREFIX}
    BUILD_COMMAND       cd ${CMAKE_SOURCE_DIR}/third_party/gperftools && make PREFIX=${TCMALLOC_BUNDLED_PREFIX}
    INSTALL_COMMAND     cd ${CMAKE_SOURCE_DIR}/third_party/gperftools && make PREFIX=${TCMALLOC_BUNDLED_PREFIX} install
  )
  set_target_properties(tcmalloc PROPERTIES EXCLUDE_FROM_ALL ON)
  add_dependencies(build_bundled_libs tcmalloc)
endmacro()


if(ENABLE_BUNDLED_TCMALLOC)
  tcmalloc_build()
endif()
