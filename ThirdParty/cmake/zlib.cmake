include (ExternalProject)

#
# Bundled zlib paths.
#
set(ZLIB_BUNDLED_PREFIX "${PROJECT_BINARY_DIR}/ThirdParty/zlib/install")
set(ZLIB_BUNDLED_LIB    "${ZLIB_BUNDLED_PREFIX}/lib/libz.a")

macro(zlib_use_bundled)
  set(ZLIB_PREFIX "${ZLIB_BUNDLED_PREFIX}")
  set(ZLIB_INCLUDE_DIRS "${ZLIB_BUNDLED_PREFIX}/include")
  set(ZLIB_LIBRARIES "${ZLIB_BUNDLED_LIB}")
  set(ENABLE_BUNDLED_ZLIB True)
endmacro()

macro(zlib_try_system)
  find_package(ZLIB)
  if(ZLIB_FOUND)
    message(STATUS "include: ${ZLIB_INCLUDE_DIRS}, lib: ${ZLIB_LIBRARIES}")
    message(STATUS "Found a system-wide zlib.")
  else()
    message(FATAL_ERROR "Not found a system zlib")
    #zlib_use_bundled()
  endif()
endmacro()

#
# zlib options.
#
option(ENABLE_BUNDLED_ZLIB "Enable building of the bundled zlib" ON)

if (NOT ENABLE_BUNDLED_ZLIB)
  zlib_try_system()
else()
  zlib_use_bundled()
endif()

message(STATUS "Use zlib includes: ${ZLIB_INCLUDE_DIRS}")
message(STATUS "Use zlib library: ${ZLIB_LIBRARIES}")

macro(zlib_build)
  ExternalProject_Add(zlib
    PREFIX              ${CMAKE_BINARY_DIR}/ThirdParty/zlib
    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/ThirdParty/zlib
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX:STRING=${CMAKE_BINARY_DIR}/ThirdParty/zlib/install
      # -DBUILD_SHARED_LIBS:BOOL=OFF
  )
  set_target_properties(zlib PROPERTIES EXCLUDE_FROM_ALL ON)
  add_dependencies(build_bundled_libs zlib)
  add_dependencies(p4est_bundled_libs zlib)
endmacro()


if(ENABLE_BUNDLED_ZLIB)
  zlib_build()
endif()
