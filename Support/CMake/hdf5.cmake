include (ExternalProject)

#
# Bundled hdf5 paths.
#
set(HDF5_BUNDLED_PREFIX "${PROJECT_BINARY_DIR}/ThirdParty/hdf5/install")
set(HDF5_BUNDLED_LIB    "${HDF5_BUNDLED_PREFIX}/lib/libhdf5-static.a")

macro(hdf5_use_bundled)
  set(HDF5_PREFIX "${HDF5_BUNDLED_PREFIX}")
  set(HDF5_INCLUDE_DIRS "${HDF5_BUNDLED_PREFIX}/include")
  set(HDF5_LIBRARIES "${HDF5_BUNDLED_LIB}")
  set(ENABLE_BUNDLED_HDF5 True)
endmacro()

macro(hdf5_try_system)
  find_package(HDF5 REQUIRED)
  message(STATUS "Found a system-wide hdf5.")
endmacro()

macro(hdf5_try_prefix)
  find_path(HDF5_INCLUDE_DIRS hdf5.h ${HDF5_PREFIX}/include NO_DEFAULT_PATH)
  find_library(HDF5_LIB libhdf5-static.a ${HDF5_PREFIX}/lib NO_DEFAULT_PATH)

  if(HDF5_INCLUDE_DIRS AND HDF5_LIB)
    set(HDF5_LIBRARIES ${HDF5_LIB})
  else()
    message(FATAL_ERROR "Couldn't find hdf5 in '${HDF5_PREFIX}'")
  endif()
endmacro()

if (HDF5_PREFIX)
  hdf5_try_prefix()
elseif (NOT ENABLE_BUNDLED_HDF5)
  hdf5_try_system()
else()
  hdf5_use_bundled()
endif()

include_directories(${HDF5_INCLUDE_DIRS})

message(STATUS "Use hdf5 includes: ${HDF5_INCLUDE_DIRS}")
message(STATUS "Use hdf5 library: ${HDF5_LIBRARIES}")

macro(hdf5_build)
  ExternalProject_Add(hdf5
    URL                 https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.1.tar.gz
    PREFIX              ${CMAKE_BINARY_DIR}/ThirdParty/hdf5
    CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:STRING=${CMAKE_BINARY_DIR}/ThirdParty/hdf5/install
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DCMAKE_CXX_FLAGS:STRING=-Wunused-but-set-variable -Wunused-variable -Wincompatible-pointer-types -Wformat-nonliteral -Wstrict-overflow -Wconversion -Wunused-but-set-variable -Wunsafe-loop-optimizations -Wunused-parameter -Wlarger-than= -Wdiscarded-qualifiers -Wfloat-equal -Wmaybe-uninitialized -Wcomment
    -DCMAKE_C_FLAGS:STRING=-Wunused-but-set-variable -Wunused-variable -Wincompatible-pointer-types -Wformat-nonliteral -Wstrict-overflow -Wconversion -Wunused-but-set-variable -Wunsafe-loop-optimizations -Wunused-parameter -Wlarger-than= -Wdiscarded-qualifiers -Wfloat-equal -Wmaybe-uninitialized -Wcomment
    )
    add_dependencies(build_bundled_libs hdf5)
endmacro()

if(ENABLE_BUNDLED_HDF5)
  hdf5_build()
endif()
