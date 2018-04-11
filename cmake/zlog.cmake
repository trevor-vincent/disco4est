include (ExternalProject)

#
# Bundled zlog paths.
#
set(ZLOG_BUNDLED_PREFIX "${PROJECT_BINARY_DIR}/third_party/zlog/install")
set(ZLOG_BUNDLED_LIB    "${ZLOG_BUNDLED_PREFIX}/lib/libzlog.a")

macro(zlog_use_bundled)
  set(ZLOG_PREFIX "${ZLOG_BUNDLED_PREFIX}")
  set(ZLOG_INCLUDE_DIRS "${ZLOG_BUNDLED_PREFIX}/include")
  set(ZLOG_LIBRARIES "${ZLOG_BUNDLED_LIB}")
  set(ENABLE_BUNDLED_ZLOG True)
endmacro()

macro(zlog_try_system)
  find_package(ZLOG REQUIRED)
  message(STATUS "Found a system-wide zlog.")
endmacro()

macro(zlog_try_prefix)
  find_path(ZLOG_INCLUDE_DIRS zlog.h ${ZLOG_PREFIX}/include NO_DEFAULT_PATH)
  find_library(ZLOG_LIB libzlog.a ${ZLOG_PREFIX}/lib NO_DEFAULT_PATH)

  if(ZLOG_INCLUDE_DIRS AND ZLOG_LIB)
    set(ZLOG_LIBRARIES ${ZLOG_LIB})
  else()
    message(FATAL_ERROR "Couldn't find zlog in '${ZLOG_PREFIX}'")
  endif()
endmacro()

if (ZLOG_PREFIX)
  zlog_try_prefix()
elseif (NOT ENABLE_BUNDLED_ZLOG)
  zlog_try_system()
else()
  zlog_use_bundled()
endif()

include_directories(${ZLOG_INCLUDE_DIRS})

message(STATUS "Use zlog includes: ${ZLOG_INCLUDE_DIRS}")
message(STATUS "Use zlog library: ${ZLOG_LIBRARIES}")

macro(zlog_build)
  ExternalProject_Add(zlog
    PREFIX              ${CMAKE_BINARY_DIR}/third_party/zlog
    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/third_party/zlog
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       cd ${CMAKE_SOURCE_DIR}/third_party/zlog && make PREFIX=${ZLOG_BUNDLED_PREFIX}
    INSTALL_COMMAND     cd ${CMAKE_SOURCE_DIR}/third_party/zlog && make PREFIX=${ZLOG_BUNDLED_PREFIX} install
  )
  set_target_properties(zlog PROPERTIES EXCLUDE_FROM_ALL ON)
  add_dependencies(build_bundled_libs zlog)
endmacro()


if(ENABLE_BUNDLED_ZLOG)
  zlog_build()
endif()
