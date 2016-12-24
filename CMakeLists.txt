cmake_minimum_required(VERSION 3.1)

include(ExternalProject)

project(d4est)

# find_package(Git)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)
find_package(ZLIB REQUIRED)

# Directory where customize cmake files reside
set(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/ThirdParty/cmake
  ${CMAKE_MODULE_PATH}
)

##
## Set build include path
##
include_directories(${CMAKE_SOURCE_DIR}/InputParser)
include_directories(${CMAKE_SOURCE_DIR}/Problems)
include_directories(${CMAKE_SOURCE_DIR}/Utilities)
include_directories(${CMAKE_SOURCE_DIR}/pXest)
include_directories(${CMAKE_BINARY_DIR})
include_directories(${MPI_C_INCLUDE_PATH})

##
## Add MPI compile flags for the whole project
##
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAG}")
SET (
  SRCS
  InputParser/ini.c
  Utilities/util.c
  )

##
## D8EST Library
##
add_library(d8est ${SRCS})
add_library(d4est ${SRCS})

add_dependencies(d8est build_bundled_libs)
add_dependencies(d4est build_bundled_libs)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")

add_custom_target(build_bundled_libs)

include(petsc)
include(p4est)

subdirs(ThirdParty)

SET (
  LIBS
  ${ZLIB_LIBRARIES}
  ${MPI_C_LIBRARIES}
  ${PETSC_LIBRARIES}
  ${P4EST_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${MATH_LIB}
  )

find_library(MATH_LIB m)

target_link_libraries(d4est ${LIBS})
set_target_properties(d4est PROPERTIES COMPILE_DEFINITIONS "DIM=2")

target_link_libraries(d8est ${LIBS})
set_target_properties(d8est PROPERTIES COMPILE_DEFINITIONS "DIM=3")

function(add_d4est_problem_3d exec src parent_dir)
  add_executable(${exec} driver.c ${CMAKE_SOURCE_DIR}/Problems/${parent_dir}/${src})
  target_link_libraries(${exec} d8est)
  set_target_properties(${exec} PROPERTIES COMPILE_DEFINITIONS "DIM=3")
  set_target_properties(${exec}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${parent_dir}"
    )
  add_custom_command(
    TARGET ${exec} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/Problems/${parent_dir}/options.input
    ${CMAKE_CURRENT_BINARY_DIR}/${parent_dir}/options.input)  
endfunction()

function(add_d4est_problem_2d exec src parent_dir)
  add_executable(${exec} driver.c ${CMAKE_SOURCE_DIR}/Problems/${parent_dir}/${src})
  target_link_libraries(${exec} d4est)
  set_target_properties(${exec} PROPERTIES COMPILE_DEFINITIONS "DIM=2")
  set_target_properties(${exec}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${parent_dir}"
    )
  add_custom_command(
    TARGET ${exec} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/Problems/${parent_dir}/options.input
    ${CMAKE_CURRENT_BINARY_DIR}/${parent_dir}/options.input)  
endfunction()

add_d4est_problem_3d(testcompile1 testcompile.c TestCompile)
add_d4est_problem_2d(testcompile2 testcompile.c TestCompile)

# add_executable(testcompile1 driver.c ./Problems/TestCompile/testcompile.c)
# target_link_libraries(testcompile1 d8est)
# set_target_properties(testcompile1 PROPERTIES COMPILE_DEFINITIONS "DIM=3")

# add_executable(testcompile2 driver.c ./Problems/TestCompile/testcompile.c)
# target_link_libraries(testcompile2 d4est)
# set_target_properties(testcompile2 PROPERTIES COMPILE_DEFINITIONS "DIM=2")
