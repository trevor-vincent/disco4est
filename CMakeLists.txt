cmake_minimum_required(VERSION 3.1)

include(ExternalProject)

project(d4est)

# Directory where customize cmake files reside
set(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/ThirdParty/cmake
  ${CMAKE_MODULE_PATH}
)

# find_package(BLAS)
# if(NOT BLAS_FOUND)
#   message("Cannot find BLAS, searching for OpenBlas")
#   find_package(OpenBLAS REQUIRED)
#   if(OpenBLAS_FOUND)
#     message("Found OpenBLAS")
#   endif()
# else()
#   find_package(LAPACK REQUIRED)
# endif()

# find_package (Threads)
find_package(MPI REQUIRED)
# find_package(ZLIB REQUIRED)
# find_package(X11 REQUIRED)
# find_package(SSL REQUIRED)

##
## Set build include path
##
include_directories(${CMAKE_SOURCE_DIR}/InputParser)
include_directories(${CMAKE_SOURCE_DIR}/EllipticSystem)
include_directories(${CMAKE_SOURCE_DIR}/Geometry)
include_directories(${CMAKE_SOURCE_DIR}/Solver)
include_directories(${CMAKE_SOURCE_DIR}/Estimators)
include_directories(${CMAKE_SOURCE_DIR}/dGMath)
include_directories(${CMAKE_SOURCE_DIR}/LinearAlgebra)
include_directories(${CMAKE_SOURCE_DIR}/GridFunctions)
include_directories(${CMAKE_SOURCE_DIR}/Problems)
include_directories(${CMAKE_SOURCE_DIR}/Utilities)
include_directories(${CMAKE_SOURCE_DIR}/pXest)
include_directories(${CMAKE_SOURCE_DIR}/hpAMR)
include_directories(${CMAKE_SOURCE_DIR}/ElementData)
include_directories(${CMAKE_SOURCE_DIR}/Flux)
include_directories(${CMAKE_SOURCE_DIR}/Operators)
include_directories(${CMAKE_SOURCE_DIR}/Visualization)
include_directories(${CMAKE_BINARY_DIR})
include_directories(${MPI_C_INCLUDE_PATH})
# include_directories(${OpenBLAS_INCLUDE_DIR})

##
## Add MPI compile flags for the whole project
##
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAG}")
SET (
  SRCS
  InputParser/ini.c
  Solver/krylov_petsc.c
  Solver/multigrid.c
  Solver/newton_petsc.c
  Solver/cg_eigs.c
  Solver/matrix_sym_tester.c
  EllipticSystem/problem_data.c
  Utilities/util.c
  dGMath/dgmath.c
  GridFunctions/grid_functions.c
  LinearAlgebra/linalg.c
  Estimators/bi_estimator.c
  Estimators/bi_estimator_flux_fcns.c
  Geometry/d4est_geometry_compact_sphere.c
  # Estimators/curved_bi_estimator.c
  # Estimators/curved_bi_estimator_flux_fcns.c
  Estimators/estimator_stats.c
  # hpAMR/curved_hp_amr.c
  # hpAMR/curved_hp_amr_smooth_pred.c
  hpAMR/hp_amr.c
  Operators/curved_poisson_operator.c
  Operators/poisson_operator.c
  Operators/poisson_debug_vecs.c
  Operators/curved_poisson_debug_vecs.c
  Operators/grad.c
  hpAMR/hp_amr_smooth_pred.c
  hpAMR/hp_amr_uniform.c
  ElementData/curved_element_data.c
  ElementData/element_data.c
  Flux/central_flux_vector_fcns.c
  Flux/compute_flux.c
  Flux/curved_compute_flux.c
  Flux/curved_Gauss_central_flux_vector_fcns.c
  Flux/curved_Gauss_sipg_flux_scalar_fcns.c
  Flux/dg_norm.c
  Flux/sipg_flux_scalar_fcns.c
  Flux/sipg_flux_vector_fcns.c
  Visualization/hacked_p4est_vtk.c
  )


# /home/tvincent/Dropbox/Research/Codes/disco4est/src_new_build/ThirdParty/petsc/install/lib/libfblas.a
##
## D8EST Library
##
add_library(d8est ${SRCS})
add_library(d4est ${SRCS})

add_dependencies(d8est build_bundled_libs)
add_dependencies(d4est build_bundled_libs)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")

add_custom_target(build_bundled_libs)
add_custom_target(p4est_bundled_libs)

include(zlib)
include_directories(${ZLIB_INCLUDE_DIRS})


include(petsc)
include(p4est)

subdirs(ThirdParty)

SET (
  LIBS
  gfortran
  # ${Threads_LIBRARIES}
  # ${ZLIB_LIBRARIES}
  ${MPI_C_LIBRARIES}
  ${PETSC_LIBRARIES}
  ${P4EST_LIBRARIES}
  ${MATH_LIB}
  ${CMAKE_BINARY_DIR}/ThirdParty/petsc/install/lib/libopenblas.a
  ${CMAKE_BINARY_DIR}/ThirdParty/zlib/install/lib/libz.a
  # ${CMAKE_BINARY_DIR}/ThirdParty/petsc/install/lib/libfblas.a
  # ${CMAKE_BINARY_DIR}/ThirdParty/petsc/install/lib/libflapack.a
  # ${LAPACK_LIBRARIES}
  # ${BLAS_LIBRARIES}

  # ${X11_LIBRARIES}
  # ${SSL_LIBRARIES}
  # X11
  # crypto
  # ssl
  )

# -llapack -lblas -lX11 -lssl -lcrypto -lm -L/usr/lib/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/5 -L/usr/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lmpi_f90 -lmpi_f77 -lgfortran -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.9 -lgfortran -lm -lquadmath -lm -lmpi_cxx -lstdc++ -L/usr/lib/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/5 -L/usr/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -ldl -lmpi -lhwloc -lgcc_s -lpthread -ldl 


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

# add_d4est_problem_3d(testcompile1 testcompile.c TestCompile)
# add_d4est_problem_2d(testcompile2 testcompile.c TestCompile)
add_d4est_problem_3d(constantdensitystar_driver constantdensitystar.c ConstantDensityStar)
add_d4est_problem_2d(stammhpmultigrid_driver stammhpmultigrid.c StammHPMultigrid)
add_d4est_problem_3d(twopuncturescube_driver twopuncturescube.c TwoPuncturesCube)
add_d4est_problem_3d(testcompactification_driver testcompactification.c TestCompactification)

# add_executable(testcompile1 driver.c ./Problems/TestCompile/testcompile.c)
# target_link_libraries(testcompile1 d8est)
# set_target_properties(testcompile1 PROPERTIES COMPILE_DEFINITIONS "DIM=3")

# add_executable(testcompile2 driver.c ./Problems/TestCompile/testcompile.c)
# target_link_libraries(testcompile2 d4est)
# set_target_properties(testcompile2 PROPERTIES COMPILE_DEFINITIONS "DIM=2")
# include("${CMAKE_BINARY_DIR}/ThirdParty/petsc/install/lib/petsc/conf/PETScBuildInternal.cmake")
