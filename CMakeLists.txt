cmake_minimum_required(VERSION 3.0)
# function(message)
#   list(GET ARGV 0 MessageType)
#   if(MessageType STREQUAL FATAL_ERROR OR
#      MessageType STREQUAL SEND_ERROR OR
#      MessageType STREQUAL WARNING OR
#      MessageType STREQUAL AUTHOR_WARNING)
#     list(REMOVE_AT ARGV 0)
#     _message(${MessageType} "${ARGV}")
#   endif()
# endfunction()


include(ExternalProject)

project(d4est)

# Directory where customize cmake files reside
set(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake
  ${CMAKE_MODULE_PATH}
)
include(machine)

find_package(MPI REQUIRED)

##
## Set build include path
##
include_directories(${CMAKE_SOURCE_DIR}/src/EllipticSystem)
include_directories(${CMAKE_SOURCE_DIR}/src/Geometry)
include_directories(${CMAKE_SOURCE_DIR}/src/Solver)
include_directories(${CMAKE_SOURCE_DIR}/src/Estimators)
include_directories(${CMAKE_SOURCE_DIR}/src/dGMath)
include_directories(${CMAKE_SOURCE_DIR}/src/LinearAlgebra)
include_directories(${CMAKE_SOURCE_DIR}/src/GridFunctions)
include_directories(${CMAKE_SOURCE_DIR}/src/Problems)
include_directories(${CMAKE_SOURCE_DIR}/src/Utilities)
include_directories(${CMAKE_SOURCE_DIR}/src/pXest)
include_directories(${CMAKE_SOURCE_DIR}/src/hpAMR)
include_directories(${CMAKE_SOURCE_DIR}/src/ElementData)
include_directories(${CMAKE_SOURCE_DIR}/src/Mesh)
include_directories(${CMAKE_SOURCE_DIR}/src/Operators)
include_directories(${CMAKE_SOURCE_DIR}/src/IO)
include_directories(${CMAKE_SOURCE_DIR}/src/Quadrature)
include_directories(${CMAKE_SOURCE_DIR}/src/Utilities)
include_directories(${CMAKE_SOURCE_DIR}/src/Kron)
include_directories(${CMAKE_SOURCE_DIR}/src/Databox)
include_directories(${CMAKE_SOURCE_DIR}/src/Dictionary)
include_directories(${CMAKE_SOURCE_DIR}/third_party/mpsort)
include_directories(${CMAKE_SOURCE_DIR}/third_party/PPF)
include_directories(${CMAKE_BINARY_DIR})
include_directories(${MPI_C_INCLUDE_PATH})

##
## Add MPI compile flags for the whole project
##
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAG}")

SET (
  SRCS_3D_ONLY
  src/Geometry/d4est_geometry_cubed_sphere.c
  src/Geometry/d4est_connectivity_cubed_sphere.c
  src/Geometry/d4est_geometry_general_wedge.c
  src/Geometry/d4est_geometry_hole_in_a_box.c
)

SET (
  SRCS_2D_ONLY
  src/Geometry/d4est_geometry_disk.c
)

SET (
  SRCS
  src/Utilities/ini.c
  src/Utilities/d4est_util.c
  src/LinearAlgebra/d4est_linalg.c
  src/EllipticSystem/d4est_elliptic_eqns.c
  src/EllipticSystem/d4est_elliptic_data.c
  src/Kron/d4est_kron.c
  src/dGMath/d4est_operators.c
  src/dGMath/d4est_gradient.c
  src/dGMath/d4est_ip_energy_norm.c
  src/dGMath/d4est_lgl.c
  src/dGMath/d4est_reference.c
  src/dGMath/d4est_poisson.c
  src/dGMath/d4est_poisson_flux.c
  src/dGMath/d4est_poisson_flux_sipg.c
  src/Geometry/d4est_geometry_brick.c
  src/Geometry/d4est_geometry.c
  src/Estimators/d4est_estimator_stats.c
  src/Estimators/d4est_estimator_bi.c
  # src/Estimators/d4est_estimator_residual.c
  src/hpAMR/d4est_amr.c
  src/hpAMR/d4est_amr_uniform.c
  src/hpAMR/d4est_amr_random.c
  src/hpAMR/d4est_amr_smooth_pred.c
  src/Mesh/d4est_element_data.c
  src/Mesh/d4est_mesh.c
  src/Mesh/d4est_xyz_functions.c
  src/Mesh/d4est_mortars.c
  src/Mesh/d4est_mortars_aux.c
  src/Mesh/d4est_ghost_data.c
  src/Mesh/d4est_ghost.c
  src/Databox/d4est_databox.c
  src/Dictionary/d4est_dictionary.c
  src/Dictionary/d4est_critbit.c
  src/Quadrature/d4est_quadrature.c
  src/Quadrature/d4est_quadrature_legendre.c
  src/Quadrature/d4est_quadrature_lobatto.c
  src/IO/d4est_h5.c
  src/IO/d4est_checkpoint.c
  src/IO/d4est_vtk.c
  src/IO/d4est_norms.c
  src/Solver/d4est_solver_cg.c
  src/Solver/d4est_solver_fcg.c
  # src/Solver/d4est_solver_newton.c
  src/Solver/d4est_solver_krylov_petsc.c
  src/Solver/d4est_solver_newton_petsc.c
  src/Solver/d4est_krylov_pc_multigrid.c
  src/Solver/d4est_solver_multigrid.c
  src/Solver/d4est_solver_multigrid_bottom_solver_cg.c
  src/Solver/d4est_solver_multigrid_bottom_solver_cheby.c
  src/Solver/d4est_solver_multigrid_bottom_solver_krylov_petsc.c
  src/Solver/d4est_solver_multigrid_element_data_updater.c
  src/Solver/d4est_solver_multigrid_logger_residual.c
  src/Solver/d4est_solver_multigrid_profiler_basic.c
  src/Solver/d4est_solver_multigrid_mesh_analyzer.c
  src/Solver/d4est_solver_multigrid_matrix_operator.c
  src/Solver/d4est_solver_multigrid_smoother_cheby.c
  src/Solver/d4est_solver_multigrid_smoother_krylov_petsc.c
  src/Solver/d4est_solver_cg_eigs.c
  # src/Solver/d4est_solver_jacobian_tester.c
  src/Solver/d4est_solver_full_matrix.c
  src/Solver/d4est_solver_matrix_symmetry.c
  src/Solver/d4est_power_method.c
  third_party/mpsort/radixsort.c
  third_party/mpsort/mpsort-mpi.c
  third_party/PPF/PPF_Print.c
  third_party/PPF/PPF_Printf.c
  third_party/PPF/ClosedSet.c
  )

add_library(d8est ${SRCS} ${SRCS_3D_ONLY})
add_library(d8est_test ${SRCS} ${SRCS_3D_ONLY})
add_library(d4est ${SRCS} ${SRCS_2D_ONLY})
add_library(d4est_test ${SRCS} ${SRCS_2D_ONLY})


add_dependencies(d8est build_bundled_libs)
add_dependencies(d8est_test build_bundled_libs)
add_dependencies(d4est_test build_bundled_libs)
add_dependencies(d4est build_bundled_libs)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
add_custom_target(build_bundled_libs)

include(openblas)
include(zlib)
include(p4est)
include(petsc)
include(zlog)
include(hdf5)

subdirs(third_party)

SET (
  LIBS
  ${MPI_C_LIBRARIES}
  ${P4EST_LIBRARIES}
  ${PETSC_LIBRARIES}
  m
  dl
  # ${CMAKE_BINARY_DIR}/ThirdParty/zlib/install/lib/libz.a
  ${ZLIB_LIBRARIES}
  ${OPENBLAS_LIBRARIES}
  pthread
  gfortran
  ${HDF5_LIBRARIES}
  ${ZLOG_LIBRARIES}
  )

if(USE_TCMALLOC)
  include(tcmalloc)
  SET (
    LIBS
    ${LIBS}
    ${TCMALLOC_LIBRARIES}
    stdc++
    )
  find_package(Libunwind)
  if(LIBUNWIND_FOUND)
    SET (
      LIBS
      ${LIBS}
      unwind
      )
  endif()
endif()

message("LIBS=${LIBS}")
target_link_libraries(d4est ${LIBS})
set_property(TARGET d4est APPEND_STRING PROPERTY
  COMPILE_FLAGS " -DDIM=2")

target_link_libraries(d8est ${LIBS})
set_property(TARGET d8est APPEND_STRING PROPERTY
  COMPILE_FLAGS " -DDIM=3")

target_link_libraries(d8est_test ${LIBS})
set_property(TARGET d8est_test APPEND_STRING PROPERTY
  COMPILE_FLAGS "-DDIM=3 -DD4EST_TEST=1")

target_link_libraries(d4est_test ${LIBS})
set_property(TARGET d4est_test APPEND_STRING PROPERTY
  COMPILE_FLAGS " -DDIM=2 -DD4EST_TEST=1")


set_property(TARGET d4est PROPERTY C_STANDARD 99)
set_property(TARGET d8est PROPERTY C_STANDARD 99)
set_property(TARGET d8est_test PROPERTY C_STANDARD 99)

function(add_d4est_problem_3d exec src parent_dir)
  add_executable(${exec} src/driver.c ${CMAKE_SOURCE_DIR}/src/Problems/${parent_dir}/${src})
  target_link_libraries(${exec} d8est)
  set_target_properties(${exec} PROPERTIES COMPILE_DEFINITIONS "DIM=3")
  target_compile_definitions(${exec} PRIVATE D4EST_PROBLEM_NAME="${parent_dir}")
  set_target_properties(${exec}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${parent_dir}"
    )
  add_custom_command(
    TARGET ${exec} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/src/Problems/${parent_dir}/*.input
    ${CMAKE_CURRENT_BINARY_DIR}/${parent_dir}/)
  add_custom_command(
    TARGET ${exec} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/src/IO/logging.conf
    ${CMAKE_CURRENT_BINARY_DIR}/${parent_dir}/logging.conf)
endfunction()


function(add_d4est_problem_2d exec src parent_dir)
  add_executable(${exec} src/driver.c ${CMAKE_SOURCE_DIR}/src/Problems/${parent_dir}/${src})
  target_link_libraries(${exec} d4est)
  set_target_properties(${exec} PROPERTIES COMPILE_DEFINITIONS "DIM=2")
  set_target_properties(${exec}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${parent_dir}"
    )
  add_custom_command(
    TARGET ${exec} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/src/Problems/${parent_dir}/*.input
    ${CMAKE_CURRENT_BINARY_DIR}/${parent_dir}/)
  add_custom_command(
    TARGET ${exec} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/src/IO/logging.conf
    ${CMAKE_CURRENT_BINARY_DIR}/${parent_dir}/logging.conf)
endfunction()


function(add_d4est_unit_test_3d exec src parent_dir test_name)
  add_executable(${exec} ${CMAKE_SOURCE_DIR}/src/Tests/${parent_dir}/${src})
  target_link_libraries(${exec} d8est_test)
  add_test(NAME ${test_name} COMMAND ${CMAKE_BINARY_DIR}/Tests/${exec} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Tests/)
  # target_compile_definitions(${exec} PRIVATE DIM=3)
  target_compile_definitions(${exec} PRIVATE DIM=3 D4EST_TEST=1)
  set_target_properties(${exec}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Tests/"
    )

  
  add_custom_command(
    TARGET ${exec} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/src/Tests/${parent_dir}/*.input
    "${CMAKE_BINARY_DIR}/Tests/")
  add_custom_command(
    TARGET ${exec} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/src/IO/logging.conf
    ${CMAKE_CURRENT_BINARY_DIR}/Tests/logging.conf)
endfunction()


function(add_d4est_unit_test_2d exec src parent_dir test_name)
  add_executable(${exec} ${CMAKE_SOURCE_DIR}/src/Tests/${parent_dir}/${src})
  target_link_libraries(${exec} d4est_test)
  add_test(NAME ${test_name} COMMAND ${CMAKE_BINARY_DIR}/Tests/${exec} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Tests/)
  target_compile_definitions(${exec} PRIVATE DIM=2 D4EST_TEST=1)
  set_target_properties(${exec}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Tests/"
    )
  # target_compile_definitions(${exec} PRIVATE D4EST_TEST=1)  
  add_custom_command(
    TARGET ${exec} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/src/Tests/${parent_dir}/*.input
    "${CMAKE_BINARY_DIR}/Tests/")
  add_custom_command(
    TARGET ${exec} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/src/IO/logging.conf
    ${CMAKE_CURRENT_BINARY_DIR}/Tests/logging.conf)
endfunction()

function(add_d4est_bash_regression_test exec parent_dir test_name)
    add_custom_command(
      TARGET d4est POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_SOURCE_DIR}/src/Tests/${parent_dir}/*.sh
      "${CMAKE_BINARY_DIR}/Tests/")

    add_test(NAME ${test_name} COMMAND bash ${CMAKE_BINARY_DIR}/Tests/${exec} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Tests/)  
endfunction()


##
## Tests
##
enable_testing()
set(CTEST_MEMORYCHECK_COMMAND /usr/bin/valgrind)
set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full --suppressions=${PROJECT_SOURCE_DIR}/suppressions.txt" )

add_d4est_unit_test_3d(test_d4est_vtk test_d4est_vtk.c Unit test_d4est_vtk)
# add_d4est_unit_test_3d(test_d4est_iterate_2 test_d4est_iterate_2.c Unit test_d4est_iterate_2)
# add_d4est_unit_test_3d(test_d4est_iterate test_d4est_iterate.c Unit test_d4est_iterate)
add_d4est_unit_test_2d(test_d4est_iterate_2d test_d4est_iterate.c Unit test_d4est_iterate_2d)

add_d4est_unit_test_3d(test_d4est_operators test_d4est_operators.c Unit test_d4est_operators)
add_d4est_unit_test_3d(test_d4est_kron test_d4est_kron.c Unit test_d4est_kron)
add_d4est_unit_test_3d(test_d4est_linalg test_d4est_linalg.c Unit test_d4est_linalg)
add_d4est_unit_test_3d(test_d4est_ghost_data test_d4est_ghost_data.c Unit test_d4est_ghost_data)

add_d4est_bash_regression_test(d4est_test_stamm.sh Regression d4est_test_stamm)
add_d4est_bash_regression_test(d4est_test_cds.sh Regression d4est_test_constant_density_star)
add_d4est_bash_regression_test(d4est_test_poisson_sinx.sh Regression d4est_test_poisson_sinx)
# add_d4est_bash_regression_test(d4est_test_poisson_lorentzian.sh Regression d4est_test_poisson_lorentzian)
add_d4est_bash_regression_test(d4est_test_twopunctures.sh Regression d4est_test_twopunctures)
add_d4est_bash_regression_test(d4est_test_twopunctures_2.sh Regression d4est_test_twopunctures_2)
add_d4est_bash_regression_test(d4est_test_mpi.sh Regression d4est_test_mpi)

add_d4est_problem_2d(stamm_driver stamm_multigrid_pc.c Stamm)
add_d4est_problem_3d(poisson_sinx_uniform_driver poisson_sinx_uniform.c Poisson)
# add_d4est_problem_3d(poisson_sinx_uniform_debug_driver poisson_sinx_uniform_debug.c Poisson)
add_d4est_problem_3d(poisson_sinx_multigrid_driver poisson_sinx_multigrid.c Poisson)
add_d4est_problem_3d(poisson_lorentzian_multigrid_driver poisson_lorentzian_multigrid.c Poisson)
# add_d4est_problem_3d(poisson_lorentzian_cg_driver poisson_lorentzian_cg.c Poisson)
# add_d4est_problem_3d(poisson_lorentzian_eigs_driver poisson_lorentzian_eigs.c Poisson)

add_d4est_problem_3d(poisson_lorentzian_amr_multigrid_driver poisson_lorentzian_amr_multigrid.c Poisson)
# add_d4est_problem_3d(poisson_1orton_multigrid_driver poisson_1orton_multigrid.c Poisson)
add_d4est_problem_3d(constant_density_star_driver constant_density_star_mgpc_newton_petsc.c ConstantDensityStar)
add_d4est_problem_3d(two_punctures_robin_diff_geom_estimator_global_driver two_punctures_cactus_diff_geom_estimator_global.c TwoPunctures)
add_d4est_problem_3d(two_punctures_robin_diff_geom_estimator_global_13tree_driver two_punctures_cactus_diff_geom_estimator_global_13tree.c TwoPunctures)
add_d4est_problem_3d(two_punctures_estimator_test_driver two_punctures_estimator_test.c TwoPunctures)
