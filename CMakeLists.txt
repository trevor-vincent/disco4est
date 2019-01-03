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
message("mpi version = ${MPI_VERSION}")
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
  src/dGMath/d4est_laplacian.c
  src/dGMath/d4est_laplacian_flux.c
  src/dGMath/d4est_laplacian_flux_sipg.c
  src/dGMath/d4est_hessian.c
  src/Mesh/d4est_mortars_with_opt.c
  # src/Mesh/d4est_mortars_aux.c
  src/dGMath/d4est_laplacian_with_opt.c
  src/dGMath/d4est_laplacian_with_opt_flux.c
  src/dGMath/d4est_laplacian_with_opt_flux_sipg.c
  src/Geometry/d4est_geometry_brick.c
  src/Geometry/d4est_geometry.c
  src/Estimators/d4est_estimator_stats.c
  src/Estimators/d4est_estimator_bi.c
  src/Estimators/d4est_estimator_bi_new.c

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
  src/Mesh/d4est_ghost_data_ext.c
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
  src/IO/p4est_vtk_ext.c
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
  src/Solver/d4est_solver_multigrid_smoother_schwarz.c
  src/Solver/d4est_solver_multigrid_smoother_krylov_petsc.c
  src/Solver/d4est_solver_cg_eigs.c
  src/Solver/d4est_solver_schwarz_metadata.c
  src/Solver/d4est_solver_schwarz_operators.c
  src/Solver/d4est_solver_schwarz_helpers.c
  src/Solver/d4est_solver_schwarz_laplacian_ext.c
  src/Solver/d4est_solver_schwarz_transfer_ghost_data.c
  src/Solver/d4est_solver_schwarz_subdomain_solver_cg.c
  src/Solver/d4est_solver_schwarz_subdomain_solver.c
  src/Solver/d4est_solver_schwarz_apply_lhs.c
  src/Solver/d4est_solver_schwarz_geometric_data.c
  src/Solver/d4est_solver_schwarz.c
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
  ${BLAS_LIBRARIES}
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

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
add_definitions("-DGIT_BRANCH=${GIT_BRANCH}")
configure_file(
  ${CMAKE_SOURCE_DIR}/version.h.in
  ${CMAKE_BINARY_DIR}/generated/version.h
  )

include_directories(${CMAKE_BINARY_DIR}/generated)

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
set_property(TARGET d4est_test PROPERTY C_STANDARD 99)

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
  # target_compile_definitions(${exec} PRIVATE --fprofile-arcs --ftest-coverage)
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

  # SET(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")
  # SET(GCC_COVERAGE_LINK_FLAGS    "-lgcov")
  target_link_libraries(${exec} d4est_test)
  add_test(NAME ${test_name} COMMAND ${CMAKE_BINARY_DIR}/Tests/${exec} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Tests/)
  target_compile_definitions(${exec} PRIVATE DIM=2 D4EST_TEST=1)
  # target_compile_definitions(${exec} PRIVATE --fprofile-arcs --ftest-coverage)


  # SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
  # SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")

  #hopefully allows for this optimization
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif()
  
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

# if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
#     find_program( CODECOV_GCOV gcov )
#     add_definitions( -fprofile-arcs -ftest-coverage )
#     link_libraries( gcov )
# endif ( CMAKE_BUILD_TYPE STREQUAL "Debug" )

##
## Tests
##
enable_testing()
set(CTEST_MEMORYCHECK_COMMAND /usr/bin/valgrind)
set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full --suppressions=${PROJECT_SOURCE_DIR}/suppressions.txt" )

# add_d4est_unit_test_3d(d4est_test_vtk d4est_test_vtk.c Unit d4est_test_vtk)
# add_d4est_unit_test_3d(d4est_test_cg_eigs d4est_test_cg_eigs.c Unit d4est_test_cg_eigs)
# add_d4est_unit_test_3d(d4est_test_kron d4est_test_kron.c Unit d4est_test_kron)
# add_d4est_unit_test_3d(d4est_test_linalg d4est_test_linalg.c Unit d4est_test_linalg)
# add_d4est_unit_test_3d(d4est_test_ghost_data d4est_test_ghost_data.c Unit d4est_test_ghost_data)
# add_d4est_unit_test_3d(d4est_test_schwarz_operators d4est_test_schwarz_operators.c Unit d4est_test_schwarz_operators)
# add_d4est_unit_test_3d(d4est_test_schwarz_cubed_sphere_new_3d d4est_test_schwarz_cubed_sphere_new.c Unit d4est_test_schwarz_cubed_sphere_new_3d)

add_d4est_unit_test_3d(d4est_test_schwarz_cubed_sphere_hanging d4est_test_schwarz_cubed_sphere_hanging.c Unit d4est_test_schwarz_cubed_sphere_hanging)
add_d4est_unit_test_3d(d4est_test_schwarz_slowerAu_cubed_sphere_hanging d4est_test_schwarz_slowerAu_cubed_sphere_hanging.c Unit d4est_test_schwarz_slowerAu_cubed_sphere_hanging)



add_d4est_unit_test_2d(d4est_test_schwarz_disk_hanging d4est_test_schwarz_disk_hanging.c Unit d4est_test_schwarz_disk_hanging)
add_d4est_unit_test_2d(d4est_test_schwarz_slowerAu_disk_hanging d4est_test_schwarz_slowerAu_disk_hanging.c Unit d4est_test_schwarz_slowerAu_disk_hanging)



# add_d4est_unit_test_2d(d4est_test_schwarz_hanging d4est_test_schwarz_hanging.c Unit d4est_test_schwarz_hanging)
# add_d4est_unit_test_2d(d4est_test_schwarz_slowerAu_hanging d4est_test_schwarz_slowerAu_hanging.c Unit d4est_test_schwarz_slowerAu_hanging)

# add_d4est_unit_test_3d(d4est_test_laplacian_symmetry d4est_test_laplacian_symmetry.c Unit d4est_test_laplacian_symmetry)
# add_d4est_unit_test_3d(d4est_test_schwarz_cubed_sphere d4est_test_schwarz_cubed_sphere.c Unit d4est_test_schwarz_cubed_sphere)
# add_d4est_unit_test_3d(d4est_test_laplacian_consistency d4est_test_laplacian_consistency.c Unit d4est_test_laplacian_consistency)
# add_d4est_unit_test_3d(d4est_test_laplacian_speedup d4est_test_laplacian_speedup.c Unit d4est_test_laplacian_speedup)
# add_d4est_bash_regression_test(d4est_test_stamm.sh Regression d4est_test_stamm)
# add_d4est_bash_regression_test(d4est_test_cds.sh Regression d4est_test_constant_density_star)
# add_d4est_bash_regression_test(d4est_test_poisson_lorentzian.sh Regression d4est_test_poisson_lorentzian)
# add_d4est_bash_regression_test(d4est_test_poisson_sinx.sh Regression d4est_test_poisson_sinx)
# add_d4est_bash_regression_test(d4est_test_twopunctures.sh Regression d4est_test_twopunctures)
# add_d4est_bash_regression_test(d4est_test_twopunctures_2.sh Regression d4est_test_twopunctures_2)
# add_d4est_bash_regression_test(d4est_test_mpi.sh Regression d4est_test_mpi)
# add_d4est_problem_2d(stamm_driver stamm_multigrid_pc.c Stamm)
# add_d4est_problem_2d(stamm_anares_driver stamm_anares.c Stamm)
# add_d4est_problem_3d(poisson_sinx_uniform_driver poisson_sinx_uniform.c Poisson)
# add_d4est_problem_3d(poisson_sinx_multigrid_driver poisson_sinx_multigrid.c Poisson)
# add_d4est_problem_3d(poisson_lorentzian_multigrid_driver poisson_lorentzian_multigrid.c Poisson)
# add_d4est_problem_3d(poisson_lorentzian_multigrid_schwarz_driver poisson_lorentzian_multigrid_schwarz.c Poisson)
# add_d4est_problem_3d(poisson_lorentzian_amr_multigrid_driver poisson_lorentzian_amr_multigrid.c Poisson)
# add_d4est_problem_3d(poisson_lorentzian_amr_multigrid_numres_driver poisson_lorentzian_amr_multigrid_numres.c Poisson)
# add_d4est_problem_3d(poisson_lorentzian_amr_multigrid_anares_driver poisson_lorentzian_amr_multigrid_anares.c Poisson)
# add_d4est_problem_3d(constant_density_star_driver constant_density_star_mgpc_newton_petsc.c ConstantDensityStar)
add_d4est_problem_3d(constant_density_star_anares_driver cds_anares.c ConstantDensityStar)
# add_d4est_problem_3d(two_punctures_robin_diff_geom_estimator_global_driver two_punctures_cactus_diff_geom_estimator_global.c TwoPunctures)
# add_d4est_problem_3d(two_punctures_robin_diff_geom_estimator_global_13tree_driver two_punctures_cactus_diff_geom_estimator_global_13tree.c TwoPunctures)
# add_d4est_problem_3d(two_punctures_robin_diff_geom_estimator_global_13tree_with_opt_driver two_punctures_cactus_diff_geom_estimator_global_13tree_with_opt.c TwoPunctures)
# add_d4est_problem_3d(two_punctures_anares_driver two_punctures_anares.c TwoPunctures)

# add_d4est_problem_3d(two_punctures_cactus_13tree_with_opt_puncture_finder_driver two_punctures_cactus_13tree_with_opt_puncture_finder.c TwoPunctures)
