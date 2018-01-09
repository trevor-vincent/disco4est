cmake_minimum_required(VERSION 2.8)

include(ExternalProject)

project(d4est)

# Directory where customize cmake files reside
set(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/ThirdParty/cmake
  ${CMAKE_MODULE_PATH}
)

if(UBUNTU)
  set(USE_OPENBLAS ON)
  set(OPENBLAS_DIR "/home/tvincent/Dropbox/Research/Codes/OpenBLAS/")
  set(OPENBLAS_ROOT "/home/tvincent/Dropbox/Research/Codes/OpenBLAS/")
  set(OPENBLAS_LIBRARIES
    "/home/tvincent/Dropbox/Research/Codes/OpenBLAS/libopenblas.a"
    )
  set(OPENBLAS_INCLUDE_DIR "/home/tvincent/Dropbox/Research/Codes/OpenBLAS/lapack-netlib/LAPACKE/include/")  
  set(USE_TCMALLOC ON)
  set(TCMALLOC_DIR "/home/tvincent/Dropbox/Research/Codes/gperftools/gperftools_install/")
  set(TCMALLOC_ROOT "/home/tvincent/Dropbox/Research/Codes/gperftools/gperftools_install/")
  set(TCMALLOC_LIBRARIES "/home/tvincent/Dropbox/Research/Codes/gperftools/gperftools_install/lib/libtcmalloc.a")
  set(TCMALLOC_INCLUDE_DIR "/home/tvincent/Dropbox/Research/Codes/gperftools/gperftools_install/include/gperftools/")

  set(USE_HDF5 ON)
  set(HDF5_LIBRARIES "/home/tvincent/Dropbox/Research/Codes/hdf5-1.10.1/install/lib/libhdf5-static.a")
  set(HDF5_INCLUDE_DIR "/home/tvincent/Dropbox/Research/Codes/hdf5-1.10.1/install/include/")
endif()

if(SCINET)
  set(USE_OPENBLAS ON)
  set(OPENBLAS_DIR "/scinet/gpc/Libraries/OpenBLAS/singlethreaded/")
  set(OPENBLAS_ROOT "/scinet/gpc/Libraries/OpenBLAS/singlethreaded/")
  set(OPENBLAS_LIBRARIES
    "/scinet/gpc/Libraries/OpenBLAS/singlethreaded/lib/libopenblas.a"
    )
  set(OPENBLAS_INCLUDE_DIR "/scinet/gpc/Libraries/OpenBLAS/singlethreaded/include")
  set(USE_HDF5 ON)
  set(HDF5_LIBRARIES "/scratch/p/pfeiffer/tvincent/hdf5-1.10.1/install/lib/libhdf5-static.a")
  set(HDF5_INCLUDE_DIR "/scratch/p/pfeiffer/tvincent/hdf5-1.10.1/install/include/") 
endif()

if(GRAHAM OR CEDAR)
  set(USE_OPENBLAS ON)
  set(OPENBLAS_DIR "/scratch/tvincent/OpenBLAS/")
  set(OPENBLAS_ROOT "/scratch/tvincent/OpenBLAS/")
  set(OPENBLAS_LIBRARIES
    "/scratch/tvincent/OpenBLAS/libopenblas.a"
    )
  set(OPENBLAS_INCLUDE_DIR "/scratch/tvincent/OpenBLAS/lapack-netlib/LAPACKE/include/")
  set(USE_HDF5 ON)
  set(HDF5_LIBRARIES "/scratch/tvincent/hdf5-1.10.1/install/lib/libhdf5-static.a")
  set(HDF5_INCLUDE_DIR "/scratch/tvincent/hdf5-1.10.1/install/include/")
endif()


if(CITA)
  set(USE_OPENBLAS ON)
  set(OPENBLAS_DIR "/opt/openblas/0.2.19/")
  set(OPENBLAS_ROOT "/opt/openblas/0.2.19/")
  set(OPENBLAS_LIBRARIES
    "/opt/openblas/0.2.19/lib/libopenblas.a"
    )
  set(OPENBLAS_INCLUDE_DIR "/opt/openblas/0.2.19/include")
endif()

if(USE_OPENBLAS)
  include_directories(${OPENBLAS_INCLUDE_DIR})
endif()

if(USE_HDF5)
  include_directories(${HDF5_INCLUDE_DIR})
endif()

find_package(MPI REQUIRED)

##
## Set build include path
##
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
include_directories(${CMAKE_SOURCE_DIR}/Mesh)
include_directories(${CMAKE_SOURCE_DIR}/Operators)
include_directories(${CMAKE_SOURCE_DIR}/IO)
include_directories(${CMAKE_SOURCE_DIR}/Quadrature)
include_directories(${CMAKE_SOURCE_DIR}/Utilities)
include_directories(${CMAKE_SOURCE_DIR}/Kron)
include_directories(${CMAKE_BINARY_DIR})
include_directories(${MPI_C_INCLUDE_PATH})

##
## Add MPI compile flags for the whole project
##
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAG}")

SET (
  SRCS_3D_ONLY
  Geometry/d4est_geometry_cubed_sphere.c
  Geometry/d4est_connectivity_cubed_sphere.c
  Geometry/d4est_geometry_general_wedge.c
  Geometry/d4est_geometry_hole_in_a_box.c
  )

SET (
  SRCS_2D_ONLY
  Geometry/d4est_geometry_disk.c
  )

SET (
  SRCS
  Utilities/ini.c
  Utilities/d4est_util.c
  LinearAlgebra/d4est_linalg.c
  EllipticSystem/d4est_elliptic_eqns.c
  EllipticSystem/d4est_elliptic_data.c
  Kron/d4est_kron.c
  dGMath/d4est_operators.c
  dGMath/d4est_gradient.c
  dGMath/d4est_ip_energy_norm.c
  dGMath/d4est_lgl.c
  dGMath/d4est_reference.c
  dGMath/d4est_poisson.c
  dGMath/d4est_poisson_flux.c
  dGMath/d4est_poisson_flux_sipg.c
  Geometry/d4est_geometry_brick.c
  Geometry/d4est_geometry.c
  Estimators/d4est_estimator_stats.c
  Estimators/d4est_estimator_bi.c
  Estimators/d4est_estimator_residual.c
  hpAMR/d4est_amr.c
  hpAMR/d4est_amr_uniform.c
  hpAMR/d4est_amr_random.c
  hpAMR/d4est_amr_smooth_pred.c
  Mesh/d4est_element_data.c
  Mesh/d4est_mesh.c
  Mesh/d4est_xyz_functions.c
  Mesh/d4est_mortars.c
  Mesh/d4est_mortars_aux.c
  Mesh/d4est_ghost_data.c
  Mesh/d4est_ghost.c
  Quadrature/d4est_quadrature.c
  Quadrature/d4est_quadrature_legendre.c
  Quadrature/d4est_quadrature_lobatto.c
  IO/d4est_h5.c
  IO/d4est_checkpoint.c
  IO/d4est_vtk.c
  IO/d4est_output.c
  Solver/d4est_solver_cg.c
  Solver/d4est_solver_newton.c
  Solver/krylov_petsc.c
  Solver/newton_petsc.c
  Solver/krylov_pc_multigrid.c
  Solver/multigrid.c
  Solver/multigrid_bottom_solver_cg.c
  Solver/multigrid_bottom_solver_cheby.c
  Solver/multigrid_bottom_solver_krylov_petsc.c
  Solver/multigrid_element_data_updater.c
  Solver/multigrid_logger_residual.c
  Solver/multigrid_matrix_operator.c
  Solver/multigrid_smoother_cheby.c
  Solver/multigrid_smoother_krylov_petsc.c
  Solver/cg_eigs.c
  Solver/d4est_solver_jacobian_tester.c
  Solver/d4est_solver_test_symmetry.c
  )


##
## D8EST Library
##
add_library(d8est ${SRCS} ${SRCS_3D_ONLY})
add_library(d4est ${SRCS} ${SRCS_2D_ONLY})

add_dependencies(d8est build_bundled_libs)
add_dependencies(d4est build_bundled_libs)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")


add_custom_target(p4est_bundled_libs)
add_custom_target(build_bundled_libs)


include(zlib)
include_directories(${ZLIB_INCLUDE_DIRS})

include(p4est)
include(petsc)

subdirs(ThirdParty)

SET (
  LIBS
  ${MPI_C_LIBRARIES}
  ${P4EST_LIBRARIES}
  ${PETSC_LIBRARIES}
  m
  dl
  ${CMAKE_BINARY_DIR}/ThirdParty/zlib/install/lib/libz.a
  )

if(USE_OPENBLAS)
  SET (
    LIBS
    ${LIBS}
    ${OPENBLAS_LIBRARIES}
    pthread
    gfortran
  )
endif()

if(USE_HDF5)
  SET (
    LIBS
    ${LIBS}
    ${HDF5_LIBRARIES}
    )
endif()

if(USE_TCMALLOC)
  SET (
    LIBS
    ${LIBS}
    ${TCMALLOC_LIBRARIES}
    stdc++
  )
endif()

if(NOT USE_OPENBLAS AND NOT DOWNLOAD_BLAS)
  SET (
    LIBS
    ${LIBS}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    pthread
    gfortran
  )
endif()

message("LIBS=${LIBS}")
target_link_libraries(d4est ${LIBS})
# set_target_properties(d4est PROPERTIES COMPILE_DEFINITIONS "DIM=2")
set_property(TARGET d4est APPEND_STRING PROPERTY
  COMPILE_FLAGS " -DDIM=2")


target_link_libraries(d8est ${LIBS})
# set_target_properties(d8est PROPERTIES COMPILE_DEFINITIONS "DIM=3")
set_property(TARGET d8est APPEND_STRING PROPERTY
  COMPILE_FLAGS " -DDIM=3")


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
    ${CMAKE_SOURCE_DIR}/Problems/${parent_dir}/*.input
    ${CMAKE_CURRENT_BINARY_DIR}/${parent_dir}/)  
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
    ${CMAKE_SOURCE_DIR}/Problems/${parent_dir}/*.input
    ${CMAKE_CURRENT_BINARY_DIR}/${parent_dir}/)  
endfunction()


function(add_d4est_test_3d exec src parent_dir test_name)
  add_executable(${exec} ${CMAKE_SOURCE_DIR}/Tests/${parent_dir}/${src})
  target_link_libraries(${exec} d8est)
  add_test(${test_name} ${exec})
  target_compile_definitions(${exec} PRIVATE DIM=3)

  add_custom_command(
    TARGET ${exec} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/Tests/${parent_dir}/*.input
    ${CMAKE_CURRENT_BINARY_DIR}/)   
endfunction()


function(add_d4est_test_2d exec src parent_dir test_name)
  add_executable(${exec} ${CMAKE_SOURCE_DIR}/Tests/${parent_dir}/${src})
  target_link_libraries(${exec} d4est)
  add_test(${test_name} ${exec})
  target_compile_definitions(${exec} PRIVATE DIM=2)

  add_custom_command(
    TARGET ${exec} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/Tests/${parent_dir}/*.input
    ${CMAKE_CURRENT_BINARY_DIR}/)   
endfunction()

##
## Tests
##
enable_testing()
set(CTEST_MEMORYCHECK_COMMAND /usr/bin/valgrind)
set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full --suppressions=${PROJECT_SOURCE_DIR}/suppressions.txt" )

add_d4est_test_3d(test_d4est_vtk test_d4est_vtk.c . test_d4est_vtk)
add_d4est_test_3d(test_d4est_operators_interpolate test_d4est_operators_interpolate.c . test_d4est_operators_interpolate)
add_d4est_test_3d(test_d4est_operators test_d4est_operators.c . test_d4est_operators)
add_d4est_problem_3d(poisson_sinx_uniform_driver poisson_sinx_uniform.c PoissonSinx)
add_d4est_problem_3d(constant_density_star_driver constant_density_star_mgpc_newton_petsc.c ConstantDensityStar)
add_d4est_problem_3d(two_punctures_driver two_punctures_robin_sphere_hp_mgpc_newton_petsc_cactus.c TwoPunctures)
