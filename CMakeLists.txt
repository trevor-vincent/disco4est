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
endif()

if(SCINET)
  set(USE_OPENBLAS ON)
  set(OPENBLAS_DIR "/scinet/gpc/Libraries/OpenBLAS/singlethreaded/")
  set(OPENBLAS_ROOT "/scinet/gpc/Libraries/OpenBLAS/singlethreaded/")
  set(OPENBLAS_LIBRARIES
    "/scinet/gpc/Libraries/OpenBLAS/singlethreaded/lib/libopenblas.a"
    )
  set(OPENBLAS_INCLUDE_DIR "/scinet/gpc/Libraries/OpenBLAS/singlethreaded/include")  
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
  )

SET (
  SRCS_2D_ONLY
  Geometry/d4est_geometry_disk.c
  )

SET (
  SRCS
  Utilities/ini.c
  Utilities/arbquad.c
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
  # Quadrature/d4est_quadrature_compactified.c
  Quadrature/d4est_quadrature.c
  Quadrature/d4est_quadrature_legendre.c
  Quadrature/d4est_quadrature_lobatto.c
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
# set( MEMORYCHECK_SUPPRESSIONS_FILE "/usr/share/openmpi/openmpi-valgrind.supp" )
# set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full --suppressions=/usr/share/openmpi/openmpi-valgrind.supp" )
set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full --suppressions=${PROJECT_SOURCE_DIR}/suppressions.txt" )

# add_executable(test_d4est_linalg Tests/test_d4est_linalg.c)
# target_link_libraries(test_d4est_linalg d4est)
# add_test(test_d4est_linalg_test test_d4est_linalg)
# target_compile_definitions(test_d4est_linalg PRIVATE NDEBUG)


# add_executable(test_d4est_operators Tests/test_d4est_operators.c)
# target_link_libraries(test_d4est_operators d8est)
# add_test(test_d4est_operators_test test_d4est_operators)
# target_compile_definitions(test_d4est_operators PRIVATE NDEBUG DIM=3)

# add_executable(test_d4est_kron Tests/test_d4est_kron.c)
# target_link_libraries(test_d4est_kron d8est)
# add_test(test_d4est_kron_test test_d4est_kron)
# target_compile_definitions(test_d4est_kron PRIVATE NDEBUG DIM=3)

# add_executable(test_d4est_amr Tests/test_d4est_amr.c)
# target_link_libraries(test_d4est_amr d8est)
# add_test(test_d4est_amr_test test_d4est_amr)
# target_compile_definitions(test_d4est_amr PRIVATE DIM=3)

# add_executable(test_d4est_poisson_1_brick_3d Tests/test_d4est_poisson_1_brick.c)
# target_link_libraries(test_d4est_poisson_1_brick_3d d8est)
# add_test(test_d4est_poisson_1_brick_3d_test test_d4est_poisson_1_brick_3d)
# target_compile_definitions(test_d4est_poisson_1_brick_3d PRIVATE DIM=3)

# add_executable(test_d4est_poisson_2_brick_2d Tests/test_d4est_poisson_2_brick.c)
# target_link_libraries(test_d4est_poisson_2_brick_2d d4est)
# add_test(test_d4est_poisson_2_brick_2d_test test_d4est_poisson_2_brick_2d)
# target_compile_definitions(test_d4est_poisson_2_brick_2d PRIVATE DIM=2)

# add_executable(test_d4est_poisson_1_brick_2d Tests/test_d4est_poisson_1_brick.c)
# target_link_libraries(test_d4est_poisson_1_brick_2d d4est)
# add_test(test_d4est_poisson_1_brick_2d_test test_d4est_poisson_1_brick_2d)
# target_compile_definitions(test_d4est_poisson_1_brick_2d PRIVATE DIM=2)

# add_executable(test_d4est_geometry Tests/test_d4est_geometry.c)
# target_link_libraries(test_d4est_geometry d8est)
# add_test(test_d4est_geometry_test test_d4est_geometry)
# target_compile_definitions(test_d4est_geometry PRIVATE NDEBUG DIM=3)

# add_executable(test_d4est_geometry_2d Tests/test_d4est_geometry.c)
# target_link_libraries(test_d4est_geometry_2d d4est)
# add_test(test_d4est_geometry_test_2d test_d4est_geometry_2d)
# target_compile_definitions(test_d4est_geometry_2d PRIVATE NDEBUG DIM=2)

# add_executable(test_d4est_mass_projection Tests/test_d4est_mass_projection.c)
# target_link_libraries(test_d4est_mass_projection d8est)
# add_test(test_d4est_mass_projection_test test_d4est_mass_projection)
# target_compile_definitions(test_d4est_mass_projection PRIVATE DIM=3)

# add_d4est_test_3d(test_d4est_boyen_york_model_fcns test_d4est_boyen_york_model_fcns.c . test_d4est_boyen_york_model_fcns)
# add_d4est_test_3d(test_d4est_okendon_fcns test_d4est_okendon_fcns.c . test_d4est_okendon_fcns)
# add_d4est_test_2d(test_d4est_poisson_3_brick test_d4est_poisson_3_brick.c . test_d4est_poisson_3_brick)
# add_d4est_test_2d(test_d4est_poisson_4_brick test_d4est_poisson_4_brick.c . test_d4est_poisson_4_brick)
# add_d4est_test_3d(test_d4est_poisson_4_cubed_sphere test_d4est_poisson_4_cubed_sphere.c . test_d4est_poisson_4_cubed_sphere)
# add_d4est_test_3d(test_d4est_poisson_consistency_curved_3d test_d4est_poisson_consistency_curved.c . test_d4est_poisson_consistency_curved_3d)
# add_d4est_test_2d(test_d4est_poisson_consistency_curved_2d test_d4est_poisson_consistency_curved.c . test_d4est_poisson_consistency_curved_2d)
# add_d4est_test_3d(test_d4est_puncture_fcns test_d4est_puncture_fcns.c . test_d4est_puncture_fcns)
# add_d4est_test_3d(test_d4est_constant_density_star_fcns test_d4est_constant_density_star_fcns.c . test_d4est_constant_density_star_fcns)
# add_d4est_test_3d(test_d4est_poisson_symmetry test_d4est_poisson_symmetry.c . test_d4est_poisson_symmetry)
# add_d4est_test_3d(test_d4est_stiffness test_d4est_stiffness.c . test_d4est_stiffness)
add_d4est_test_3d(test_d4est_brick_fcns test_d4est_brick_fcns.c . test_d4est_brick_fcns)

# add_executable(test_d4est_quadrature_compactified Tests/test_d4est_quadrature_compactified.c)
# target_link_libraries(test_d4est_quadrature_compactified d8est)
# add_test(test_d4est_quadrature_compactified_test test_d4est_quadrature_compactified)
# target_compile_definitions(test_d4est_quadrature_compactified PRIVATE DIM=3)


# add_executable(test_d4est_quadrature_compactified_2 Tests/test_d4est_quadrature_compactified_2.c)
# target_link_libraries(test_d4est_quadrature_compactified_2 d8est)
# add_test(test_d4est_quadrature_compactified_2_test test_d4est_quadrature_compactified_2)
# target_compile_definitions(test_d4est_quadrature_compactified_2 PRIVATE DIM=3)



# add_executable(test_d4est_mortar_geometric_terms Tests/test_d4est_mortar_geometric_terms.c)
# target_link_libraries(test_d4est_mortar_geometric_terms d8est)
# add_test(test_d4est_mortar_geometric_terms_test test_d4est_mortar_geometric_terms)
# target_compile_definitions(test_d4est_mortar_geometric_terms PRIVATE DIM=3)

# add_executable(test_d4est_dictionary Tests/test_d4est_dictionary.c)
# target_link_libraries(test_d4est_dictionary d4est)
# add_test(test_d4est_dictionary_test test_d4est_dictionary)
# target_compile_definitions(test_d4est_dictionary PRIVATE D4EST_DIM=3 D4EST_SERIAL D4EST_DEBUG)


# get_directory_property( DirDefs COMPILE_DEFINITIONS )
# message( "DUCK COMPILE_DEFINITIONS = ${DirDefs}" )
# add_d4est_problem_3d(testcompile1 testcompile.c TestCompile)
# add_d4est_problem_2d(testcompile2 testcompile.c TestCompile)
# add_d4est_problem_3d(constantdensitystar_driver constantdensitystar.c ConstantDensityStar)
add_d4est_problem_3d(poisson_sinx_uniform_driver poisson_sinx_uniform.c PoissonSinx)
add_d4est_problem_3d(poisson_rsquared_driver poisson_rsquared.c PoissonRSquared)


add_test(test_d4est_brick_fcns test_d4est_brick_fcns)

add_d4est_problem_2d(stamm_krylov_driver stamm_krylov.c Stamm)
add_d4est_problem_2d(stamm_multigrid_driver stamm_multigrid.c Stamm)
add_d4est_problem_2d(stamm_multigrid_pc_driver stamm_multigrid_pc.c Stamm)
add_d4est_problem_3d(stamm_multigrid_pc_sphere_driver stamm_multigrid_pc_sphere.c Stamm)
add_d4est_problem_2d(poisson_exp_robin_driver poisson_exp_robin.c PoissonExpRobin)
add_d4est_problem_3d(constant_density_star_driver constant_density_star.c ConstantDensityStar)
add_d4est_problem_3d(constant_density_star_mgpc_driver constant_density_star_mgpc.c ConstantDensityStar)
# add_d4est_problem_3d(constant_density_star_inbin_driver constant_density_star_inbin.c ConstantDensityStar)
# add_d4est_problem_3d(two_punctures_robin_brick_driver two_punctures_robin_brick.c TwoPunctures)
# add_d4est_problem_3d(two_punctures_robin_sphere_driver two_punctures_robin_sphere.c TwoPunctures)
add_d4est_problem_3d(two_punctures_robin_sphere_hp_driver two_punctures_robin_sphere_hp.c TwoPunctures)
add_d4est_problem_3d(two_punctures_robin_sphere_hp_mgpc_driver two_punctures_robin_sphere_hp_mgpc.c TwoPunctures)
add_d4est_problem_3d(two_punctures_robin_sphere_hp_mgpc_noinbin_driver two_punctures_robin_sphere_hp_mgpc.c TwoPunctures)
add_d4est_problem_3d(two_punctures_robin_brick_hp_mgpc_driver two_punctures_robin_brick_hp_mgpc_noinbin.c TwoPunctures)
add_d4est_problem_3d(two_punctures_robin_sphere_hp_mgpc_noinbin_tree12_driver two_punctures_robin_sphere_hp_mgpc_noinbin_tree12.c TwoPunctures)
add_d4est_problem_3d(two_punctures_robin_sphere_hp_mgpc_tree12_driver two_punctures_robin_sphere_hp_mgpc_tree12.c TwoPunctures)
add_d4est_problem_3d(random_punctures_robin_sphere_hp_driver random_punctures_robin_sphere_hp.c TwoPunctures)
# add_d4est_problem_3d(two_punctures_driver two_punctures.c TwoPunctures)
# add_d4est_problem_2d(constant_density_star_driver stamm_krylov.c Stamm)
# add_d4est_problem_3d(okendon_uniform_driver okendon_uniform.c Okendon)
# add_d4est_problem_3d(boyen_york_model_driver boyen_york_model_uniform.c BoyenYorkModel)
# add_d4est_problem_2d(stammhpmultigrid_driver stammhpmultigrid.c StammHPMultigrid)
# add_d4est_problem_2d(test_uniform_wedge_driver test_uniform_wedge.c UniformExamples/Wedge)
# add_d4est_problem_2d(test_uniform_wedge_multigrid_driver test_uniform_wedge_multigrid.c UniformExamples/WedgeMultigrid)
# add_d4est_problem_3d(test_uniform_cubedsphere_poisson_driver test_uniform_cubedsphere_poisson.c UniformExamples/CubedSpherePoisson)
# add_d4est_problem_3d(test_uniform_cubedsphere_helmholtz_driver test_uniform_cubedsphere_helmholtz.c UniformExamples/CubedSphereHelmholtz)
# add_d4est_problem_3d(test_uniform_cubedsphere_helmholtz_with_hole_driver test_uniform_cubedsphere_with_hole_helmholtz.c UniformExamples/CubedSphereWithHoleHelmholtz)
# add_d4est_problem_3d(test_uniform_cubedsphere_helmholtz_with_hole_multigrid_driver test_uniform_cubed_sphere_with_hole_helmholtz_multigrid.c UniformExamples/CubedSphereWithHoleHelmholtzMultigrid)
# add_d4est_problem_3d(test_uniform_kidder_driver test_uniform_kidder.c UniformExamples/Kidder)
# add_d4est_problem_3d(test_jacobian_driver test_jacobian.c TestJacobian)
# add_d4est_problem_2d(testmultigridmatrixop_driver testmultigridmatrixop.c TestMultigridMatrixOp)
# add_d4est_problem_3d(testcurvedipflux_driver testcurvedipflux.c TestCurvedIPFlux)
# add_d4est_problem_2d(testcurvedipflux_driver_2d testcurvedipflux.c TestCurvedIPFlux)
# add_d4est_problem_2d(testcurvedprimal_driver_2d testcurvedprimal.c TestCurvedPrimal)
# add_d4est_problem_3d(testcurvedprimal_driver_3d testcurvedprimal.c TestCurvedPrimal)
# add_d4est_problem_3d(testcubeprimal_driver_3d testcubeprimal.c TestCubePrimal)
# add_d4est_problem_3d(testhelmholtzmultigrid_driver_3d testhelmholtzmultigrid.c TestHelmholtzMultigrid)
# add_d4est_problem_3d(testhelmholtzmultigridsphere_driver_3d testhelmholtzmultigridsphere.c TestHelmholtzMultigrid)
# add_d4est_problem_3d(testhelmholtzmultigridsphere_withhole_driver_3d testhelmholtzmultigridsphere_withhole.c TestHelmholtzMultigrid)
# add_d4est_problem_3d(testhelmholtzmultigridsphere_withhole_tp_driver_3d testhelmholtzmultigridsphere_withhole_tp.c TestHelmholtzMultigrid)
#add_d4est_problem_3d(constantdensitystarprimal_driver constantdensitystarprimal.c ConstantDensityStarPrimal)
# add_d4est_problem_3d(cdssphere_driver cdssphere.c ConstantDensityStarSphere)
# add_d4est_problem_3d(testcurvedipflux_driver_3d testcurvedipflux.c TestCurvedIPFlux)
# add_d4est_problem_3d(testcurvedshell_driver_3d testcurvedshell.c TestCurvedShell)
# add_d4est_problem_3d(twopuncturessphere_driver twopuncturessphere.c TwoPuncturesSphere)
# add_d4est_problem_3d(twopuncturessphere_mg_driver twopuncturessphere_mg.c TwoPuncturesSphere)
# add_d4est_problem_3d(twopunctures_testfcns_driver twopunctures_testfcns.c TwoPuncturesSphere)
# add_d4est_problem_3d(twopunctures_estonly_driver twopuncturessphere_estonly.c TwoPuncturesSphere)
# add_d4est_problem_3d(twopunctures_nomg_driver twopunctures_nomg.c TwoPuncturesSphere)
# add_d4est_problem_3d(twopunctures_test_driver twopuncturessphere_test.c TwoPuncturesSphere)
# add_d4est_problem_2d(twopunctures_estonly_driver_2d twopuncturessphere_estonly.c TwoPuncturesSphere)
# add_d4est_problem_3d(testmortarjacobianterms_driver_3d testmortarjacobianterms.c TestMortarJacobianTerms)
# add_d4est_problem_2d(testmortarjacobianterms_driver_2d testmortarjacobianterms.c TestMortarJacobianTerms)
# add_d4est_problem_3d(testcurvedipfluxoncube_driver testcurvedipfluxoncube.c TestCurvedIPFluxOnCube)
# add_d4est_problem_2d(testcurvedipfluxoncube2d_driver testcurvedipfluxoncube.c TestCurvedIPFluxOnCube)
# add_d4est_problem_3d(testnoncurvedipfluxoncube_driver testnoncurvedipfluxoncube.c TestNonCurvedIPFluxOnCube)
# add_d4est_problem_3d(twopuncturescube_driver twopuncturescube.c TwoPuncturesCube)
# add_d4est_curved_problem_3d(testcompactification_driver testcompactification.c TestCompactification)
# add_d4est_problem_3d(testcompactification_driver testcompactification.c TestCompactification)

# add_executable(testcompile1 driver.c ./Problems/TestCompile/testcompile.c)
# target_link_libraries(testcompile1 d8est)
# set_target_properties(testcompile1 PROPERTIES COMPILE_DEFINITIONS "DIM=3")

# add_executable(testcompile2 driver.c ./Problems/TestCompile/testcompile.c)
# target_link_libraries(testcompile2 d4est)
# set_target_properties(testcompile2 PROPERTIES COMPILE_DEFINITIONS "DIM=2")
# include("${CMAKE_BINARY_DIR}/ThirdParty/petsc/install/lib/petsc/conf/PETScBuildInternal.cmake"
