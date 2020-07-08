# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project( Polygon_mesh_processing_Examples )

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()

# include for local package
find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)

# Ajouts pour In-meshing
set(source_dir "/home/felix/code/source")
include_directories(${source_dir})

# Options
# =====================================================================================================================
option(WITH_POISSON_REF   "Compile the reference Screened Poisson Reconstruction" OFF)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "cmake")

# Dependencies & utility
# =====================================================================================================================
if(WITH_POISSON_REF)
  find_package(OpenMP REQUIRED)
endif()

function(target_compile_to_bin target)
  set_target_properties(${target} PROPERTIES
          VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin
          RUNTIME_OUTPUT_DIRECTORY      $<1:${CMAKE_CURRENT_LIST_DIR}/bin>
          ARCHIVE_OUTPUT_DIRECTORY      $<1:${CMAKE_CURRENT_LIST_DIR}/bin>
          LIBRARY_OUTPUT_DIRECTORY      $<1:${CMAKE_CURRENT_LIST_DIR}/bin>)
endfunction()
function(adjust_exe_target target)
  source_group(${source_dir} FILES ${ARGN})
  target_compile_to_bin(${target})
  target_link_libraries(${target} PUBLIC core)
endfunction()

# Poisson (reference)
# =====================================================================================================================
if(WITH_POISSON_REF)
  file(GLOB REF_SOURCES ${source_dir}/poisson_recon/*.h ${source_dir}/poisson_recon/*.inl)
  source_group("poisson_recon" FILES ${REF_SOURCES} ${source_dir}/poisson_recon/PoissonRecon.cpp ${source_dir}/poisson_recon/PlyFile.cpp)
  add_executable(poisson_recon_ref ${REF_SOURCES} ${source_dir}/poisson_recon/PoissonRecon.cpp ${source_dir}/poisson_recon/PlyFile.cpp)
  target_link_libraries(poisson_recon_ref PUBLIC OpenMP::OpenMP_CXX)
  if(MSVC)
    target_compile_options(poisson_recon_ref PUBLIC /arch:AVX2)
  endif()
  target_compile_to_bin(poisson_recon_ref)
endif()

# Poisson: custom implementation
# =====================================================================================================================
file(GLOB_RECURSE SOURCES ${source_dir}/poisson_surface/*.h ${source_dir}/poisson_surface/*.cpp ${source_dir}/poisson_surface/*.inl)
source_group(TREE ${source_dir} FILES ${SOURCES})
add_library(poisson ${SOURCES})
CGAL_target_use_Eigen(poisson)
target_include_directories(poisson PUBLIC ${source_dir})
target_compile_features(poisson PUBLIC cxx_std_17)
set_target_properties(poisson PROPERTIES FOLDER "core")
if(MSVC)
  target_compile_definitions(poisson PUBLIC _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS)
endif()

# In-meshing: core
# =====================================================================================================================
file(GLOB_RECURSE SOURCES ${source_dir}/core/*.h ${source_dir}/core/*.cpp ${source_dir}/core/*.inl)
set(SOURCES ${SOURCES} ${source_dir}/poisson_recon/Bridge.cpp ${source_dir}/poisson_recon/Bridge.h ${source_dir}/base.h ${source_dir}/debugbreak.h)
source_group(TREE ${source_dir} FILES ${SOURCES})
message("SOURCES : ")
message(${SOURCES})
add_library(core ${SOURCES})
target_link_libraries(core PUBLIC poisson)
set_target_properties(core PROPERTIES FOLDER "core")

# In-meshing: command line
# =====================================================================================================================
create_single_source_cgal_program( "in_meshing.cpp" )
CGAL_target_use_Eigen(in_meshing)
#target_link_libraries(in_meshing PUBLIC core)


adjust_exe_target(in_meshing in_meshing.cpp)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main)

# Ajouts pour In-meshing - fin

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()


find_package( OpenMesh QUIET )

if ( OpenMesh_FOUND )
include( UseOpenMesh )
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()

# include for local directory


# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program( "hausdorff_distance_remeshing_example.cpp")

if (EIGEN3_FOUND)
  create_single_source_cgal_program( "hole_filling_example.cpp" )
  CGAL_target_use_Eigen(hole_filling_example)
  create_single_source_cgal_program( "hole_filling_example_SM.cpp" )
  CGAL_target_use_Eigen(hole_filling_example_SM)
  create_single_source_cgal_program( "refine_fair_example.cpp")
  CGAL_target_use_Eigen(refine_fair_example)
  create_single_source_cgal_program( "shape_smoothing_example.cpp")
  CGAL_target_use_Eigen(shape_smoothing_example)
  create_single_source_cgal_program( "hole_filling_example_LCC.cpp" )
  CGAL_target_use_Eigen(hole_filling_example_LCC)
endif(EIGEN3_FOUND)

create_single_source_cgal_program( "self_intersections_example.cpp" )
create_single_source_cgal_program( "stitch_borders_example.cpp" )
create_single_source_cgal_program( "compute_normals_example_Polyhedron.cpp"  CXX_FEATURES cxx_range_for )
create_single_source_cgal_program( "compute_normals_example.cpp" CXX_FEATURES cxx_range_for cxx_auto_type )
create_single_source_cgal_program( "point_inside_example.cpp")
create_single_source_cgal_program( "triangulate_faces_example.cpp")
create_single_source_cgal_program( "connected_components_example.cpp")
create_single_source_cgal_program( "face_filtered_graph_example.cpp")
create_single_source_cgal_program( "orient_polygon_soup_example.cpp")
create_single_source_cgal_program( "triangulate_polyline_example.cpp")
create_single_source_cgal_program( "mesh_slicer_example.cpp")
#create_single_source_cgal_program( "remove_degeneracies_example.cpp")
create_single_source_cgal_program( "isotropic_remeshing_example.cpp")
create_single_source_cgal_program( "isotropic_remeshing_of_patch_example.cpp")
create_single_source_cgal_program( "surface_mesh_intersection.cpp")
create_single_source_cgal_program( "corefinement_SM.cpp")
create_single_source_cgal_program( "corefinement_consecutive_bool_op.cpp" )
create_single_source_cgal_program( "corefinement_difference_remeshed.cpp" )
create_single_source_cgal_program( "corefinement_mesh_union.cpp" )
create_single_source_cgal_program( "corefinement_mesh_union_and_intersection.cpp" )
create_single_source_cgal_program( "corefinement_mesh_union_with_attributes.cpp" )
create_single_source_cgal_program( "corefinement_polyhedron_union.cpp" )
create_single_source_cgal_program( "random_perturbation_SM_example.cpp" )
create_single_source_cgal_program( "corefinement_LCC.cpp")
create_single_source_cgal_program( "detect_features_example.cpp" )
create_single_source_cgal_program( "volume_connected_components.cpp" )
create_single_source_cgal_program( "manifoldness_repair_example.cpp" )
create_single_source_cgal_program( "repair_polygon_soup_example.cpp" )
create_single_source_cgal_program( "mesh_smoothing_example.cpp")
create_single_source_cgal_program( "locate_example.cpp")
create_single_source_cgal_program( "orientation_pipeline_example.cpp")
#create_single_source_cgal_program( "self_snapping_example.cpp")
#create_single_source_cgal_program( "snapping_example.cpp")


if(OpenMesh_FOUND)

create_single_source_cgal_program( "compute_normals_example_OM.cpp" )
target_link_libraries( compute_normals_example_OM PRIVATE ${OPENMESH_LIBRARIES} )

if (EIGEN3_FOUND)
  create_single_source_cgal_program( "hole_filling_example_OM.cpp" )
  target_link_libraries( hole_filling_example_OM PRIVATE ${OPENMESH_LIBRARIES} )
  CGAL_target_use_Eigen( hole_filling_example_OM )
endif()

create_single_source_cgal_program( "point_inside_example_OM.cpp")
target_link_libraries( point_inside_example_OM PRIVATE ${OPENMESH_LIBRARIES} )

create_single_source_cgal_program( "stitch_borders_example_OM.cpp" )
target_link_libraries( stitch_borders_example_OM PRIVATE ${OPENMESH_LIBRARIES} )

#create_single_source_cgal_program( "remove_degeneracies_example.cpp")
#target_link_libraries( remove_degeneracies_example PRIVATE ${OPENMESH_LIBRARIES} )

create_single_source_cgal_program( "triangulate_faces_example_OM.cpp")
target_link_libraries( triangulate_faces_example_OM PRIVATE ${OPENMESH_LIBRARIES} )
endif(OpenMesh_FOUND)

find_package( TBB )
if( TBB_FOUND )
  CGAL_target_use_TBB(self_intersections_example)
  CGAL_target_use_TBB(hausdorff_distance_remeshing_example)
else()
  message( STATUS "NOTICE: Intel TBB was not found. Sequential code will be used." )
endif()

find_package(Ceres QUIET)
if(TARGET ceres)
  target_compile_definitions( mesh_smoothing_example PRIVATE CGAL_PMP_USE_CERES_SOLVER )
  target_link_libraries( mesh_smoothing_example PRIVATE ceres )
endif(TARGET ceres)
