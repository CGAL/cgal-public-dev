# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.12...3.31)
project(Splat_surface_reconstruction_3_Examples)

# Find CGAL
find_package(CGAL REQUIRED)


# Find Eigen3
find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  # Executables that require Eigen 3
  create_single_source_cgal_program("splat_reconstruction_function.cpp")
  target_link_libraries(splat_reconstruction_function PRIVATE CGAL::Eigen3_support)
else()
  message("NOTICE: The examples require Eigen 3.1 (or greater) will not be compiled.")
endif()

# Use TBB for parallelism
if(CGAL_ACTIVATE_CONCURRENT_MESH_3)
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)
  if(TARGET CGAL::TBB_support)
    target_link_libraries(splat_reconstruction_function PRIVATE CGAL::TBB_support)
  endif()#Find TBB
endif()
