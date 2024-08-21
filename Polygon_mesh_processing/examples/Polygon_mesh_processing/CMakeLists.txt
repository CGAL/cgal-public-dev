# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)
project(Polygon_mesh_processing_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("sample_example.cpp" )




find_package(Eigen3 3.2.0 QUIET) #(requires 3.2.0 or greater)
include(CGAL_Eigen3_support)

create_single_source_cgal_program( "poisson_sampling_sm_example.cpp" )
target_link_libraries(poisson_sampling_sm_example PUBLIC CGAL::Eigen3_support)