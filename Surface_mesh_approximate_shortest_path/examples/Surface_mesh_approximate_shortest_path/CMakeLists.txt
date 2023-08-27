# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Surface_mesh_approximate_shortest_path_Examples)

find_package(CGAL REQUIRED)

#find_package(Eigen3 3.3.0 QUIET)
find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: These examples require the Eigen library (3.3 or greater), and will not be compiled.")
  return()
endif()

# include for local directory
include_directories(BEFORE include)

#create_single_source_cgal_program("approximate_shortest_path_example.cpp")
create_single_source_cgal_program("Surface_mesh_approximate_shortest_path_custom_enqueue_strategy.cpp")
#create_single_source_cgal_program("Surface_mesh_approximate_shortest_path_anisotropic_example.cpp")
#create_single_source_cgal_program("Surface_mesh_approximate_shortest_path_isotropic_anisotropic_parts.cpp")
#create_single_source_cgal_program("Shortest_path_algorithms_benchmark.cpp")
