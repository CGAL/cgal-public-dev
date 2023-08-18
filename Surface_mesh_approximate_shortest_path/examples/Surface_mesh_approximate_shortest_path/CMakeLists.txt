# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Surface_mesh_approximate_shortest_path_Examples)

#find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(CGAL REQUIRED)

create_single_source_cgal_program("approximate_shortest_path_example.cpp")
create_single_source_cgal_program("Surface_mesh_approximate_shortest_path_multiple_sources.cpp")
create_single_source_cgal_program("Surface_mesh_approximate_shortest_path_anisotropic_example.cpp")
create_single_source_cgal_program("Surface_mesh_approximate_shortest_path_isotropic_anisotropic_parts.cpp")
create_single_source_cgal_program("Shortest_path_algorithms_benchmark.cpp")
