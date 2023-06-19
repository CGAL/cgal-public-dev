cmake_minimum_required(VERSION 3.1...3.23)
project(Surface_mesh_approximate_shortest_path_Tests)

#include_directories(BEFORE "include")

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core)

create_single_source_cgal_program("Surface_mesh_approximate_shortest_path_traits_test.cpp")
