# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Alpha_wrap_3_Examples)

find_package(CGAL REQUIRED)

# create a target per cppfile
create_single_source_cgal_program("triangle_soup_wrap.cpp")
create_single_source_cgal_program("point_set_wrap.cpp")
create_single_source_cgal_program("wrap_from_cavity.cpp")
create_single_source_cgal_program("mixed_inputs_wrap.cpp")
create_single_source_cgal_program("volumetric_wrap.cpp")
create_single_source_cgal_program("successive_wraps.cpp")
create_single_source_cgal_program("pause_and_resume_wrapping.cpp")

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("triangle_mesh_wrap.cpp")
  target_link_libraries(triangle_mesh_wrap PRIVATE CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: Some examples use Eigen, and will not be compiled.")
endif()
