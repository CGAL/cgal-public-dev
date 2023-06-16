# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.25)
project( sandbox CXX )

find_package(CGAL REQUIRED)

# create a target per cppfile
file(
  GLOB cppfiles ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program("${cppfile}")
endforeach()
