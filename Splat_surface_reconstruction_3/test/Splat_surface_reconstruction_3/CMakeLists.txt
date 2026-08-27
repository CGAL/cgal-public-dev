# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project( Splat_surface_reconstruction_3 )

find_package(CGAL REQUIRED)

create_single_source_cgal_program( "test.cpp" )

