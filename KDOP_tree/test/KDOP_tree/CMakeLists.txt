# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( KDOP_tree )


cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory
include_directories( BEFORE include )

# include for local package


# Creating entries for all C++ files with "main" routine
# ##########################################################


create_single_source_cgal_program( "sixDOP.cpp" )
create_single_source_cgal_program( "kdop_test_kdop_tree.cpp" )
create_single_source_cgal_program( "kdop_test_do_intersect.cpp" )
create_single_source_cgal_program( "kdop_test_ray_query.cpp" )
create_single_source_cgal_program( "kdop_test_correctness.cpp" )
create_single_source_cgal_program( "kdop_test_plate.cpp" )
create_single_source_cgal_program( "kdop_test_distance_query.cpp" )
