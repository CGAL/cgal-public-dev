# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Shape_regularization_Examples )


cmake_minimum_required(VERSION 3.1)

# CGAL and its components
find_package( CGAL QUIET  )
if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
if ( NOT EIGEN3_FOUND )
  message(STATUS "This project requires the Eigen library, and will not be compiled.")
  return() 
endif()

include( ${EIGEN3_USE_FILE} )

# Boost and its components
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()


# Creating entries for all C++ files with "main" routine
# ##########################################################

create_single_source_cgal_program( "segment_regularization.cpp" )


