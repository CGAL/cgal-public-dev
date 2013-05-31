# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( Symbolic )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()


# CGAL and its components
find_package( CGAL QUIET COMPONENTS Symbolic Core )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package
include_directories( BEFORE ../../include )


# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "flag.cpp" )


