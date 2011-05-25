# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( Lines_through_segments_3 )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS Qt4 Core )

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

# include for local package
include_directories( BEFORE ../../include )
include_directories (/home/asaf/trees/cgal/branches/experimental-packages/AOS_2_rational_arc_traits_d_1/include)
include_directories (/home/asaf/trees/cgal/trunk/Algebraic_kernel_d/include)

# include helper file
include( ${CGAL_USE_FILE} )

# Qt4
set( QT_USE_QTXML     true )
set( QT_USE_QTMAIN    true )
set( QT_USE_QTSCRIPT  true )
set( QT_USE_QTOPENGL  true )

find_package(Qt4)  

if ( NOT QT_FOUND )

  message(STATUS "This project requires the Qt4 library, and will not be compiled.")
  return()  

endif()


# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )


if ( CGAL_Qt4_FOUND AND QT_FOUND )

  include( ${QT_USE_FILE} )
  include_directories( ${QT_INCLUDE_DIR} )  

endif()

#create_single_source_cgal_program( "lines_through_4_polytopes.cpp" )
create_single_source_cgal_program( "lines_through_4_segments.cpp" )
#create_single_source_cgal_program( "lines_through_4_segments_rat_traits_debug.cpp" )
#create_single_source_cgal_program( "lines_through_4_segments_rat_traits_debug1.cpp" )
#create_single_source_cgal_program( "lines_through_4_segments_rat_traits_debug2.cpp" )

