# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Algebraic_curve_kernel_2_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  include( CGAL_VersionUtils ) 
  
  find_package( MPFI )
  if( MPFI_FOUND ) 
   include( ${MPFI_USE_FILE} )
  endif( MPFI_FOUND )


  include_directories (BEFORE include)

  include_directories (BEFORE ../../include)

  include_directories(BEFORE $ENV{CGAL}/trunk/Arrangement_on_surface_2/include/)
  include_directories(BEFORE $ENV{CGAL}/candidate-packages/Algebraic_kernel_d/include/)
  include_directories(BEFORE $ENV{CGAL}/experimental-packages/Curved_kernel_via_analysis_2/include/)
  include_directories(BEFORE $ENV{CGAL}/candidate-packages/Arrangement_on_surface_2_extension/include/)

  include_directories(BEFORE $ENV{CGAL}/trunk/STL_Extension/include/)
  include_directories(BEFORE $ENV{CGAL}/trunk/Testsuite/include/)
  include_directories(BEFORE $ENV{CGAL}/trunk/Number_types/include/)
  include_directories(BEFORE $ENV{CGAL}/trunk/Arithmetic_kernel/include/)
  include_directories(BEFORE $ENV{CGAL}/trunk/Interval_support/include/)
  include_directories(BEFORE $ENV{CGAL}/trunk/Algebraic_foundations/include/)
  include_directories(BEFORE $ENV{CGAL}/trunk/Modular_arithmetic/include/)
  include_directories(BEFORE $ENV{CGAL}/trunk/Polynomial/include/)


  create_single_source_cgal_program( "arrangements_2.cpp" )
  create_single_source_cgal_program( "curve_analyser_2.cpp" )
  create_single_source_cgal_program( "curve_pair_analyser_2.cpp" )
  create_single_source_cgal_program( "rotated_arrangements_2.cpp" )
  create_single_source_cgal_program( "rotated_curve_analyser_2.cpp" )
  create_single_source_cgal_program( "Bezier_curves_with_ak_2.cpp" )
  create_single_source_cgal_program( "generate_bezier_curves.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

