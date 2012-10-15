# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Balls_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )
   
  include( ${CGAL_USE_FILE} )
  find_package( MPFI )          
  if( MPFI_FOUND ) 
    include( ${MPFI_USE_FILE} )
  endif( MPFI_FOUND )  
  
  include_directories(BEFORE ../../include)
  include_directories(BEFORE include)

  include_directories(BEFORE $ENV{CGAL_TRUNK}/Number_types/include/)
  include_directories(BEFORE $ENV{CGAL_TRUNK}/Algebraic_foundations/include/)
  include_directories(BEFORE $ENV{CGAL_TRUNK}/Arithmetic_kernel/include/)
  include_directories(BEFORE $ENV{CGAL_TRUNK}/Modular_arithmetic/include/)
  include_directories(BEFORE $ENV{CGAL_TRUNK}/Polynomial/include/)
  include_directories(BEFORE $ENV{CGAL_TRUNK}/Interval_support/include/)


  include_directories(BEFORE $ENV{CGAL_TRUNK}/Arrangement_on_surface_2/include/)
  include_directories(BEFORE $ENV{CGAL_TRUNK}/Envelope_3/include/)
  include_directories(BEFORE $ENV{CGAL_TRUNK}/Mesh_3/include/)
  include_directories(BEFORE $ENV{CGAL_TRUNK}/Intersections_3/include/)
  
  include_directories(BEFORE $ENV{CGAL_CAN}/Algebraic_kernel_d/include/)
  include_directories(BEFORE $ENV{CGAL_CAN}/Arrangement_on_surface_2_extension/include/)
  include_directories(BEFORE $ENV{CGAL_EXP}/Algebraic_kernel_d/include/)
  include_directories(BEFORE $ENV{CGAL_EXP}/Curved_kernel_via_analysis_2/include/)
  

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../../include)
  include_directories (BEFORE ../../test/Voronoi_diagram_of_lines_3/include)
  include_directories (BEFORE ./include)

  create_single_source_cgal_program( "Voronoi_cell_of_line_3.cpp" )
  create_single_source_cgal_program( "VDOL_generate_mesh_3.cpp" )
  create_single_source_cgal_program( "test_mesh_io.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

