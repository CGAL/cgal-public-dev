# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Voronoi_diagram_of_lines_3_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../../include)
  include_directories (BEFORE include)

  include_directories (BEFORE $ENV{CGAL_CAN}/Number_types/include)
  include_directories (BEFORE $ENV{CGAL_CAN}/Algebraic_foundations/include)          
  include_directories (BEFORE $ENV{CGAL_CAN}/Arithmetic_kernel/include)           
  include_directories (BEFORE $ENV{CGAL_CAN}/Circular_kernel_3/include)          
  include_directories (BEFORE $ENV{CGAL_CAN}/Number_types/include) 
  include_directories (BEFORE $ENV{CGAL_CAN}/Algebraic_kernel_d/include)        
  include_directories (BEFORE $ENV{CGAL_CAN}/Combinatorial_map/include)     
  include_directories (BEFORE $ENV{CGAL_CAN}/Algebraic_kernel_for_circles/include)    
  include_directories (BEFORE $ENV{CGAL_CAN}/Boolean_set_operations_2/include)     
  include_directories (BEFORE $ENV{CGAL_CAN}/Segment_Delaunay_graph_2/include) 
  include_directories (BEFORE $ENV{CGAL_CAN}/Algebraic_kernel_for_spheres/include)    
  include_directories (BEFORE $ENV{CGAL_CAN}/Circular_kernel_2/include)          
  include_directories (BEFORE $ENV{CGAL_CAN}/Kernel_d/include)                   
  include_directories (BEFORE $ENV{CGAL_CAN}/Spatial_sorting/include) 

 # include_directories(BEFORE $ENV{CGAL_EXP}/Algebraic_kernel_d/include/)

  create_single_source_cgal_program( "envelope_3.cpp" )
  create_single_source_cgal_program( "insertion_bug.cpp" )
  create_single_source_cgal_program( "test_envelope_3.cpp" )
  create_single_source_cgal_program( "Voronoi_diagram_of_lines_3.cpp" )
  create_single_source_cgal_program( "Voronoi_diagram_of_lines_hierarchy_3.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

