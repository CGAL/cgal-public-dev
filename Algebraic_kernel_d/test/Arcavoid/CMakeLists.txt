# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( test_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  find_package( MPFI )
  if( MPFI_FOUND )
    include( ${MPFI_USE_FILE} )
  endif( MPFI_FOUND )

  include_directories (BEFORE ../include)

  include_directories (BEFORE include)

  #add_definitions( -DCGAL_GMPFR_NO_REFCOUNT=1 )

  create_single_source_cgal_program( "double.cpp" )
  create_single_source_cgal_program( "fesignal.cpp" )
  create_single_source_cgal_program( "test_Arcavoid_kernel.cpp" )
  create_single_source_cgal_program( "test_Cartesian_complex.cpp" )
  create_single_source_cgal_program( "test_Double_with_exponent.cpp" )
  create_single_source_cgal_program( "test_Gmpfi.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

