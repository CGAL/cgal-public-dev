# based on a script created by cgal_create_cmake_script

project( Lifting_kernel_d_test )

cmake_minimum_required(VERSION 2.8.11)

find_package( CGAL QUIET )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories( BEFORE ../../include include )

  create_single_source_cgal_program( "lifted_orientation_d.cpp" )

else( CGAL_FOUND )

  message(STATUS
    "These tests require the CGAL library and will not be compiled.")

endif( CGAL_FOUND )
