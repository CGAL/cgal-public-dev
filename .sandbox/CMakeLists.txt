# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required( VERSION 3.1...3.23 )
project( gsoc2023 )

set( CMAKE_CXX_EXTENSIONS OFF )

include( CheckIPOSupported )
check_ipo_supported( RESULT ipo_supported )
if( ipo_supported )
  set( CMAKE_INTERPROCEDURAL_OPTIMIZATION True )
endif()

add_subdirectory(src build)
