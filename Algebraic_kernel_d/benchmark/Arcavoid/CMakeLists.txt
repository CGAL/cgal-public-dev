# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Arcavoid_example )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../include)

  include_directories (BEFORE include)

  create_single_source_cgal_program( "Abs_void.cpp" )
  create_single_source_cgal_program( "Absolute_void.cpp" )
  create_single_source_cgal_program( "Arcavoid.cpp" )
  create_single_source_cgal_program( "Bitstream_descartes.cpp" )
  create_single_source_cgal_program( "Descartes.cpp" )
  create_single_source_cgal_program( "RS.cpp" )
  create_single_source_cgal_program( "convert_polynomial.cpp" )
  create_single_source_cgal_program( "generators.cpp" )
  create_single_source_cgal_program( "void.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

