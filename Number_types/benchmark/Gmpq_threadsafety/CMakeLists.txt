# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Gmpq_threadsafety_ )

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

  include_directories (BEFORE "../../include")

  create_single_source_cgal_program( "dt3_bench.cpp" )
  create_single_source_cgal_program( "dt3_bench_gmpxx.cpp" )
  
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    create_single_source_cgal_program( "showing_the_bug.cpp" )
  else()
    message(STATUS "OpenMP not found, showing_the_bug.cpp will not be compiled.")
  endif()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

