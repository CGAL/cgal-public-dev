# This is the CMake script for compiling this folder.

project( Implicit_surface_reconstruction_3_Demo)

cmake_minimum_required(VERSION 2.8.11)

# Test if compiler has C++11
list(FIND CMAKE_CXX_COMPILE_FEATURES cxx_generalized_initializers has_cpp11)
if(has_cpp11 LESS 0)
  message(STATUS "NOTICE: These tests require a C++11 compiler and will not be compiled.")
  return()
endif()

# Use C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Find CGAL
find_package(CGAL QUIET)

# Complier flags
#add_definitions( -Wall )
#add_definitions( -g )
#add_definitions( -O0 )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  include_directories(BEFORE ../../include)
  
  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
  if(EIGEN3_FOUND)
    message(STATUS "Eigen = ${EIGEN3_USE_FILE}")
    include( ${EIGEN3_USE_FILE} )
  else()
    message(STATUS "NOTICE: Eigen 3.1.0 (or greater) was not found, tests won't be compiled.")
  endif()
  
  find_package(Spectra REQUIRED)
  if(SPECTRA_FOUND)
    include_directories("${SPECTRA_INCLUDE_DIR}")
  else()
    message(STATUS "NOTICE: Spectra was not found, tests won't be compiled.")
  endif()

  if(NOT EIGEN3_FOUND OR NOT SPECTRA_FOUND)
    return()
  endif()

  # VisualC++ optimization for applications dealing with large data
  if (MSVC)
    # Allow Windows 32bit applications to use up to 3GB of RAM
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

    # Prints new compilation options
    message( STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'" )
    message( STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'" )
    message( STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'" )
    message( STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'" )
  endif()

  # Temporary debugging stuff

  # Link with Boost.ProgramOptions
  find_package(Boost QUIET COMPONENTS program_options)
  if(Boost_PROGRAM_OPTIONS_FOUND)
    if( CGAL_AUTO_LINK_ENABLED )
      message( STATUS "Boost.ProgramOptions library: found" )
    else()
      message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
    endif()
    add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
    add_definitions( "-DDEBUG_TRACE")
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
    create_single_source_cgal_program( "implicit_reconstruction_new.cpp" )
  else()
    message(STATUS "NOTICE: Test implicit_reconstruction_new.cpp requires boost program_option and will not be compiled.")
  endif()

else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

