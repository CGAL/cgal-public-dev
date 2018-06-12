# This is the CMake script for compiling this folder.

project( Spectral_surface_reconstruction_3_Tests )

cmake_minimum_required(VERSION 2.8.11)


# Find CGAL
find_package(CGAL QUIET)

if ( CGAL_FOUND )

  find_package(Spectra REQUIRED)
  include_directories("${SPECTRA_INCLUDE_DIR}")
  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  include_directories(BEFORE ../../include)
  #include_directories(BEFORE ../../include/CGAL/Spectra/include)

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

  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
  if(EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
    # Executables that require Eigen 3.1
    

    # Link with Boost.ProgramOptions (optional)
    find_package(Boost QUIET COMPONENTS program_options)
    if(Boost_PROGRAM_OPTIONS_FOUND)
      if( CGAL_AUTO_LINK_ENABLED )
        message( STATUS "Boost.ProgramOptions library: found" )
      else()
        message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
      endif()
      add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
      list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
      create_single_source_cgal_program( "spectral_reconstruction_test.cpp" )
    else()
      message(STATUS "NOTICE: Example spectral_reconstruction_test.cpp requires boost program_option and will not be compiled.")
    endif()
  else()

    message(STATUS "NOTICE: Some of the executables in this directory need Eigen 3.1 (or greater) and will not be compiled.")

  endif()

else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

