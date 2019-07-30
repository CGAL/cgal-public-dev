# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Implicit_surface_reconstruction_3_Tests )

cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_CXX_STANDARD 14)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

    include_directories (BEFORE "include")

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

  find_package(Boost COMPONENTS system filesystem REQUIRED)
  add_definitions( "-BOOST_FILESYSTEM_LIBRARY" )
  add_definitions( "-BOOST_SYSTEM_LIBRARY" )
  add_definitions( "-DDEBUG_TRACE")
  list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_FILESYSTEM_LIBRARY})
  list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_SYSTEM_LIBRARY})
  create_single_source_cgal_program( "isr_test_mean_distance_mtp.cpp" )
  create_single_source_cgal_program( "isr_test_mean_distance_ptm.cpp" )
  create_single_source_cgal_program( "isr_test_hausdorff_ptm.cpp" )
  create_single_source_cgal_program( "isr_test_hausdorff_mtp.cpp" )
  create_single_source_cgal_program( "isr_test_mean_ang_dev.cpp" )
  create_single_source_cgal_program( "isr_test_check_genus.cpp" )
  create_single_source_cgal_program( "isr_test_hole_filling.cpp" )
  create_single_source_cgal_program( "isr_test_check_cc.cpp" )
  create_single_source_cgal_program( "isr_test_octree.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()


