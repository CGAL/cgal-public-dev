# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Algebraic_kernel_d_test )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

# CGAL and its components
find_package(CGAL QUIET COMPONENTS Core MPFI RS RS3 )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )

include( CGAL_CreateSingleSourceCGALProgram )

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory
include_directories( BEFORE include )

# include for local package
include_directories( BEFORE ../../include )

# TODO 2012 these includes will disappear when we merge
# Numerical_algebraic_kernel_d, Algebraic_kernel_2,
# exp/Algebraic_kernel_d and Arcavoid in here
include_directories( BEFORE $ENV{HOME}/CGAL/SVN/branches/experimental-packages/Algebraic_kernel_2/include)
include_directories( BEFORE $ENV{HOME}/CGAL/SVN/branches//unsorted-branches/eric/Numerical_algebraic_kernel_d/include)
include_directories( BEFORE $ENV{HOME}/CGAL/SVN/branches/experimental-packages/Algebraic_kernel_d/include)
include_directories( BEFORE $ENV{HOME}/CGAL/SVN/branches/unsorted-branches/akobel/Arcavoid/include)

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

### 1d ###
create_single_source_cgal_program( "Descartes.cpp" )
create_single_source_cgal_program( "Bitstream_descartes.cpp" )
create_single_source_cgal_program( "isolators.cpp" )
create_single_source_cgal_program( "Real_embeddable_traits_extension.cpp" )
create_single_source_cgal_program( "Algebraic_real_d_1.cpp" )
create_single_source_cgal_program( "Algebraic_kernel_d_1.cpp" )
if(RS_FOUND)
  create_single_source_cgal_program( "Algebraic_kernel_rs_gmpq_d_1.cpp" )
  create_single_source_cgal_program( "Algebraic_kernel_rs_gmpz_d_1.cpp" )
  create_single_source_cgal_program( "rs_isolator.cpp" )
else()
  message(STATUS "NOTICE: Some tests require the RS library, and will not be compiled.")
endif()
 
### 2d ###
create_single_source_cgal_program( "Algebraic_kernel_d_2.cpp" )
create_single_source_cgal_program( "Algebraic_kernel_d_2_v1.cpp" )
create_single_source_cgal_program( "Algebraic_kernel_d_2_wo_teissier.cpp" )
create_single_source_cgal_program( "Curve_analysis_2.cpp" )
create_single_source_cgal_program( "Curve_analysis_2_v1.cpp" )
create_single_source_cgal_program( "Curve_analysis_2_wo_teissier.cpp" )
create_single_source_cgal_program( "Curve_analysis_2_mkbs_lifter.cpp" )
create_single_source_cgal_program( "Curve_pair_analysis_2.cpp" )
create_single_source_cgal_program( "Curve_pair_analysis_2_v1.cpp" )
create_single_source_cgal_program( "Curve_pair_analysis_2_wo_teissier.cpp" )
create_single_source_cgal_program( "Algebraic_curve_kernel_2.cpp" )
create_single_source_cgal_program( "Algebraic_curve_kernel_2_v1.cpp" )
create_single_source_cgal_program( "Algebraic_curve_kernel_2_wo_teissier.cpp" )
create_single_source_cgal_program( "algebraic_curve_kernel_2_tools.cpp" )
