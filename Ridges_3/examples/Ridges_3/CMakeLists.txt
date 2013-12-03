if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Ridges_3
  REQUIRED_DEPENDENCIES CGAL Eigen3
  OPTIONAL_DEPENDENCIES Boost_program_options
  SOURCES Compute_Ridges_Umbilics.cpp PolyhedralSurf.cpp)

# Add the same targets for a LAPACK executable.
CGAL_example(SOURCE Compute_Ridges_Umbilics.cpp TARGET Compute_Ridges_Umbilics_LAPACK
  REQUIRED_DEPENDENCIES CGAL LAPACK
  OPTIONAL_DEPENDENCIES Boost_program_options)
CGAL_example(SOURCE PolyhedralSurf.cpp TARGET PolyhedralSurf_LAPACK
  REQUIRED_DEPENDENCIES CGAL LAPACK
  OPTIONAL_DEPENDENCIES Boost_program_options)
