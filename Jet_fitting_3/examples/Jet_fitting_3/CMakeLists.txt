if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(NAME Jet_fitting_3_examples)

CGAL_example(SOURCE Mesh_estimation.cpp
  REQUIRED_DEPENDENCIES CGAL Eigen3
  OPTIONAL_DEPENDENCIES Boost_program_options)
CGAL_example(SOURCE Single_estimation.cpp
  REQUIRED_DEPENDENCIES CGAL Eigen3
  OPTIONAL_DEPENDENCIES Boost_program_options)

CGAL_example(TARGET Mesh_estimation_LAPACK
  SOURCE Mesh_estimation.cpp
  REQUIRED_DEPENDENCIES CGAL LAPACK
  OPTIONAL_DEPENDENCIES Boost_program_options)
CGAL_example(TARGET Single_estimation_LAPACK
  SOURCE Single_estimation.cpp
  REQUIRED_DEPENDENCIES CGAL LAPACK
  OPTIONAL_DEPENDENCIES Boost_program_options)
