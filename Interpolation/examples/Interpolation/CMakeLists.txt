if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Interpolation_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES interpolation_2.cpp linear_interpolation_2.cpp nn_coordinates_2.cpp nn_coordinates_3.cpp rn_coordinates_2.cpp sibson_interpolation_2.cpp surface_neighbor_coordinates_3.cpp)
