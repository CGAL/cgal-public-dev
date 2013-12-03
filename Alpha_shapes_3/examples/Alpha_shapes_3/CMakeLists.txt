if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Alpha_shapes_3_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES ex_alpha_shapes_3.cpp ex_alpha_shapes_exact_alpha.cpp ex_alpha_shapes_with_fast_location_3.cpp ex_fixed_weighted_alpha_shapes_3.cpp ex_periodic_alpha_shapes_3.cpp ex_weighted_alpha_shapes_3.cpp)
