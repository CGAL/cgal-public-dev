if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Combinatorial_map_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES map_3_dynamic_onmerge.cpp map_3_foreach.cpp map_3_marks.cpp map_3_operations.cpp map_3_simple_example.cpp map_3_with_colored_facets.cpp map_4_simple_example.cpp)
