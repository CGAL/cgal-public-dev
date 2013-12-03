if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Boolean_set_operations_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES bezier_traits_adapter2.cpp bezier_traits_adapter.cpp circle_segment.cpp conic_traits_adapter.cpp connect_polygon.cpp do_intersect.cpp dxf_union.cpp sequence.cpp set_union.cpp simple_join_intersect.cpp symmetric_difference.cpp)
