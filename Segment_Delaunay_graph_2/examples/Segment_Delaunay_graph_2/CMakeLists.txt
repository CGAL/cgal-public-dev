if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Segment_Delaunay_graph_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES sdg-count-sites.cpp sdg-fast-sp.cpp sdg-fast-sp-polygon.cpp sdg-filtered-traits.cpp sdg-info-set.cpp sdg-red-blue-info.cpp sdg-voronoi-edges.cpp)
