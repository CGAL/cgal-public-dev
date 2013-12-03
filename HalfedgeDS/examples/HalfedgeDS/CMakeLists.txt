if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME HalfedgeDS_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES hds_prog_color.cpp hds_prog_compact2.cpp hds_prog_compact.cpp hds_prog_default.cpp hds_prog_edge_iterator.cpp hds_prog_graph2.cpp hds_prog_graph.cpp hds_prog_halfedge_iterator.cpp hds_prog_vector.cpp)
