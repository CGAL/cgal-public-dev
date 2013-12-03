if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME AABB_tree_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES AABB_custom_example.cpp AABB_custom_indexed_triangle_set_array_example.cpp AABB_custom_indexed_triangle_set_example.cpp AABB_custom_triangle_soup_example.cpp AABB_face_graph_triangle_example.cpp AABB_halfedge_graph_edge_example.cpp AABB_insertion_example.cpp AABB_polyhedron_edge_example.cpp AABB_polyhedron_facet_distance_example.cpp AABB_polyhedron_facet_intersection_example.cpp AABB_segment_3_example.cpp AABB_triangle_3_example.cpp)
