if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Polyhedron_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES polyhedron_prog_color.cpp polyhedron_prog_cube.cpp polyhedron_prog_cut_cube.cpp polyhedron_prog_incr_builder.cpp polyhedron_prog_normals.cpp polyhedron_prog_off.cpp polyhedron_prog_planes.cpp polyhedron_prog_simple.cpp polyhedron_prog_subdiv.cpp polyhedron_prog_subdiv_with_boundary.cpp polyhedron_prog_tetra.cpp polyhedron_prog_traits.cpp polyhedron_prog_vector.cpp polyhedron_prog_vertex_color.cpp polyhedron_prog_vertex_normal.cpp polyhedron_self_intersection.cpp)
