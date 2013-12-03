if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Polyhedron_IO_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES iv2off.cpp off2iv.cpp off2off.cpp off2stl.cpp off2vrml.cpp off2wav.cpp off_bbox.cpp off_glue.cpp off_transform.cpp polyhedron2vrml.cpp polyhedron_copy.cpp terr_trian.cpp)
