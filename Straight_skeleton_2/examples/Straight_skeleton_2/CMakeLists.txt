if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Straight_skeleton_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES Create_offset_polygons_2.cpp Create_saop_from_polygon_with_holes_2.cpp Create_skeleton_and_offset_polygons_2.cpp Create_skeleton_and_offset_polygons_with_holes_2.cpp Create_straight_skeleton_2.cpp Create_straight_skeleton_from_polygon_with_holes_2.cpp Low_level_API.cpp Show_offset_polygon.cpp Show_straight_skeleton.cpp)
