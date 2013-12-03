if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Triangulation_3_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES adding_handles_3.cpp color.cpp copy_triangulation_3.cpp fast_location_3.cpp find_conflicts_3.cpp info_insert_with_pair_iterator.cpp info_insert_with_pair_iterator_regular.cpp info_insert_with_transform_iterator.cpp info_insert_with_zip_iterator.cpp linking_2d_and_3d.cpp regular_3.cpp simple_triangulation_3.cpp simplex.cpp tds.cpp)
