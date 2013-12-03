if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Triangulation_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES adding_handles.cpp colored_face.cpp constrained.cpp constrained_hierarchy_plus.cpp constrained_plus.cpp copy_triangulation_2.cpp hierarchy.cpp info_insert_with_pair_iterator_2.cpp info_insert_with_pair_iterator_regular_2.cpp info_insert_with_transform_iterator_2.cpp info_insert_with_zip_iterator_2.cpp polygon_triangulation.cpp print_cropped_voronoi.cpp regular.cpp terrain.cpp triangulation_prog1.cpp voronoi.cpp)
