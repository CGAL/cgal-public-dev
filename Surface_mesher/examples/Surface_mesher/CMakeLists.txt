if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Surface_mesher_examples
  REQUIRED_DEPENDENCIES CGAL CGAL_ImageIO
  OPTIONAL_DEPENDENCIES
  SOURCES mesh_a_3d_gray_image.cpp mesh_an_implicit_function.cpp)
