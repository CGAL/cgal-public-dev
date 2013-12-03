if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

# If you want to visualize a linear cell complex, there are 2 viewers 
# based on qt and vtk. Just uncomment the corresponding line
# include("CMakeLCCViewerQt.inc")
# include("CMakeLCCViewerVtk.inc")

CGAL_example_suite(
  NAME Linear_cell_complex_examples
  REQUIRED_DEPENDENCIES CGAL CGAL_Qt4
  OPTIONAL_DEPENDENCIES
  SOURCES linear_cell_complex_3.cpp linear_cell_complex_4.cpp
  linear_cell_complex_3_with_colored_vertices.cpp linear_cell_complex_3_with_mypoint.cpp
  plane_graph_to_lcc_2.cpp linear_cell_complex_3_triangulation.cpp voronoi_2.cpp voronoi_3.cpp)

# Only add additional link targets, if we those targets are actually
# built.
if(TARGET linear_cell_complex_3_triangulation)
  target_link_libraries(linear_cell_complex_3_triangulation ${MAP_VIEWER_LIBRARIES})
endif()
if(TARGET voronoi_2)
  target_link_libraries(voronoi_2 ${MAP_VIEWER_LIBRARIES})
endif()
if(TARGET voronoi_3)
  target_link_libraries(voronoi_3 ${MAP_VIEWER_LIBRARIES})
endif()