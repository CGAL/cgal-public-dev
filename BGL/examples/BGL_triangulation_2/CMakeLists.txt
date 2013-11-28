include(CGAL_suite)
CGAL_example_suite(
  NAME BGL_triangulation_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES dijkstra.cpp dijkstra_with_internal_properties.cpp emst.cpp)
