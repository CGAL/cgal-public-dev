include(CGAL_suite)
CGAL_example_suite(
  NAME Polygon_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES Example.cpp polygon_algorithms.cpp Polygon.cpp projected_polygon.cpp)
