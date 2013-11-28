include(CGAL_suite)
CGAL_example_suite(
  NAME Core_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES delaunay.cpp)
