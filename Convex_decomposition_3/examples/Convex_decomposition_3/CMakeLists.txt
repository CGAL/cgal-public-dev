include(CGAL_suite)
CGAL_example_suite(
  NAME Convex_decomposition_3_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES list_of_convex_parts.cpp)
