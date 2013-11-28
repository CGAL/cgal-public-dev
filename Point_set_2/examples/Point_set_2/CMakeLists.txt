include(CGAL_suite)
CGAL_example_suite(
  NAME Point_set_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES nearest_neighbor.cpp range_search.cpp)
