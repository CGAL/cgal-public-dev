include(CGAL_suite)
CGAL_example_suite(
  NAME Interval_skip_list_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES intervals.cpp isl_terrain.cpp)
