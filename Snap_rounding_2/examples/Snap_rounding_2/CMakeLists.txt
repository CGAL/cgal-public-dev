include(CGAL_suite)
CGAL_example_suite(
  NAME Snap_rounding_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES snap_rounding.cpp snap_rounding_to_integer.cpp)
