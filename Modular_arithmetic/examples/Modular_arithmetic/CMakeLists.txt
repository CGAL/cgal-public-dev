include(CGAL_suite)
CGAL_example_suite(
  NAME Modular_arithmetic_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES modular_filter.cpp)
