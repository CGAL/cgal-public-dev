include(CGAL_suite)
CGAL_example_suite(
  NAME Filtered_kernel_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES Filtered_predicate.cpp)
