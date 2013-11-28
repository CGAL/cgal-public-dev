include(CGAL_suite)
CGAL_example_suite(
  NAME Partition_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES approx_convex_partition_2.cpp greene_approx_convex_partition_2.cpp optimal_convex_partition_2.cpp y_monotone_partition_2.cpp)
