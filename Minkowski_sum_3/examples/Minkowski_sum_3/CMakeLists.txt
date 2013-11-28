include(CGAL_suite)
CGAL_example_suite(
  NAME Minkowski_sum_3_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES cube_offset.cpp glide.cpp)
