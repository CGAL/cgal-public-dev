include(CGAL_suite)
CGAL_example_suite(
  NAME Stream_lines_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES stl_regular_field.cpp stl_triangular_field.cpp)
