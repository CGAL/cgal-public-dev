include(CGAL_suite)
CGAL_example_suite(
  NAME Alpha_shapes_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES ex_alpha_shapes_2.cpp ex_weighted_alpha_shapes_2.cpp)
