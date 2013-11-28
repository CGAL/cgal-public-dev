include(CGAL_suite)
CGAL_example_suite(
  NAME Min_quadrilateral_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES minimum_enclosing_parallelogram_2.cpp minimum_enclosing_rectangle_2.cpp minimum_enclosing_strip_2.cpp)
