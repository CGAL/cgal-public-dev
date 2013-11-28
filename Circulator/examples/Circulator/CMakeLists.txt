include(CGAL_suite)
CGAL_example_suite(
  NAME Circulator_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES circulator_prog1.cpp circulator_prog2.cpp circulator_prog3.cpp)
