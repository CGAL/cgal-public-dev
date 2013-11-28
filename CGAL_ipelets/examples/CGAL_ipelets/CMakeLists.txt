include(CGAL_suite)
CGAL_example_suite(
  NAME CGAL_ipelets_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES test_grabbers.cpp)
