include(CGAL_suite)
CGAL_example_suite(
  NAME Profiling_tools_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES Profile_counter.cpp Profile_histogram_counter.cpp Profile_timer.cpp)
