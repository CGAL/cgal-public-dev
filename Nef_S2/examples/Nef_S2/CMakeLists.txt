include(CGAL_suite)
CGAL_example_suite(
  NAME Nef_S2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES nef_s2_construction.cpp nef_s2_exploration.cpp nef_s2_point_location.cpp nef_s2_simple.cpp)
