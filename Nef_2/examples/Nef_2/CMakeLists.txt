include(CGAL_suite)
CGAL_example_suite(
  NAME Nef_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES nef_2_construction.cpp nef_2_intersection.cpp nef_2_polylines.cpp)
