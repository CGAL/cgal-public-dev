include(CGAL_suite)
CGAL_example_suite(
  NAME Envelope_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES convex_hull.cpp envelope_circles.cpp envelope_segments.cpp)
