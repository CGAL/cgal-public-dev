include(CGAL_suite)
CGAL_example_suite(
  NAME Approximate_min_ellipsoid_d_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES ellipsoid.cpp ellipsoid_for_maple.cpp)
