include(CGAL_suite)
CGAL_example_suite(
  NAME Min_annulus_d_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES min_annulus_d.cpp min_annulus_d_fast_exact.cpp)
