include(CGAL_suite)
CGAL_example_suite(
  NAME Convex_hull_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES array_convex_hull_2.cpp ch_from_cin_to_cout.cpp ch_graham_anderson.cpp ch_timing.cpp convex_hull_yz.cpp iostream_convex_hull_2.cpp vector_convex_hull_2.cpp)
