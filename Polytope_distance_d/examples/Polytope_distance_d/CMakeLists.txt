include(CGAL_suite)
CGAL_example_suite(
  NAME Polytope_distance_d_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES all_furthest_neighbors_2.cpp polytope_distance_d.cpp polytope_distance_d_fast_exact.cpp width_simplex.cpp)
