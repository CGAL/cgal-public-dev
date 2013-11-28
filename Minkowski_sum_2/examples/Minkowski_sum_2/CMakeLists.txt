include(CGAL_suite)
CGAL_example_suite(
  NAME Minkowski_sum_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES approx_inset.cpp approx_offset.cpp exact_inset.cpp exact_offset.cpp sum_by_decomposition.cpp sum_triangles.cpp sum_with_holes.cpp)
