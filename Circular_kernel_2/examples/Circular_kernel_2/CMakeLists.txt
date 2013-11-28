include(CGAL_suite)
CGAL_example_suite(
  NAME Circular_kernel_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES functor_has_on_2.cpp intersecting_arcs.cpp)
