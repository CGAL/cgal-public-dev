include(CGAL_suite)
CGAL_example_suite(
  NAME Circular_kernel_3_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES functor_compare_theta_3.cpp functor_has_on_3.cpp intersecting_spheres.cpp)
