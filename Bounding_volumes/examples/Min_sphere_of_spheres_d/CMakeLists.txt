include(CGAL_suite)
CGAL_example_suite(
  NAME Min_sphere_of_spheres_d_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES benchmark.cpp min_sphere_of_spheres_d_2.cpp min_sphere_of_spheres_d_3.cpp min_sphere_of_spheres_d_d.cpp)
