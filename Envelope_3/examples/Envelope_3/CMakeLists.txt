include(CGAL_suite)
CGAL_example_suite(
  NAME Envelope_3_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES envelope_planes.cpp envelope_spheres.cpp envelope_triangles.cpp)
