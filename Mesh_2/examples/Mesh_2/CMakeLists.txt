include(CGAL_suite)
CGAL_example_suite(
  NAME Mesh_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES conforming.cpp mesh_class.cpp mesh_global.cpp mesh_with_seeds.cpp)
