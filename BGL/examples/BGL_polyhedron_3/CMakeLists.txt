include(CGAL_suite)
CGAL_example_suite(
  NAME BGL_polyhedron_3_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES distance.cpp kruskal.cpp kruskal_with_stored_id.cpp)
