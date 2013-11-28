include(CGAL_suite)
CGAL_example_suite(
  NAME Surface_mesh_simplification_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES edge_collapse_constrained_polyhedron.cpp edge_collapse_enriched_polyhedron.cpp edge_collapse_polyhedron.cpp)
