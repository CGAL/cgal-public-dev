include(CGAL_suite)
CGAL_example_suite(
  NAME Periodic_3_triangulation_3_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES colored_vertices.cpp covering.cpp find_conflicts.cpp geometric_access.cpp large_point_set.cpp periodic_adding_handles.cpp simple_example.cpp)
