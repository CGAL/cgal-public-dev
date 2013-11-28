include(CGAL_suite)
CGAL_example_suite(
  NAME Nef_3_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES comparison.cpp complex_construction.cpp exploration_SM.cpp extended_kernel.cpp handling_double_coordinates.cpp interface_polyhedron.cpp nef_3_construction.cpp nef_3_point_location.cpp nef_3_simple.cpp nefIO.cpp offIO.cpp point_set_operations.cpp polygon_construction.cpp polyline_construction.cpp set_operations.cpp shell_exploration.cpp topological_operations.cpp transformation.cpp)
