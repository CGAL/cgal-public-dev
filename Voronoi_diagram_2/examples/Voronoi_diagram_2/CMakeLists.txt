include(CGAL_suite)
CGAL_example_suite(
  NAME Voronoi_diagram_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES vd_2_point_location.cpp)
