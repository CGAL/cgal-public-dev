include(CGAL_suite)
CGAL_example_suite(
  NAME Inscribed_areas_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES extremal_polygon_2_area.cpp extremal_polygon_2_perimeter.cpp largest_empty_rectangle.cpp)
