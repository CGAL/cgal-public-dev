include(CGAL_suite)
CGAL_example_suite(
  NAME Box_intersection_d_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES box_grid.cpp custom_box_grid.cpp minimal.cpp minimal_self.cpp proximity_custom_box_traits.cpp triangle_self_intersect.cpp triangle_self_intersect_pointers.cpp)
