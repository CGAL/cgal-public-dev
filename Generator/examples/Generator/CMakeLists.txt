include(CGAL_suite)
CGAL_example_suite(
  NAME Generator_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES ball_d.cpp combination_enumerator.cpp cube_d.cpp grid_d.cpp name_pairs.cpp random_convex_set.cpp random_degenerate_point_set.cpp random_grid.cpp random_polygon2.cpp random_polygon.cpp random_segments1.cpp random_segments2.cpp sphere_d.cpp)
