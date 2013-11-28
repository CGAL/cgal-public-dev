include(CGAL_suite)
CGAL_example_suite(
  NAME Principal_component_analysis_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES barycenter.cpp bounding_box.cpp centroid.cpp linear_least_squares_fitting_points_2.cpp linear_least_squares_fitting_triangles_3.cpp)
