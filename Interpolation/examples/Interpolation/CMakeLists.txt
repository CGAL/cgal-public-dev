include(CGAL_suite)
CGAL_example_suite(
  NAME Interpolation_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES interpolation_2.cpp linear_interpolation_2.cpp nn_coordinates_2.cpp nn_coordinates_3.cpp rn_coordinates_2.cpp sibson_interpolation_2.cpp surface_neighbor_coordinates_3.cpp)
