include(CGAL_suite)
CGAL_example_suite(
  NAME Surface_mesher_examples
  REQUIRED_DEPENDENCIES CGAL CGAL_ImageIO
  OPTIONAL_DEPENDENCIES
  SOURCES mesh_a_3d_gray_image.cpp mesh_an_implicit_function.cpp)
