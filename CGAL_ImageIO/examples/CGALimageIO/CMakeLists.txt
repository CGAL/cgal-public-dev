include(CGAL_suite)
CGAL_example_suite(
  NAME CGAL_ImageIO_examples
  REQUIRED_DEPENDENCIES CGAL CGAL_ImageIO
  SOURCES convert_raw_image_to_inr.cpp test_imageio.cpp)
