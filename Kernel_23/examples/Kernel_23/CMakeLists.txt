include(CGAL_suite)
CGAL_example_suite(
  NAME Kernel_23_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES cartesian_converter.cpp intersections.cpp MyKernel.cpp)
