include(CGAL_suite)
CGAL_example_suite(
  NAME STL_Extension_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES Default.cpp Dispatch_output_iterator.cpp in_place_list_prog.cpp min_element_if_example.cpp min_max_element_example.cpp)
