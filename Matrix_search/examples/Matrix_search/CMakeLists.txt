include(CGAL_suite)
CGAL_example_suite(
  NAME Matrix_search_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES sorted_matrix_search.cpp)
