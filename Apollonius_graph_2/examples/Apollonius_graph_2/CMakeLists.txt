include(CGAL_suite)
CGAL_example_suite(
  NAME Apollonius_graph_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES ag2_exact_traits.cpp ag2_exact_traits_sqrt.cpp ag2_filtered_traits_no_hidden.cpp ag2_hierarchy.cpp)
