include(CGAL_suite)
CGAL_example_suite(
  NAME Spatial_sorting_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES example_delaunay_2.cpp hilbert.cpp hilbert_policies.cpp myPoint.cpp small_example_delaunay_2.cpp sort_indices.cpp sp_sort_using_property_map_2.cpp sp_sort_using_property_map_3.cpp sp_sort_using_property_map_d.cpp)
