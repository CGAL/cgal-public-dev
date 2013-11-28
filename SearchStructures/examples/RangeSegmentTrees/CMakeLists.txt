include(CGAL_suite)
CGAL_example_suite(
  NAME RangeSegmentTrees_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES range_tree_1.cpp range_tree_2.cpp range_tree_3.cpp range_tree_4.cpp range_tree_map_2.cpp range_tree_set_2.cpp segment_tree_1.cpp segment_tree_2.cpp segment_tree_3.cpp segment_tree_4.cpp segment_tree_map_2.cpp segment_tree_set_2.cpp segment_tree_set_3.cpp software_design_rgs.cpp)
