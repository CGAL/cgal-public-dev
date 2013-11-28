include(CGAL_suite)
CGAL_example_suite(
  NAME Convex_hull_3_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES dynamic_hull_3.cpp incremental_hull_3.cpp incremental_hull_class_3.cpp quickhull_3.cpp)
