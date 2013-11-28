include(CGAL_suite)
CGAL_example_suite(
  NAME Kinetic_data_structures_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES Kinetic_Delaunay_triangulation_2.cpp Kinetic_Delaunay_triangulation_3.cpp Kinetic_insert_event.cpp Kinetic_regular_triangulation_3.cpp Kinetic_sort.cpp Kinetic_sweepline.cpp)
