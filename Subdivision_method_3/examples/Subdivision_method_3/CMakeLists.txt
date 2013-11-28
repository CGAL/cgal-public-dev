include(CGAL_suite)
CGAL_example_suite(
  NAME Subdivision_method_3_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES CatmullClark_subdivision.cpp Customized_subdivision.cpp DooSabin_subdivision.cpp Loop_subdivision.cpp)
