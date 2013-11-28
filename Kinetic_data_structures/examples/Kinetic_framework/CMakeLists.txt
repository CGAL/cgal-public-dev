include(CGAL_suite)
CGAL_example_suite(
  NAME Kinetic_framework_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES adding_a_certificate.cpp defining_a_simulation_traits.cpp listener.cpp pointer_queue.cpp ref_counted.cpp trivial_kds.cpp)
