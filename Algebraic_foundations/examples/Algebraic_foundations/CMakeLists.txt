include(CGAL_suite)
CGAL_example_suite(
  NAME Algebraic_foundations_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES algebraic_structure_dispatch.cpp fraction_traits.cpp implicit_interoperable_dispatch.cpp integralize.cpp interoperable.cpp)
