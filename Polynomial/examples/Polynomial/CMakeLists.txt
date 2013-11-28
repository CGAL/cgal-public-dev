include(CGAL_suite)
CGAL_example_suite(
  NAME Polynomial_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES coefficient_access.cpp construction.cpp degree.cpp gcd_up_to_constant_factor.cpp subresultants.cpp substitute.cpp swap_move.cpp)
