include(CGAL_suite)
CGAL_example_suite(
  NAME Algebraic_kernel_d_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES CGAL_Core MPFI
  SOURCES Compare_1.cpp Construct_algebraic_real_1.cpp Isolate_1.cpp Sign_at_1.cpp Solve_1.cpp)
