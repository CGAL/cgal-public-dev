if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME QP_solver_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES convex_hull_containment2.cpp convex_hull_containment_benchmarks.cpp convex_hull_containment.cpp cycling.cpp first_lp.cpp first_lp_from_iterators.cpp first_lp_from_mps.cpp first_nonnegative_lp.cpp first_nonnegative_lp_from_iterators.cpp first_nonnegative_lp_from_mps.cpp first_nonnegative_qp.cpp first_nonnegative_qp_from_iterators.cpp first_nonnegative_qp_from_mps.cpp first_qp_basic_constraints.cpp first_qp.cpp first_qp_from_iterators.cpp first_qp_from_mps.cpp important_variables.cpp infeasibility_certificate.cpp invert_matrix.cpp optimality_certificate.cpp print_first_lp.cpp print_first_nonnegative_lp.cpp print_first_nonnegative_qp.cpp print_first_qp.cpp unboundedness_certificate.cpp)
