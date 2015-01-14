#include <cmath>
#include <cassert>
#include <limits>
#include <vector>
#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/properties/triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

using namespace CGAL::Properties::Triangulation_2;

bool approx_eq(double a, double b) {
  static double epsilon = 1e-3;
  std::cout << "Got: " << a << ", expected: " << b << "." << std::endl;
  return std::fabs(a - b) <= epsilon;
}

template <typename T>
bool exact_eq(T a, T b) {
  std::cout << "Got: " << a << ", expected: " << b << std::endl;
  return a == b;
}

void run_tests_epic() {
  std::cout << "[Testing with EPIC kernel]" << std::endl;

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_2 Point;
  typedef CGAL::Delaunay_triangulation_2<Kernel> Delaunay;
  typedef Delaunay::Vertex_handle Vertex_handle;
  typedef Delaunay::Face_handle Face_handle;
  typedef Delaunay::Edge Edge;

  double infinity = std::numeric_limits<double>::infinity();

  // Triangulations for testing. We consider three test cases,
  // an empty triangulation, a triangulation with one point, and a
  // triangulation with five points.
  Delaunay dt_0;
  Delaunay dt_1;
  Delaunay dt_2;
  Delaunay dt_5;

  // Create the test triangulations, and give short-hand names for
  // the vertices of interest.

  // 0 Dimensional examples.
  Vertex_handle v_0 = dt_0.infinite_vertex();
  Vertex_handle v_1 = dt_1.insert(Point(0.0, 0.0));

  // 1 Dimensional example.
  Vertex_handle v_2_1 = dt_2.insert(Point(0.0, 0.0));
  Vertex_handle v_2_2 = dt_2.insert(Point(1.0, 0.0));

  // 2 Dimensional example.
  std::vector<Vertex_handle> v_5(6);
  v_5[0] = dt_5.insert(Point(0.0, 0.0));
  v_5[1] = dt_5.insert(Point(-0.5, 0.5));
  v_5[2] = dt_5.insert(Point(1.0, 1.0));
  v_5[3] = dt_5.insert(Point(1.0, -2.0));
  v_5[4] = dt_5.insert(Point(-1.0, -2.0));
  v_5[5] = dt_5.infinite_vertex();

  // Check helper functions compile as expected.
  make_degree(dt_0);

  make_dual_area(dt_0);
  make_dual_area(dt_0, CGAL::Finite_test_tag());
  make_dual_area(dt_0, CGAL::No_finite_test_tag());

  make_link_length(dt_0);
  make_link_length(dt_0, CGAL::Finite_test_tag());
  make_link_length(dt_0, CGAL::No_finite_test_tag());

  make_max_star_angle(dt_0);
  make_max_star_angle(dt_0, CGAL::Finite_test_tag());
  make_max_star_angle(dt_0, CGAL::No_finite_test_tag());

  make_min_star_angle(dt_0);
  make_min_star_angle(dt_0, CGAL::Finite_test_tag());
  make_min_star_angle(dt_0, CGAL::No_finite_test_tag());

  //-- Degree ----------------------------------------------------------------//

  std::cout << "Testing 'Degree'" << std::endl;

  Degree<Delaunay> degree;

  assert(exact_eq<int>(degree(v_0), 0));
  assert(exact_eq<int>(degree(v_1), 0));
  assert(exact_eq<int>(degree(v_5[0]), 4));
  assert(exact_eq<int>(degree(v_5[1]), 4));
  assert(exact_eq<int>(degree(v_5[5]), 4));

  //-- Dual_area -------------------------------------------------------------//

  std::cout << "Testing 'Dual_area'" << std::endl;

  Dual_area<Delaunay, CGAL::No_finite_test_tag> dual_area_0_a(dt_0),
      dual_area_1_a(dt_1), dual_area_2_a(dt_2), dual_area_5_a(dt_5);
  Dual_area<Delaunay, CGAL::Finite_test_tag> dual_area_0_b(dt_0),
      dual_area_1_b(dt_1), dual_area_2_b(dt_2), dual_area_5_b(dt_5);

  Dual_area<Delaunay> dual_area_0_c(dt_0), dual_area_1_c(dt_1),
      dual_area_2_c(dt_2), dual_area_5_c(dt_5);

  // For lower-dimensional triangulations, the dual area is zero.
  assert(approx_eq(dual_area_0_a(v_0),   0));
  assert(approx_eq(dual_area_0_b(v_0),   0));
  assert(approx_eq(dual_area_0_c(v_0),   0));
  assert(approx_eq(dual_area_1_a(v_1),   0));
  assert(approx_eq(dual_area_1_b(v_1),   0));
  assert(approx_eq(dual_area_1_c(v_1),   0));
  assert(approx_eq(dual_area_2_a(v_2_1), 0));
  assert(approx_eq(dual_area_2_b(v_2_1), 0));
  assert(approx_eq(dual_area_2_c(v_2_1), 0));

  // The Voronoi cell is unbounded in this these cases.
  assert(exact_eq<double>(dual_area_5_c(v_5[1]), infinity));
  assert(exact_eq<double>(dual_area_5_c(v_5[5]), infinity));

  // Compare dual area with manually computed value.
  assert(approx_eq(dual_area_5_a(v_5[0]), 2.64583));
  assert(approx_eq(dual_area_5_b(v_5[0]), 2.64583));

  //-- Star_area -------------------------------------------------------------//

  std::cout << "Testing 'Star_area'" << std::endl;

  Star_area<Delaunay, CGAL::No_finite_test_tag> star_area_0_a(dt_0),
      star_area_1_a(dt_1), star_area_2_a(dt_2), star_area_5_a(dt_5);
  Star_area<Delaunay, CGAL::Finite_test_tag> star_area_0_b(dt_0),
      star_area_1_b(dt_1), star_area_2_b(dt_2), star_area_5_b(dt_5);
  Star_area<Delaunay> star_area_0_c(dt_0), star_area_1_c(dt_1),
      star_area_2_c(dt_2), star_area_5_c(dt_5);


  assert(approx_eq(star_area_0_a(v_0), 0));
  assert(approx_eq(star_area_0_b(v_0), 0));
  assert(approx_eq(star_area_1_a(v_1), 0));
  assert(approx_eq(star_area_1_b(v_1), 0));
  assert(approx_eq(star_area_2_a(v_2_1), 0));
  assert(approx_eq(star_area_2_b(v_2_1), 0));

  // Unbounded tests.
  assert(exact_eq(star_area_5_b(v_5[2]), infinity));

  // Bounded tests.
  assert(approx_eq(star_area_5_b(v_5[0]), 4.75));
  assert(approx_eq(star_area_5_a(v_5[0]), 4.75));

  //-- Link_length -----------------------------------------------------------//

  std::cout << "Testing 'Link_length'" << std::endl;

  Link_length<Delaunay, CGAL::No_finite_test_tag> link_length_0_a(dt_0),
      link_length_1_a(dt_1), link_length_2_a(dt_2), link_length_5_a(dt_5);
  Link_length<Delaunay, CGAL::Finite_test_tag> link_length_0_b(dt_0),
      link_length_1_b(dt_1), link_length_2_b(dt_2), link_length_5_b(dt_5);
  Link_length<Delaunay> link_length_0_c(dt_0), link_length_1_c(dt_1),
      link_length_2_c(dt_2), link_length_5_c(dt_5);

  // The link length will give zero for degenerate triangulations.
  assert(approx_eq(link_length_0_a(v_0), 0));
  assert(approx_eq(link_length_0_b(v_0), 0));

  assert(approx_eq(link_length_1_a(v_1), 0));
  assert(approx_eq(link_length_1_b(v_1), 0));

  assert(approx_eq(link_length_2_a(v_2_1), 0));
  assert(approx_eq(link_length_2_b(v_2_1), 0));

  // The link length is infinite if the link contains the infinite vertex.
  assert(exact_eq<double>(link_length_5_b(v_5[1]), infinity));

  // The link length of the infinite vertex is the same as the perimeter
  // of the convex hull in this case.
  assert(approx_eq(link_length_5_b(v_5[0]), 9.130648586));
  assert(approx_eq(link_length_5_a(v_5[5]), 9.130648586));

  //-- Max_star_angle --------------------------------------------------------//

  std::cout << "Testing 'Max_star_angle'" << std::endl;

  Max_star_angle<Delaunay, CGAL::No_finite_test_tag> max_star_angle_0_a(dt_0),
      max_star_angle_1_a(dt_1), max_star_angle_2_a(dt_2),
      max_star_angle_5_a(dt_5);
  Max_star_angle<Delaunay, CGAL::Finite_test_tag> max_star_angle_0_b(dt_0),
      max_star_angle_1_b(dt_1), max_star_angle_2_b(dt_2),
      max_star_angle_5_b(dt_5);
  Max_star_angle<Delaunay> max_star_angle_0_c(dt_0), max_star_angle_2_c(dt_2),
      max_star_angle_1_c(dt_1), max_star_angle_5_c(dt_5);

  assert(approx_eq(max_star_angle_0_a(v_0), 0));
  assert(approx_eq(max_star_angle_0_b(v_0), 0));
  assert(approx_eq(max_star_angle_0_c(v_0), 0));

  assert(approx_eq(max_star_angle_1_a(v_1), 0));
  assert(approx_eq(max_star_angle_1_b(v_1), 0));
  assert(approx_eq(max_star_angle_1_c(v_1), 0));

  assert(approx_eq(max_star_angle_2_a(v_2_1), 0));
  assert(approx_eq(max_star_angle_2_b(v_2_1), 0));
  assert(approx_eq(max_star_angle_2_c(v_2_1), 0));

  // Compare with manually computed values.
  assert(approx_eq(max_star_angle_5_a(v_5[0]), 1.89254));
  assert(approx_eq(max_star_angle_5_b(v_5[0]), 1.89254));

  //-- Min_star_angle --------------------------------------------------------//

  std::cout << "Testing 'Min_star_angle'" << std::endl;

  Min_star_angle<Delaunay, CGAL::No_finite_test_tag> min_star_angle_0_a(dt_0),
      min_star_angle_1_a(dt_1), min_star_angle_5_a(dt_5);
  Min_star_angle<Delaunay, CGAL::Finite_test_tag> min_star_angle_0_b(dt_0),
      min_star_angle_1_b(dt_1), min_star_angle_5_b(dt_5);
  Min_star_angle<Delaunay> min_star_angle_0_c(dt_0), min_star_angle_1_c(dt_1),
      min_star_angle_5_c(dt_5);

  // We have chosen to treat an infinite triangle as having an angle of zero
  // adjacent to the infinite point.
  assert(approx_eq(min_star_angle_0_a(v_0), 0));      
  assert(approx_eq(min_star_angle_0_b(v_0), 0));
  assert(approx_eq(min_star_angle_1_a(v_1), 0));
  assert(approx_eq(min_star_angle_1_c(v_1), 0));

  // Compare with manually computed values.
  assert(approx_eq(min_star_angle_5_a(v_5[0]), 0.92729));
  assert(approx_eq(min_star_angle_5_b(v_5[0]), 0.92729));
}

void run_tests_epec() {
  std::cout << "[Testing with EPEC kernel]" << std::endl;

  typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
  typedef Kernel::FT FT;
  typedef Kernel::Point_2 Point;
  typedef CGAL::Delaunay_triangulation_2<Kernel> Delaunay;
  typedef Delaunay::Vertex_handle Vertex_handle;
  typedef Delaunay::Face_handle Face_handle;
  typedef Delaunay::Edge Edge;

  Delaunay dt_0;
  Delaunay dt_1;
  Delaunay dt_2;
  Delaunay dt_5;

  // Create the test triangulations, and give short-hand names for
  // the vertices of interest.
  Vertex_handle v_0 = dt_0.infinite_vertex();
  Vertex_handle v_1 = dt_1.insert(Point(0.0, 0.0));

  // 1 Dimensional example.
  Vertex_handle v_2_1 = dt_2.insert(Point(0.0, 0.0));
  Vertex_handle v_2_2 = dt_2.insert(Point(1.0, 0.0));

  // All vertices in the triangulation with five points.
  std::vector<Vertex_handle> v_5(6);
  v_5[0] = dt_5.insert(Point(0.0, 0.0));
  v_5[1] = dt_5.insert(Point(-0.5, 0.5));
  v_5[2] = dt_5.insert(Point(1.0, 1.0));
  v_5[3] = dt_5.insert(Point(1.0, -2.0));
  v_5[4] = dt_5.insert(Point(-1.0, -2.0));
  v_5[5] = dt_5.infinite_vertex();

  // Check helper functions compile as expected.
  make_dual_area(dt_5, CGAL::No_finite_test_tag());
  make_link_length(dt_5, CGAL::No_finite_test_tag());
  make_max_star_angle(dt_5, CGAL::No_finite_test_tag());
  make_min_star_angle(dt_5, CGAL::No_finite_test_tag());

  //-- Degree ----------------------------------------------------------------//

  std::cout << "Testing 'Degree'" << std::endl;

  Degree<Delaunay> degree;

  assert(exact_eq<int>(degree(v_1), 0));
  assert(exact_eq<int>(degree(v_5[0]), 4));
  assert(exact_eq<int>(degree(v_5[1]), 4));
  assert(exact_eq<int>(degree(v_5[5]), 4));

  //-- Dual_area -------------------------------------------------------------//

  std::cout << "Testing 'Dual_area'" << std::endl;

  Dual_area<Delaunay, CGAL::No_finite_test_tag> dual_area_0(dt_0),
      dual_area_1(dt_1), dual_area_2(dt_2), dual_area_5(dt_5);

  // For lower-dimensional triangulations, the dual area is zero.
  assert(exact_eq<FT>(dual_area_0(v_0),   0));
  assert(exact_eq<FT>(dual_area_1(v_1),   0));
  assert(exact_eq<FT>(dual_area_2(v_2_1), 0));

  // Compare dual area with manually computed value.
  assert(approx_eq(CGAL::to_double(dual_area_5(v_5[0])), 2.64583));

  //-- Star_area -------------------------------------------------------------//

  std::cout << "Testing 'Star_area'" << std::endl;

  Star_area<Delaunay, CGAL::No_finite_test_tag> star_area_5_a(dt_5);

  // Bounded tests.
  assert(exact_eq<FT>(star_area_5_a(v_5[0]), 4.75));

  //-- Link_length -----------------------------------------------------------//

  std::cout << "Testing 'Link_length'" << std::endl;

  Link_length<Delaunay, CGAL::No_finite_test_tag> link_length_1_a(dt_1),
      link_length_5_a(dt_5);

  // The link length around the single vertex is zero.
  assert(exact_eq<FT>(link_length_1_a(v_1), 0));

  // The link length of the infinite vertex is the same as the perimeter
  // of the convex hull in this case.
  assert(approx_eq(CGAL::to_double(link_length_5_a(v_5[5])), 9.130648586));

  //-- Max_star_angle --------------------------------------------------------//

  std::cout << "Testing 'Max_star_angle'" << std::endl;

  Max_star_angle<Delaunay, CGAL::No_finite_test_tag> max_star_angle_5_a(dt_5);

  // Compare with manually computed values.
  assert(approx_eq(CGAL::to_double(max_star_angle_5_a(v_5[0])), 1.89254));

  //-- Min_star_angle --------------------------------------------------------//

  std::cout << "Testing 'Min_star_angle'" << std::endl;

  Min_star_angle<Delaunay, CGAL::No_finite_test_tag> min_star_angle_5_a(dt_5);

  // Compare with manually computed values.
  assert(approx_eq(CGAL::to_double(min_star_angle_5_a(v_5[0])), 0.92729));
}

int main() {
  std::cout << "[Testing Triangulation_2 vertex properties]" << std::endl;

  run_tests_epic();
  run_tests_epec();

  std::cout << "[Sucess]" << std::endl;
}
