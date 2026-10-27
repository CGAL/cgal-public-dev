#pragma once
#include "common.h"

bool run_intersect_2d_test(Point_2 p, Point_2 q, Point_2 r, Point_2 s, bool expected_with_overlap, bool expected_without_overlap) {
  Polyhedron mesh;

  double edge_length = 0.5;

  // const FT padding = 1e-8;
  // Build the grid and insert points + normals.
  CGAL::Box_grid<Kernel> grid{FT(edge_length),
                              FT(-1.0), FT(-1.0), FT(-1.0),
                              FT(3.0), FT(3.0), FT(3.0)};

  // Construct the reconstruction after initializing the mesh.
  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, std::vector<Vector_3>{});

  TEST_ASSERT(expected_without_overlap == Access::intersect_2d(reconstruction, p, q, r, s, false), "Without overlap incorrect");
  TEST_ASSERT(expected_with_overlap == Access::intersect_2d(reconstruction, p, q, r, s, true), "With overlap incorrect");

  return true;
}

bool test_intersect_2d_intersecting() {
  Point_2 p = Point_2(0, 0);
  Point_2 q = Point_2(1, 1);
  Point_2 r = Point_2(0, 1);
  Point_2 s = Point_2(1, 0);

  return run_intersect_2d_test(p, q, r, s, true, true);
}

bool test_intersect_2d_non_intersecting_parallel() {
  Point_2 p = Point_2(0, 0);
  Point_2 q = Point_2(0, 1);
  Point_2 r = Point_2(1, 0);
  Point_2 s = Point_2(1, 1);

  return run_intersect_2d_test(p, q, r, s, false, false);
}

bool test_intersect_2d_close_parallel_lines() {
  Point_2 p = Point_2(0, 0);
  Point_2 q = Point_2(0, 1);
  Point_2 r = Point_2(0.01, 0);
  Point_2 s = Point_2(0.01, 1);

  return run_intersect_2d_test(p, q, r, s, false, false);
}

bool test_intersect_2d_non_intersecting() {
  Point_2 p = Point_2(0, 0);
  Point_2 q = Point_2(0, 1);
  Point_2 r = Point_2(0.5, 0);
  Point_2 s = Point_2(0.3, 2);

  return run_intersect_2d_test(p, q, r, s, false, false);
}

bool test_intersect_2d_same() {
  Point_2 p = Point_2(0, 0);
  Point_2 q = Point_2(0, 1);
  Point_2 r = Point_2(0, 0);
  Point_2 s = Point_2(0, 1);

  return run_intersect_2d_test(p, q, r, s, false, true);
}

bool test_intersect_2d_shared_endpoint() {
  Point_2 p = Point_2(0, 0);
  Point_2 q = Point_2(0, 1);
  Point_2 r = Point_2(0, 0);
  Point_2 s = Point_2(1, 0);

  return run_intersect_2d_test(p, q, r, s, false, false);
}

bool test_intersect_2d_line_overlap() {
  Point_2 p = Point_2(0, 0);
  Point_2 q = Point_2(0, 1);
  Point_2 r = Point_2(0, 0.5);
  Point_2 s = Point_2(0, 1.5);

  return run_intersect_2d_test(p, q, r, s, false, true);
}

bool test_intersect_2d_endpoint_on_interior() {
  Point_2 p = Point_2(0, 0);
  Point_2 q = Point_2(0, 1);
  Point_2 r = Point_2(0, 0.5);
  Point_2 s = Point_2(1, 1);

  return run_intersect_2d_test(p, q, r, s, false, false);
}

void run_intersect_2d_tests(TestData& test_data) {
  test(test_intersect_2d_intersecting, "test_intersect_2d_intersecting", test_data);
  test(test_intersect_2d_non_intersecting_parallel, "test_intersect_2d_non_intersecting_parallel", test_data);
  test(test_intersect_2d_close_parallel_lines, "test_intersect_2d_close_parallel_lines", test_data);
  test(test_intersect_2d_non_intersecting, "test_intersect_2d_non_intersecting", test_data);
  test(test_intersect_2d_same, "test_intersect_2d_same", test_data);
  test(test_intersect_2d_shared_endpoint, "test_intersect_2d_shared_endpoint", test_data);
  test(test_intersect_2d_line_overlap, "test_intersect_2d_line_overlap", test_data);
  test(test_intersect_2d_endpoint_on_interior, "test_intersect_2d_endpoint_on_interior", test_data);
}