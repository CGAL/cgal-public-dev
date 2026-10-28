#pragma once
#include "common.h"

inline CGAL::Box_grid<Kernel> make_box_grid(std::vector<Point> points, std::vector<Vector_3> normals, FT d) {

  double edge_length = d;

  const auto bbox = CGAL::bounding_box(points.begin(), points.end()); // recompute bbox after centering and scaling


  const FT padding = FT(1.1) * edge_length;

  // const FT padding = 1e-8;
  // Build the grid and insert points + normals.
  CGAL::Box_grid<Kernel> grid{d,
                              -2.0, -2.0, -2.0,
                              2.0, 2.0, 2.0 }; // initialize grid with cell size equal to average spacing and bounding box [-1,1]^3
  grid.build(points, normals);
  grid.estimate_individual_splat_sizes();

  return grid;
}

inline bool test_generate_candidates_from_parents_basic_case() {
  std::vector<Point_3> splat_points = {
    Point_3(0.5, 0.0, 0.0)
  };
  std::vector<Vector_3> splat_normals = {
    Vector_3(0.0, 0.0, 1.0)
  };

  CGAL::Box_grid<Kernel> grid = make_box_grid(splat_points, splat_normals, FT(1.0));

  std::vector<Point_3> points = {
    Point_3(0.0, 0.0, 0.0),
    Point_3(1.0, 0.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  Polyhedron mesh;

  std::vector<vertex_descriptor> point_descriptors;

  for (const Point& p : points)
    point_descriptors.push_back(mesh.add_vertex(p));


  // Construct the reconstruction after initializing the mesh.
  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, normals);

  std::vector<Candidate> candidates = TestAccess::generate_candidates_from_parents(
    reconstruction, point_descriptors[0], point_descriptors[1]);

  TEST_ASSERT(candidates.size() == 2, "Expected 2 candidates from basic case");
  TEST_ASSERT(approxEq(candidates[0].position, Point_3(0.5, std::sqrt(3.0)/2.0, 0.0)) ||
              approxEq(candidates[0].position, Point_3(0.5, -std::sqrt(3.0)/2.0, 0.0)), "Expected candidate at midpoint");
  TEST_ASSERT(approxEq(candidates[1].position, Point_3(0.5, std::sqrt(3.0)/2.0, 0.0)) ||
              approxEq(candidates[1].position, Point_3(0.5, -std::sqrt(3.0)/2.0, 0.0)), "Expected candidate at midpoint");

  return true;
}

inline bool test_generate_candidates_from_parents_too_close() {
  std::vector<Point_3> splat_points = {
    Point_3(0.5, 0.0, 0.0)
  };
  std::vector<Vector_3> splat_normals = {
    Vector_3(0.0, 0.0, 1.0)
  };

  CGAL::Box_grid<Kernel> grid = make_box_grid(splat_points, splat_normals, FT(1.0));

  std::vector<Point_3> points = {
    Point_3(0.0, 0.0, 0.0),
    Point_3(0.99, 0.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  Polyhedron mesh;

  std::vector<vertex_descriptor> point_descriptors;

  for (const Point& p : points)
    point_descriptors.push_back(mesh.add_vertex(p));


  // Construct the reconstruction after initializing the mesh.
  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, normals);

  std::vector<Candidate> candidates = TestAccess::generate_candidates_from_parents(
    reconstruction, point_descriptors[0], point_descriptors[1]);

  TEST_ASSERT(candidates.size() == 0, "Expected 0 candidates from too close invalid case");

  return true;
}

inline bool test_generate_candidates_from_parents_too_far() {
  std::vector<Point_3> splat_points = {
    Point_3(0.0, 0.0, 0.0)
  };
  std::vector<Vector_3> splat_normals = {
    Vector_3(0.0, 0.0, 1.0)
  };

  CGAL::Box_grid<Kernel> grid = make_box_grid(splat_points, splat_normals, FT(1.0));

  std::vector<Point_3> points = {
    Point_3(-1.0, 0.0, 0.0),
    Point_3(1.01, 0.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  Polyhedron mesh;

  std::vector<vertex_descriptor> point_descriptors;

  for (const Point& p : points)
    point_descriptors.push_back(mesh.add_vertex(p));


  // Construct the reconstruction after initializing the mesh.
  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, normals);

  std::vector<Candidate> candidates = TestAccess::generate_candidates_from_parents(
    reconstruction, point_descriptors[0], point_descriptors[1]);

  TEST_ASSERT(candidates.size() == 0, "Expected 0 candidates from too far invalid case");

  return true;
}

inline void run_generate_candidates_from_parents_tests(TestData& test_data) {
  test(test_generate_candidates_from_parents_basic_case, "test_generate_candidates_from_parents_basic_case", test_data);
  test(test_generate_candidates_from_parents_too_close, "test_generate_candidates_from_parents_too_close", test_data);
  test(test_generate_candidates_from_parents_too_far, "test_generate_candidates_from_parents_too_far", test_data);
}
