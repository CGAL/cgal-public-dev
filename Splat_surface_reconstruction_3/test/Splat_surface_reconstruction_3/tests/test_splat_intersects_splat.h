#pragma once
#include "common.h"

inline TestGrid make_test_grid(FT d, std::vector<Point> points, std::vector<Vector_3> normals) {
  TestGrid grid = GridTestAccess::make_grid(d, -2.0, -2.0, -2.0, 2.0, 2.0, 2.0);
  grid.build(points, normals);
  return grid;
}

inline bool test_splat_intersects_splat_edge_on_edge() {
  std::vector<Point> points = {
    Point(0.0, 0.0, 0.0),
    Point(0.0, 0.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(1.0, 0.0, 0.0),
    Vector_3(0.0, 1.0, 0.0)
  };
  TestGrid grid = make_test_grid(1.0, points, normals);
  grid.set_all_splat_sizes(1.0);

  std::vector<FT> splat_sizes = GridTestAccess::splat_sizes(grid);

  Polyhedron mesh;

  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, std::vector<Vector_3>{});

  Point_3 intersection_point_1, intersection_point_2;
  bool intersects = TestAccess::splat_intersects_splat(reconstruction, 0, 1, intersection_point_1, intersection_point_2);

  std::vector<Point_3> intersection_points = {
    intersection_point_1,
    intersection_point_2
  };

  std::vector<Point_3> expected_intersection_points = {
    Point_3(0.0, 0.0, -1.0),
    Point_3(0.0, 0.0, 1.0)
  };

  for (std::size_t i = 0; i < intersection_points.size(); ++i) {
    for (std::size_t j = 0; j < expected_intersection_points.size(); ++j) {
      if (approxEq(intersection_points[i], expected_intersection_points[j])) {
        expected_intersection_points.erase(expected_intersection_points.begin() + j);
        break;
      }
    }
  }

  TEST_ASSERT(expected_intersection_points.empty(), "All expected intersection points should be found");


  return true;
}

inline bool test_splat_intersects_splat_part_overlap() {
  std::vector<Point> points = {
    Point(0.0, 0.0, -0.5),
    Point(0.0, 0.0, 0.5)
  };
  std::vector<Vector_3> normals = {
    Vector_3(1.0, 0.0, 0.0),
    Vector_3(0.0, 1.0, 0.0)
  };
  TestGrid grid = make_test_grid(1.0, points, normals);
  grid.set_all_splat_sizes(1.0);

  std::vector<FT> splat_sizes = GridTestAccess::splat_sizes(grid);

  Polyhedron mesh;

  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, std::vector<Vector_3>{});

  Point_3 intersection_point_1, intersection_point_2;
  bool intersects = TestAccess::splat_intersects_splat(reconstruction, 0, 1, intersection_point_1, intersection_point_2);

  std::vector<Point_3> intersection_points = {
    intersection_point_1,
    intersection_point_2
  };

  std::vector<Point_3> expected_intersection_points = {
    Point_3(0.0, 0.0, -0.5),
    Point_3(0.0, 0.0, 0.5)
  };

  for (std::size_t i = 0; i < intersection_points.size(); ++i) {
    for (std::size_t j = 0; j < expected_intersection_points.size(); ++j) {
      if (approxEq(intersection_points[i], expected_intersection_points[j])) {
        expected_intersection_points.erase(expected_intersection_points.begin() + j);
        break;
      }
    }
  }
  TEST_ASSERT(expected_intersection_points.empty(), "All expected intersection points should be found");


  return true;
}

inline bool test_splat_intersects_splat_both_intersections_on_one_border() {
  double half_sqrt_3 = 0.5 * std::sqrt(3.0);
  std::vector<Point> points = {
    Point(0.0, 0.0, 0.0),
    Point(half_sqrt_3, 0.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(1.0, 0.0, 0.0),
    Vector_3(0.0, 1.0, 0.0)
  };
  TestGrid grid = make_test_grid(1.0, points, normals);
  grid.set_all_splat_sizes(1.0);

  std::vector<FT> splat_sizes = GridTestAccess::splat_sizes(grid);

  Polyhedron mesh;

  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, std::vector<Vector_3>{});

  Point_3 intersection_point_1, intersection_point_2;
  bool intersects = TestAccess::splat_intersects_splat(reconstruction, 0, 1, intersection_point_1, intersection_point_2);

  std::vector<Point_3> intersection_points = {
    intersection_point_1,
    intersection_point_2
  };

  std::vector<Point_3> expected_intersection_points = {
    Point_3(0.0, 0.0, -0.5),
    Point_3(0.0, 0.0, 0.5)
  };

  for (std::size_t i = 0; i < intersection_points.size(); ++i) {
    for (std::size_t j = 0; j < expected_intersection_points.size(); ++j) {
      if (approxEq(intersection_points[i], expected_intersection_points[j])) {
        expected_intersection_points.erase(expected_intersection_points.begin() + j);
        break;
      }
    }
  }


  TEST_ASSERT(expected_intersection_points.empty(), "All expected intersection points should be found");


  return true;
}

inline bool test_splat_intersects_splat_single_intersection() {
  std::vector<Point> points = {
    Point(0.0, 0.0, 0.0),
    Point(1.0, 0.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(1.0, 0.0, 0.0),
    Vector_3(0.0, 1.0, 0.0)
  };
  TestGrid grid = make_test_grid(1.0, points, normals);
  grid.set_all_splat_sizes(1.0);

  std::vector<FT> splat_sizes = GridTestAccess::splat_sizes(grid);

  Polyhedron mesh;

  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, std::vector<Vector_3>{});

  Point_3 intersection_point_1, intersection_point_2;
  bool intersects = TestAccess::splat_intersects_splat(reconstruction, 0, 1, intersection_point_1, intersection_point_2);

  TEST_ASSERT(!intersects, "Splats should not intersect");

  return true;
}



inline void run_splat_intersects_splats_tests(TestData& test_data) {
  test(test_splat_intersects_splat_edge_on_edge, "test_splat_intersects_splat_edge_on_edge", test_data);
  test(test_splat_intersects_splat_part_overlap, "test_splat_intersects_splat_part_overlap", test_data);
  test(test_splat_intersects_splat_both_intersections_on_one_border, "test_splat_intersects_splat_both_intersections_on_one_border", test_data);
  test(test_splat_intersects_splat_single_intersection, "test_splat_intersects_splat_single_intersection", test_data);
}
