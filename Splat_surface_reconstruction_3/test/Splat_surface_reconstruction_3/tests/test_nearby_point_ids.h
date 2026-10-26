#pragma once
#include "common.h"

inline bool nearby_point_ids_contains_all_points_within_radius(
    TestGrid& grid,
    Point_3 query,
    double radius)
{
  const std::vector<std::size_t> ids =
      GridTestAccess::nearby_point_ids(grid, query, radius);

  const std::vector<Point_3>& points = grid.points();
  const FT radius_squared = FT(radius) * FT(radius);

  for (std::size_t i = 0; i < points.size(); ++i) {
    if (CGAL::squared_distance(query, points[i]) <= radius_squared) {
      if (std::find(ids.begin(), ids.end(), i) == ids.end()) {
        return false;
      }
    }
  }

  return true;
}

inline bool test_nearby_point_ids_small_radius() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0);

  std::vector<Point_3> points = {
    Point_3(0.5, 0.5, 0.5),
    Point_3(1.5, 0.5, 0.5),
    Point_3(2.5, 0.5, 0.5),
    Point_3(3.5, 0.5, 0.5)
  };

  std::vector<Vector_3> normals(
      points.size(), Vector_3(0.0, 0.0, 1.0));

  grid.build(points, normals);

  Point_3 query(1.5, 0.5, 0.5);

  TEST_ASSERT(
      nearby_point_ids_contains_all_points_within_radius(
          grid, query, 0.1),
      "A small-radius query must include every point within its radius");

  return true;
}

inline bool test_nearby_point_ids_neighboring_cells() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0);

  std::vector<Point_3> points = {
    Point_3(0.95, 1.5, 1.5),
    Point_3(1.05, 1.5, 1.5),
    Point_3(2.05, 1.5, 1.5),
    Point_3(3.05, 1.5, 1.5)
  };

  std::vector<Vector_3> normals(
      points.size(), Vector_3(0.0, 0.0, 1.0));

  grid.build(points, normals);

  Point_3 query(1.0, 1.5, 1.5);

  TEST_ASSERT(
      nearby_point_ids_contains_all_points_within_radius(
          grid, query, 1.1),
      "The query must include nearby points across cell boundaries");

  return true;
}

inline bool test_nearby_point_ids_large_radius() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0);

  std::vector<Point_3> points = {
    Point_3(0.5, 0.5, 0.5),
    Point_3(3.5, 0.5, 0.5),
    Point_3(0.5, 3.5, 3.5),
    Point_3(3.5, 3.5, 3.5)
  };

  std::vector<Vector_3> normals(
      points.size(), Vector_3(0.0, 0.0, 1.0));

  grid.build(points, normals);

  Point_3 query(2.0, 2.0, 2.0);

  TEST_ASSERT(
      nearby_point_ids_contains_all_points_within_radius(
          grid, query, 4.0),
      "A large-radius query must include all points within range");

  return true;
}

inline bool test_nearby_point_ids_zero_radius() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0);

  std::vector<Point_3> points = {
    Point_3(1.5, 1.5, 1.5),
    Point_3(1.5, 1.5, 1.5),
    Point_3(2.5, 1.5, 1.5)
  };

  std::vector<Vector_3> normals(
      points.size(), Vector_3(0.0, 0.0, 1.0));

  grid.build(points, normals);

  Point_3 query(1.5, 1.5, 1.5);

  TEST_ASSERT(
      nearby_point_ids_contains_all_points_within_radius(
          grid, query, 0.0),
      "A zero-radius query must include points exactly at the query");

  const std::vector<std::size_t> ids =
      GridTestAccess::nearby_point_ids(grid, query, 0.0);

  TEST_ASSERT(
      std::is_sorted(ids.begin(), ids.end()),
      "Returned point IDs should be sorted");

  TEST_ASSERT(
      std::adjacent_find(ids.begin(), ids.end()) == ids.end(),
      "Returned point IDs should not contain duplicates");

  return true;
}

inline bool test_nearby_point_ids_near_grid_boundary() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0);

  std::vector<Point_3> points = {
    Point_3(0.1, 0.5, 0.5),
    Point_3(0.5, 0.5, 0.5),
    Point_3(1.5, 0.5, 0.5),
    Point_3(3.9, 3.5, 3.5)
  };

  std::vector<Vector_3> normals(
      points.size(), Vector_3(0.0, 0.0, 1.0));

  grid.build(points, normals);

  TEST_ASSERT(
      nearby_point_ids_contains_all_points_within_radius(
          grid, Point_3(0.05, 0.5, 0.5), 0.5),
      "A query near the lower grid boundary must find nearby points");

  TEST_ASSERT(
      nearby_point_ids_contains_all_points_within_radius(
          grid, Point_3(3.95, 3.5, 3.5), 0.2),
      "A query near the upper grid boundary must find nearby points");

  return true;
}

inline bool test_nearby_point_ids_query_outside_grid() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0);

  std::vector<Point_3> points = {
    Point_3(0.2, 1.5, 1.5),
    Point_3(1.5, 1.5, 1.5)
  };

  std::vector<Vector_3> normals(
      points.size(), Vector_3(0.0, 0.0, 1.0));

  grid.build(points, normals);

  Point_3 query(-0.1, 1.5, 1.5);

  // The query is outside the grid, but point 0 is only 0.3 units away.
  TEST_ASSERT(
      GridTestAccess::nearby_point_ids(grid, query, 0.4).empty(),
      "A query outside the grid should fail and return no points");

  return true;
}


// random test
inline bool test_nearby_point_ids_random() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0);

  std::vector<Point_3> points = rand_points_3(1000, 0.0, 4.0, 0.0, 4.0, 0.0, 4.0);

  std::vector<Vector_3> normals(
      points.size(), Vector_3(0.0, 0.0, 1.0));

  grid.build(points, normals);

  Point_3 query(1.5, 1.5, 1.5);

  TEST_ASSERT(
      nearby_point_ids_contains_all_points_within_radius(
          grid, Point_3(0.05, 0.5, 0.5), 0.5),
      "A query must find nearby points");

  TEST_ASSERT(
      nearby_point_ids_contains_all_points_within_radius(
          grid, Point_3(3.95, 3.5, 3.5), 0.2),
      "A query must find nearby points");

  return true;
}



inline void run_nearby_point_ids_tests(TestData& test_data) {
  test(test_nearby_point_ids_small_radius, "test_nearby_point_ids_small_radius", test_data);
  test(test_nearby_point_ids_neighboring_cells, "test_nearby_point_ids_neighboring_cells", test_data);
  test(test_nearby_point_ids_large_radius, "test_nearby_point_ids_large_radius", test_data);
  test(test_nearby_point_ids_zero_radius, "test_nearby_point_ids_zero_radius", test_data);
  test(test_nearby_point_ids_near_grid_boundary, "test_nearby_point_ids_near_grid_boundary", test_data);
  test(test_nearby_point_ids_query_outside_grid, "test_nearby_point_ids_query_outside_grid", test_data);
  test(test_nearby_point_ids_random, "test_nearby_point_ids_random", test_data);
}