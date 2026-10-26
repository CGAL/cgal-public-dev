#pragma once
#include "common.h"



inline bool test_build_cell_point_associations() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0);

  std::vector<Point_3> points = {
    Point_3(1.5, 1.5, 1.5), // 0: cell (1,1,1)
    Point_3(2.5, 1.5, 1.5), // 1: cell (2,1,1)
    Point_3(3.5, 1.5, 1.5), // 2: cell (3,1,1)
    Point_3(2.0, 1.5, 1.5), // 3: on a cell boundary
    Point_3(0.5, 0.5, 0.5)  // 4: cell (0,0,0)
  };

  std::vector<Vector_3> normals(
      points.size(), Vector_3(0.0, 0.0, 1.0));

  grid.build(points, normals);

  // Containment: each point belongs to exactly one cell.
  TEST_ASSERT(
      GridTestAccess::point_ids_in_box(grid, 1, 1, 1)
          == std::vector<std::size_t>{0},
      "Cell (1,1,1) should contain point 0");

  TEST_ASSERT(
      GridTestAccess::point_ids_in_box(grid, 2, 1, 1)
          == (std::vector<std::size_t>{1, 3}),
      "Cell (2,1,1) should contain points 1 and 3");

  TEST_ASSERT(
      GridTestAccess::point_ids_in_box(grid, 3, 1, 1)
          == std::vector<std::size_t>{2},
      "Cell (3,1,1) should contain point 2");

  TEST_ASSERT(
      GridTestAccess::point_ids_in_box(grid, 0, 0, 0)
          == std::vector<std::size_t>{4},
      "Cell (0,0,0) should contain point 4");

  // Point 0 is within d of face-, edge-, and corner-adjacent boxes.
  TEST_ASSERT(
      GridTestAccess::point_ids_within_distance_d_of_box(grid, 2, 1, 1)
          .size() > 0,
      "Face-adjacent box should associate nearby points");

  TEST_ASSERT(
      GridTestAccess::point_ids_within_distance_d_of_box(grid, 2, 2, 1)
          .size() > 0,
      "Edge-adjacent box should associate nearby points");

  TEST_ASSERT(
      GridTestAccess::point_ids_within_distance_d_of_box(grid, 2, 2, 2)
          .size() > 0,
      "Corner-adjacent box should associate nearby points");

  // Cell (3,1,1) is 1.5 units from point 0, but may still associate
  // with point 1 or point 2. Check membership of point 0 specifically.
  const auto& far_cell =
      GridTestAccess::point_ids_within_distance_d_of_box(grid, 3, 1, 1);

  TEST_ASSERT(
      std::find(far_cell.begin(), far_cell.end(), 0) == far_cell.end(),
      "Cell (3,1,1) should not associate point 0");

  // Cell (3,1,1) is the containing cell of point 2.
  TEST_ASSERT(
      std::find(far_cell.begin(), far_cell.end(), 2) != far_cell.end(),
      "A point should be associated with its own cell");

  return true;
}


inline bool test_build_cell_point_associations_distance_check() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0);

  std::vector<Point_3> points = {
    Point_3(1.01, 1.5, 1.5), // 0: near left edge of cell (1,1,1)
    Point_3(1.99, 1.5, 1.5), // 1: near right edge of cell (1,1,1)
  };

  std::vector<Vector_3> normals(
      points.size(), Vector_3(0.0, 0.0, 1.0));

  grid.build(points, normals);

  // Both points are contained in cell (1,1,1).
  TEST_ASSERT(
      GridTestAccess::point_ids_in_box(grid, 1, 1, 1)
          == (std::vector<std::size_t>{0, 1}),
      "Cell (1,1,1) should contain both points");

  const auto& left_neighbor =
      GridTestAccess::point_ids_within_distance_d_of_box(grid, 0, 1, 1);
  const auto& right_neighbor =
      GridTestAccess::point_ids_within_distance_d_of_box(grid, 2, 1, 1);

  // Point 0 is 0.01 from the left face, so it is within d of the left neighbor.
  TEST_ASSERT(
      std::find(left_neighbor.begin(), left_neighbor.end(), 0)
          != left_neighbor.end(),
      "The left neighbor should associate point 0");

  // Point 1 is 0.01 from the right face, so it is within d of the right neighbor.
  TEST_ASSERT(
      std::find(right_neighbor.begin(), right_neighbor.end(), 1)
          != right_neighbor.end(),
      "The right neighbor should associate point 1");

  // With d = 1, these points are not more than d away from their adjacent
  // boxes. To test a point more than d away, use a larger offset in a
  // separate scenario or test a cell two steps away.

  const auto& two_cells_away =
      GridTestAccess::point_ids_within_distance_d_of_box(grid, 3, 1, 1);

  TEST_ASSERT(
      std::find(two_cells_away.begin(), two_cells_away.end(), 0)
          == two_cells_away.end(),
      "A box beyond distance d should not associate point 0");

  return true;
}





inline bool test_build_point_associations_multiple_points() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 3.0, 1.0, 1.0);

  std::vector<Point_3> points = {
    Point_3(0.5, 0.5, 0.5),
    Point_3(1.5, 0.5, 0.5),
    Point_3(2.5, 0.5, 0.5)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  grid.build(points, normals);



  // Each point is contained in exactly one cell.
  TEST_ASSERT(
      (GridTestAccess::point_ids_in_box(grid, 0, 0, 0)
          == std::vector<std::size_t>{0}),
      "The first cell should contain only point 0");

  TEST_ASSERT(
      (GridTestAccess::point_ids_in_box(grid, 1, 0, 0)
          == std::vector<std::size_t>{1}),
      "The second cell should contain only point 1");

  TEST_ASSERT(
      (GridTestAccess::point_ids_in_box(grid, 2, 0, 0)
          == std::vector<std::size_t>{2}),
      "The third cell should contain only point 2");

  // Each cell should also associate nearby points from adjacent cells.
  TEST_ASSERT(
      (GridTestAccess::point_ids_within_distance_d_of_box(grid, 0, 0, 0)
          == std::vector<std::size_t>{0, 1}),
      "The first cell should associate points 0 and 1");

  TEST_ASSERT(
      (GridTestAccess::point_ids_within_distance_d_of_box(grid, 1, 0, 0)
          == std::vector<std::size_t>{0, 1, 2}),
      "The middle cell should associate all three points");

  TEST_ASSERT(
      (GridTestAccess::point_ids_within_distance_d_of_box(grid, 2, 0, 0)
          == std::vector<std::size_t>{1, 2}),
      "The last cell should associate points 1 and 2");

  return true;
}


inline bool test_build_normal_filtering_removes_opposing_normal() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);

  std::vector<Point_3> points = {
    Point_3(0.5, 0.5, 0.5),
    Point_3(0.5, 0.5, 0.5),
    Point_3(0.5, 0.5, 0.5)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, -1.0)
  };

  grid.build(points, normals);

  // The unnormalized normal sum is (0,0,1), with length 1.
  // The negative-z normal should therefore be filtered out.
  TEST_ASSERT(
      (GridTestAccess::point_ids_within_distance_d_of_box(grid, 0, 0, 0)
          == std::vector<std::size_t>{0, 1}),
      "Filtering should remove the point with an opposing normal");

  // Filtering should not alter cell containment.
  TEST_ASSERT(
      (GridTestAccess::point_ids_in_box(grid, 0, 0, 0)
          == std::vector<std::size_t>{0, 1, 2}),
      "Filtering should not change the contained-point list");

  return true;
}


inline bool test_build_normal_filtering_keeps_points_when_sum_is_small() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);

  std::vector<Point_3> points = {
    Point_3(0.5, 0.5, 0.5),
    Point_3(0.5, 0.5, 0.5)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, -1.0)
  };

  grid.build(points, normals);

  // The sum is zero, below the threshold of 0.1.
  // The filtering rule must preserve all points in this case.
  TEST_ASSERT(
      (GridTestAccess::point_ids_within_distance_d_of_box(grid, 0, 0, 0)
          == std::vector<std::size_t>{0, 1}),
      "Filtering should retain all points when the normal sum is too small");

  return true;
}



inline void run_build_cell_point_associations_tests(TestData& test_data) {
  test(test_build_cell_point_associations, "test_build_cell_point_associations", test_data);
  test(test_build_cell_point_associations_distance_check, "test_build_cell_point_associations_distance_check", test_data);
  test(test_build_point_associations_multiple_points, "test_build_point_associations_multiple_points", test_data);
  test(test_build_normal_filtering_removes_opposing_normal, "test_build_normal_filtering_removes_opposing_normal", test_data);
  test(test_build_normal_filtering_keeps_points_when_sum_is_small, "test_build_normal_filtering_keeps_points_when_sum_is_small", test_data);
}
