#pragma once
#include "common.h"

inline bool test_compute_block_normals_single_direction() {
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
    Vector_3(0.0, 0.0, 1.0)
  };

  grid.build(points, normals);

  std::vector<Vector_3> block_normals = grid.block_normals();

  TEST_ASSERT(block_normals.size() == 1,
              "There should be one block normal for a single cell");

  TEST_ASSERT(
      CGAL::to_double(block_normals[0].squared_length()) > 0.99,
      "The computed block normal should be a unit vector");

  TEST_ASSERT(
      CGAL::to_double(block_normals[0] * normals[0]) > 0.99,
      "The block normal should align with all associated point normals");

  return true;
}


inline bool test_compute_block_normals_opposing_directions() {
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

  // build() filters out the opposing normal, leaving two positive-z normals.
  std::vector<Vector_3> block_normals = grid.block_normals();

  TEST_ASSERT(block_normals.size() == 1,
              "There should be one block normal for a single cell");

  TEST_ASSERT(
      CGAL::to_double(block_normals[0] * Vector_3(0.0, 0.0, 1.0)) > 0.99,
      "The block normal should align with the filtered positive-z normals");

  return true;
}



// Grid tests

inline void run_compute_block_normals_tests(TestData& test_data) {
  test(test_compute_block_normals_single_direction, "test_compute_block_normals_single_direction", test_data);
  test(test_compute_block_normals_opposing_directions, "test_compute_block_normals_opposing_directions", test_data);
}
