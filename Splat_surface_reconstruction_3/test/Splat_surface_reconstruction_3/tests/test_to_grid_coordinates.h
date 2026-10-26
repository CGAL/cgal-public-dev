#pragma once
#include "common.h"

inline bool test_to_grid_coordinates_lowest_corner() {
  TestGrid grid = GridTestAccess::make_grid(1.0, -2.0, -2.0, -2.0, 2.0, 2.0, 2.0);
  int ix, iy, iz;
  Point_3 p(-2.0, -2.0, -2.0);
  bool result = GridTestAccess::to_grid_coordinates(grid, p, ix, iy, iz);
  TEST_ASSERT(result, "Failed to convert lowest corner to grid coordinates");
  TEST_ASSERT(ix == 0 && iy == 0 && iz == 0, "Lowest corner grid coordinates are incorrect");
  return true;
}

inline bool test_to_grid_coordinates_highest_corner() {
  TestGrid grid = GridTestAccess::make_grid(1.0, -2.0, -2.0, -2.0, 2.0, 2.0, 2.0);
  int ix, iy, iz;
  Point_3 p(2.0, 2.0, 2.0);
  bool result = GridTestAccess::to_grid_coordinates(grid, p, ix, iy, iz);
  TEST_ASSERT(result, "Failed to convert highest corner to grid coordinates");
  TEST_ASSERT(ix == 3 && iy == 3 && iz == 3, "Highest corner grid coordinates are incorrect");
  return true;
}

inline bool test_to_grid_coordinates_box_corner() {
  TestGrid grid = GridTestAccess::make_grid(1.0, -2.0, -2.0, -2.0, 2.0, 2.0, 2.0);
  int ix, iy, iz;
  Point_3 p(1.0, 1.0, 1.0);
  bool result = GridTestAccess::to_grid_coordinates(grid, p, ix, iy, iz);
  TEST_ASSERT(result, "Failed to convert box corner to grid coordinates");
  TEST_ASSERT(ix == 3 && iy == 3 && iz == 3, "Box corner grid coordinates are incorrect");
  return true;
}

inline bool test_to_grid_coordinates_all_different() {
  TestGrid grid = GridTestAccess::make_grid(1.0, -2.0, -2.0, -2.0, 2.0, 2.0, 2.0);
  int ix, iy, iz;
  Point_3 p(-1.5, 0.5, 1.5);
  bool result = GridTestAccess::to_grid_coordinates(grid, p, ix, iy, iz);
  TEST_ASSERT(result, "Failed to convert all different point to grid coordinates");
  TEST_ASSERT(ix == 0 && iy == 2 && iz == 3, "All different point grid coordinates are incorrect");
  return true;
}

inline bool test_to_grid_coordinates_out_of_bounds_1() {
  TestGrid grid = GridTestAccess::make_grid(1.0, -2.0, -2.0, -2.0, 2.0, 2.0, 2.0);
  int ix, iy, iz;
  Point_3 p(1.0, 1.0, 2.001);
  bool result = GridTestAccess::to_grid_coordinates(grid, p, ix, iy, iz);
  TEST_ASSERT(!result, "Out of bounds point should not convert to grid coordinates");
  return true;
}

inline bool test_to_grid_coordinates_out_of_bounds_2() {
  TestGrid grid = GridTestAccess::make_grid(1.0, -2.0, -2.0, -2.0, 2.0, 2.0, 2.0);
  int ix, iy, iz;
  Point_3 p(1.0, 1.0, -2.001);
  bool result = GridTestAccess::to_grid_coordinates(grid, p, ix, iy, iz);
  TEST_ASSERT(!result, "Out of bounds point should not convert to grid coordinates");
  return true;
}

inline bool test_to_grid_coordinates_random() {
  TestGrid grid = GridTestAccess::make_grid(1.0, -2.0, -2.0, -2.0, 2.0, 2.0, 2.0);
  int ix, iy, iz;
  int num_iterations = 1000;
  for (int i = 0; i < num_iterations; ++i) {
    int chosen_x = rand_int(0, 3);
    int chosen_y = rand_int(0, 3);
    int chosen_z = rand_int(0, 3);

    FT min_x = -2.0 + chosen_x * 1.0;
    FT min_y = -2.0 + chosen_y * 1.0;
    FT min_z = -2.0 + chosen_z * 1.0;

    FT x = rand_float_exclusive(min_x, min_x + 1.0);
    FT y = rand_float_exclusive(min_y, min_y + 1.0);
    FT z = rand_float_exclusive(min_z, min_z + 1.0);

    Point_3 p(x, y, z);
    bool result = GridTestAccess::to_grid_coordinates(grid, p, ix, iy, iz);
    TEST_ASSERT(result, "Failed to convert all different point to grid coordinates");
    TEST_ASSERT(ix == chosen_x && iy == chosen_y && iz == chosen_z, "All different point grid coordinates are incorrect");
  }
  return true;
}

inline bool test_to_grid_coordinates_random_different_box_size() {
  TestGrid grid = GridTestAccess::make_grid(0.2, -2.0, -4.0, -2.0, 2.0, 2.0, 6.0);
  int ix, iy, iz;
  int num_iterations = 1000;
  for (int i = 0; i < num_iterations; ++i) {
    int chosen_x = rand_int(0, 19);
    int chosen_y = rand_int(0, 29);
    int chosen_z = rand_int(0, 39);

    FT min_x = -2.0 + chosen_x * 0.2;
    FT min_y = -4.0 + chosen_y * 0.2;
    FT min_z = -2.0 + chosen_z * 0.2;

    FT x = rand_float_exclusive(min_x, min_x + 0.2);
    FT y = rand_float_exclusive(min_y, min_y + 0.2);
    FT z = rand_float_exclusive(min_z, min_z + 0.2);

    Point_3 p(x, y, z);
    bool result = GridTestAccess::to_grid_coordinates(grid, p, ix, iy, iz);
    TEST_ASSERT(result, "Failed to convert all different point to grid coordinates");
    TEST_ASSERT(ix == chosen_x && iy == chosen_y && iz == chosen_z, "All different point grid coordinates are incorrect");
  }
  return true;
}



inline void run_to_grid_coordinates_tests(TestData& test_data) {
  test(test_to_grid_coordinates_lowest_corner, "test_to_grid_coordinates_lowest_corner", test_data);
  test(test_to_grid_coordinates_highest_corner, "test_to_grid_coordinates_highest_corner", test_data);
  test(test_to_grid_coordinates_box_corner, "test_to_grid_coordinates_box_corner", test_data);
  test(test_to_grid_coordinates_all_different, "test_to_grid_coordinates_all_different", test_data);
  test(test_to_grid_coordinates_out_of_bounds_1, "test_to_grid_coordinates_out_of_bounds_1", test_data);
  test(test_to_grid_coordinates_out_of_bounds_2, "test_to_grid_coordinates_out_of_bounds_2", test_data);
  test(test_to_grid_coordinates_random, "test_to_grid_coordinates_random", test_data);
  test(test_to_grid_coordinates_random_different_box_size, "test_to_grid_coordinates_random_different_box_size", test_data);
}
