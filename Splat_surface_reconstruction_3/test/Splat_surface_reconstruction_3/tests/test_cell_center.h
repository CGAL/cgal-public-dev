#pragma once
#include "common.h"

inline bool test_cell_center() {
  TestGrid grid = GridTestAccess::make_grid(0.2, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
  int ix = 3, iy = 4, iz = 5;
  Point_3 center = GridTestAccess::cell_center(grid, ix, iy, iz);
  FT expected_x = 0.7;
  FT expected_y = 0.9;
  FT expected_z = 1.1;
  TEST_ASSERT(approxEq(center, Point_3(expected_x, expected_y, expected_z)), "Cell center calculation is incorrect");
  return true;
}

inline bool test_cell_center_with_negative_mins() {
  TestGrid grid = GridTestAccess::make_grid(0.2, -1.0, -2.0, -3.0, 2.0, 2.0, 2.0);
  int ix = 3, iy = 4, iz = 5;
  Point_3 center = GridTestAccess::cell_center(grid, ix, iy, iz);
  FT expected_x = -0.3;
  FT expected_y = -1.1;
  FT expected_z = -1.9;
  TEST_ASSERT(approxEq(center, Point_3(expected_x, expected_y, expected_z)), "Cell center calculation is incorrect");
  return true;
}



inline void run_cell_center_tests(TestData& test_data) {
  test(test_cell_center, "test_cell_center", test_data);
  test(test_cell_center_with_negative_mins, "test_cell_center_with_negative_mins", test_data);
}
