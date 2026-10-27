#pragma once
#include "common.h"

inline bool test_grid(TestGrid& grid, std::size_t expected_nx, std::size_t expected_ny, std::size_t expected_nz) {
  TEST_ASSERT(GridTestAccess::nx(grid) == expected_nx, "nx is not initialized correctly");
  TEST_ASSERT(GridTestAccess::ny(grid) == expected_ny, "ny is not initialized correctly");
  TEST_ASSERT(GridTestAccess::nz(grid) == expected_nz, "nz is not initialized correctly");

  const std::vector<Cell>& cells = GridTestAccess::cells(grid);
  TEST_ASSERT(cells.size() == expected_nx * expected_ny * expected_nz, "Cells are not initialized correctly");

  return true;
}

inline bool test_initialize_grid() {
  // Initialize grid is called with the constructor
  TestGrid grid = GridTestAccess::make_grid(0.2, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

  return test_grid(grid, 10, 10, 10);
}

inline bool test_initialize_grid_different_sizes() {
  // Initialize grid is called with the constructor
  TestGrid grid = GridTestAccess::make_grid(0.2, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0);

  return test_grid(grid, 5, 10, 15);
}

inline bool test_nitialize_grid_odd_box_size() {
  // Initialize grid is called with the constructor
  TestGrid grid = GridTestAccess::make_grid(0.8, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

  return test_grid(grid, 3, 3, 3);
}

inline bool test_nitialize_grid_minimum_size() {
  // Initialize grid is called with the constructor
  TestGrid grid = GridTestAccess::make_grid(4.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

  return test_grid(grid, 1, 1, 1);
}

inline bool test_nitialize_grid_zero_size() {
  // Initialize grid is called with the constructor
  TestGrid grid = GridTestAccess::make_grid(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  return test_grid(grid, 1, 1, 1);
}



inline void run_initialize_grid_tests(TestData& test_data) {
  test(test_initialize_grid, "test_initialize_grid", test_data);
  test(test_initialize_grid_different_sizes, "test_initialize_grid_different_sizes", test_data);
  test(test_nitialize_grid_odd_box_size, "test_nitialize_grid_odd_box_size", test_data);
  test(test_nitialize_grid_minimum_size, "test_nitialize_grid_minimum_size", test_data);
  test(test_nitialize_grid_zero_size, "test_nitialize_grid_zero_size", test_data);
}
