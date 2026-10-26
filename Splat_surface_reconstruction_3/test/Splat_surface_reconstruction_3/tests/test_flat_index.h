#pragma once
#include "common.h"

inline bool test_flat_index_all_equal() {
  TestGrid grid = GridTestAccess::make_grid(0.2, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
  int ix = 1, iy = 2, iz = 3;
  std::size_t index = GridTestAccess::flat_index(grid, ix, iy, iz);
  std::size_t expected_index = 100 + 20 + 3;
  TEST_ASSERT(index == expected_index, "Flat index calculation is incorrect");
  return true;
}

inline bool test_flat_index_different_sizes() {
  TestGrid grid = GridTestAccess::make_grid(0.1, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0);
  int ix = 4, iy = 2, iz = 3;
  std::size_t index = GridTestAccess::flat_index(grid, ix, iy, iz);
  std::size_t expected_index = 2400 + 60 + 3;
  TEST_ASSERT(index == expected_index, "Flat index calculation is incorrect");
  return true;
}


inline void run_flat_index_tests(TestData& test_data) {
  test(test_flat_index_all_equal, "test_flat_index_all_equal", test_data);
  test(test_flat_index_different_sizes, "test_flat_index_different_sizes", test_data);
}
