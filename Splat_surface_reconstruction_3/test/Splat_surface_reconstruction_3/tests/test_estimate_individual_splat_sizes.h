#pragma once
#include "common.h"

inline TestGrid make_test_grid(std::vector<Point> points, std::vector<Vector_3> normals) {
  TestGrid grid = GridTestAccess::make_grid(2.0, -2.0, -2.0, -2.0, 2.0, 2.0, 2.0);
  grid.build(points, normals);
  return grid;
}

inline bool test_estimate_individual_splat_sizes_three_splats() {
  FT sqrt_3 = std::sqrt(3.0);
  // equilateral triangle
  std::vector<Point> points = {
    Point(-1.0, 0.0, 0.0),
    Point(0.0, sqrt_3, 0.0),
    Point(1.0, 0.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };
  TestGrid grid = make_test_grid(points, normals);

  std::vector<FT> splat_sizes = GridTestAccess::estimate_individual_splat_sizes(grid);

  TEST_ASSERT(splat_sizes.size() == points.size(), "Splat sizes vector has incorrect size");

  std::cout << "Splat sizes: ";
  for (const auto& size : splat_sizes) {
    std::cout << size << " ";
  }
  std::cout << std::endl;

  std::cout << "Expected splat size: " << 2.0 / sqrt_3 << std::endl;

  TEST_ASSERT(approxEq(splat_sizes[0], 2.0 / sqrt_3), "First splat size incorrect");
  TEST_ASSERT(approxEq(splat_sizes[1], 2.0 / sqrt_3), "Second splat size incorrect");
  TEST_ASSERT(approxEq(splat_sizes[2], 2.0 / sqrt_3), "Third splat size incorrect");

  return true;
}

inline bool test_estimate_individual_splat_sizes_three_splats_uneven_triangle() {
  FT sqrt_3 = std::sqrt(3.0);
  // equilateral triangle
  std::vector<Point> points = {
    Point(-1.0, 0.0, 0.0),
    Point(0.0, 1.0, 0.0),
    Point(1.0, 0.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };
  TestGrid grid = make_test_grid(points, normals);

  std::vector<FT> splat_sizes = GridTestAccess::estimate_individual_splat_sizes(grid);

  TEST_ASSERT(splat_sizes.size() == points.size(), "Splat sizes vector has incorrect size");

  std::cout << "Splat sizes: ";
  for (const auto& size : splat_sizes) {
    std::cout << size << " ";
  }
  std::cout << std::endl;

  std::cout << "Expected splat size: " << 1.0 << std::endl;

  TEST_ASSERT(approxEq(splat_sizes[0], 1.0), "First splat size incorrect");
  TEST_ASSERT(approxEq(splat_sizes[1], 1.0), "Second splat size incorrect");
  TEST_ASSERT(approxEq(splat_sizes[2], 1.0), "Third splat size incorrect");

  return true;
}

inline bool test_estimate_individual_splat_sizes_uneven_spacing() {
  std::vector<Point> points = {
    Point(-1.0, 0.0, 0.0),
    Point(-0.8, 0.0, 0.0),
    Point( 0.0, 0.0, 0.0),
    Point( 1.0, 0.0, 0.0)
  };

  std::vector<Vector_3> normals(points.size(),
                                Vector_3(0.0, 0.0, 1.0));

  TestGrid grid = make_test_grid(points, normals);
  std::vector<FT> sizes =
      GridTestAccess::estimate_individual_splat_sizes(grid);

  TEST_ASSERT(sizes.size() == points.size(),
              "Incorrect number of splat sizes");

  for (const FT size : sizes) {
    TEST_ASSERT(std::isfinite(CGAL::to_double(size)) &&
                size >= FT(0),
                "Splat size must be finite and non-negative");
  }

  return true;
}


inline void run_estimate_individual_splat_sizes_tests(TestData& test_data) {
  test(test_estimate_individual_splat_sizes_three_splats, "test_estimate_individual_splat_sizes_three_splats", test_data);
  test(test_estimate_individual_splat_sizes_three_splats_uneven_triangle, "test_estimate_individual_splat_sizes_three_splats_uneven_triangle", test_data);
  test(test_estimate_individual_splat_sizes_uneven_spacing, "test_estimate_individual_splat_sizes_uneven_spacing", test_data);
}
