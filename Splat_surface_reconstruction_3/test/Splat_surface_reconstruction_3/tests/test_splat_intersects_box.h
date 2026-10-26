#pragma once
#include "common.h"

// Tests for splat_intersects_box
inline void test_splat_intersects_box_basic_case() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0);

  Point_3 point(1.0, 0.5, 0.5);
  Vector_3 normal(0.0, 0.0, 1.0);
  FT size = 0.5;

  grid.build(std::vector<Point_3>{point},
             std::vector<Vector_3>{normal});
  grid.set_all_splat_sizes(size);

  TEST_ASSERT(
      GridTestAccess::splat_intersects_box(grid, 0, 0, 0, 0),
      "Disk centered on the box boundary should intersect the box");
}

inline void test_splat_intersects_box_center_inside() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0);

  Point_3 point(0.5, 0.5, 0.5);
  Vector_3 normal(0.0, 0.0, 1.0);

  grid.build(std::vector<Point_3>{point},
             std::vector<Vector_3>{normal});
  grid.set_all_splat_sizes(0.1);

  TEST_ASSERT(
      GridTestAccess::splat_intersects_box(grid, 0, 0, 0, 0),
      "Disk centered inside the box should intersect it");
}

inline void test_splat_intersects_box_tangent_to_face() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0);

  Point_3 point(1.5, 0.5, 0.5);
  Vector_3 normal(0.0, 0.0, 1.0);

  grid.build(std::vector<Point_3>{point},
             std::vector<Vector_3>{normal});
  grid.set_all_splat_sizes(0.5);

  TEST_ASSERT(
      GridTestAccess::splat_intersects_box(grid, 0, 0, 0, 0),
      "Disk tangent to the box face should count as intersecting");
}

inline void test_splat_does_not_intersect_box_outside_radius() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0);

  Point_3 point(1.6, 0.5, 0.5);
  Vector_3 normal(0.0, 0.0, 1.0);

  grid.build(std::vector<Point_3>{point},
             std::vector<Vector_3>{normal});
  grid.set_all_splat_sizes(0.5);

  TEST_ASSERT(
      !GridTestAccess::splat_intersects_box(grid, 0, 0, 0, 0),
      "Disk separated from the box by more than its radius should not intersect");
}

inline void test_splat_does_not_intersect_box_wrong_plane() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0);

  // The disk lies in the plane z = 1.1, entirely above the box.
  Point_3 point(0.5, 0.5, 1.1);
  Vector_3 normal(0.0, 0.0, 1.0);

  grid.build(std::vector<Point_3>{point},
             std::vector<Vector_3>{normal});
  grid.set_all_splat_sizes(0.05);

  TEST_ASSERT(
      !GridTestAccess::splat_intersects_box(grid, 0, 0, 0, 0),
      "Disk outside the box's plane and radius should not intersect");
}

inline void test_splat_disk_does_not_intersect_box_near_sphere_only() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0);

  // Disk lies in z = 0.5, with radius 0.2.
  // The neighboring cell is separated in z, so the disk cannot reach it.
  Point_3 point(0.5, 0.5, 0.5);
  Vector_3 normal(0.0, 0.0, 1.0);

  grid.build(std::vector<Point_3>{point},
             std::vector<Vector_3>{normal});
  grid.set_all_splat_sizes(0.2);

  // Cell (0, 0, 1) occupies [0,1] x [0,1] x [1,2].
  // The disk at z=0.5 cannot intersect this cell.
  TEST_ASSERT(
      !GridTestAccess::splat_intersects_box(grid, 0, 0, 0, 1),
      "A disk must not intersect a cell outside its plane and radius");
}

inline void test_splat_intersects_box_with_tilted_normal() {
  TestGrid grid = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0);

  Point_3 point(1.0, 0.5, 0.5);
  Vector_3 normal(1.0, 0.0, 0.0);

  grid.build(std::vector<Point_3>{point},
             std::vector<Vector_3>{normal});
  grid.set_all_splat_sizes(0.5);

  // Disk lies in x = 1 and extends in y and z.
  TEST_ASSERT(
      GridTestAccess::splat_intersects_box(grid, 0, 0, 0, 0),
      "Disk in a vertical plane should intersect the box");
}

inline void test_splat_intersection_depends_on_normal() {
  TestGrid grid_z = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0);

  TestGrid grid_x = GridTestAccess::make_grid(
      1.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0);

  Point_3 point(0.5, 0.5, 0.5);

  grid_z.build(std::vector<Point_3>{point},
               std::vector<Vector_3>{Vector_3(0.0, 0.0, 1.0)});
  grid_x.build(std::vector<Point_3>{point},
               std::vector<Vector_3>{Vector_3(1.0, 0.0, 0.0)});

  grid_z.set_all_splat_sizes(0.5);
  grid_x.set_all_splat_sizes(0.5);

  // Both disks intersect the cell containing their common center.
  TEST_ASSERT(
      GridTestAccess::splat_intersects_box(grid_z, 0, 0, 0, 0),
      "Horizontal disk should intersect the center cell");
  TEST_ASSERT(
      GridTestAccess::splat_intersects_box(grid_x, 0, 0, 0, 0),
      "Vertical disk should intersect the center cell");
}




inline void run_splat_intersects_box_tests(TestData& test_data) {
  test(test_splat_intersects_box_basic_case, "test_splat_intersects_box_basic_case", test_data);
  test(test_splat_intersects_box_center_inside, "test_splat_intersects_box_center_inside", test_data);
  test(test_splat_intersects_box_tangent_to_face, "test_splat_intersects_box_tangent_to_face", test_data);
  test(test_splat_does_not_intersect_box_outside_radius, "test_splat_does_not_intersect_box_outside_radius", test_data);
  test(test_splat_does_not_intersect_box_wrong_plane, "test_splat_does_not_intersect_box_wrong_plane", test_data);
  test(test_splat_disk_does_not_intersect_box_near_sphere_only, "test_splat_disk_does_not_intersect_box_near_sphere_only", test_data);
  test(test_splat_intersects_box_with_tilted_normal, "test_splat_intersects_box_with_tilted_normal", test_data);
  test(test_splat_intersection_depends_on_normal, "test_splat_intersection_depends_on_normal", test_data);
}
