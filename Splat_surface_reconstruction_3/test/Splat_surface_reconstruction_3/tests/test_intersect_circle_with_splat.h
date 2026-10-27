#pragma once
#include "common.h"

bool is_on_splat(const Point_3& point, const Point_3& splat_center, const Vector_3& splat_normal, FT splat_radius) {
  Vector_3 diff = point - splat_center;
  FT dist_to_plane = CGAL::abs(diff * splat_normal);
  return dist_to_plane <= 1e-6;
}

bool test_intersect_circle_with_splat_simple() {
  // all that matters here is the d (circle is the set of points exactly d away from both pa and pb)
  FT d = 1.0;
  TestGrid grid = GridTestAccess::make_grid(d, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

  Point_3 parent_a(0.0, 0.0, 0.0);
  Point_3 parent_b(1.0, 0.0, 0.0);
  Point_3 splat_center(0.5, 0.0, 0.0);
  Vector_3 splat_normal(0.0, 1.0, 0.0);
  FT splat_radius = 1.0;

  std::size_t num_expected_points = 2;
  
  std::vector<Point_3> actual_points = GridTestAccess::intersect_circle_with_splat(
    grid, parent_a, parent_b, splat_center, splat_normal, splat_radius);

  TEST_ASSERT(actual_points.size() == num_expected_points, "Unexpected number of intersection points");

  for (const Point_3& p : actual_points) {
    TEST_ASSERT(is_on_splat(p, splat_center, splat_normal, splat_radius), "Point is not on splat");

    FT distance_from_p_a = CGAL::sqrt(CGAL::squared_distance(p, parent_a));
    FT distance_from_p_b = CGAL::sqrt(CGAL::squared_distance(p, parent_b));
    
    TEST_ASSERT(CGAL::abs(distance_from_p_a - d) <= 1e-6, "Point is not distance d from parent points");
    TEST_ASSERT(CGAL::abs(distance_from_p_b - d) <= 1e-6, "Point is not distance d from parent points");
  }

  // They must be different points
  if (actual_points.size() == 2) {
    TEST_ASSERT(CGAL::squared_distance(actual_points[0], actual_points[1]) > 1e-6, "Intersection points are not different");
  }

  return true;
}

bool test_intersect_circle_with_splat_other_normal() {
  // all that matters here is the d (circle is the set of points exactly d away from both pa and pb)
  FT d = 1.0;
  TestGrid grid = GridTestAccess::make_grid(d, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

  Point_3 parent_a(0.0, 0.0, 0.0);
  Point_3 parent_b(1.0, 0.0, 0.0);
  Point_3 splat_center(0.5, 0.0, 0.0);
  Vector_3 splat_normal(0.0, 0.0, 1.0);
  FT splat_radius = 1.0;

  std::size_t num_expected_points = 2;
  
  std::vector<Point_3> actual_points = GridTestAccess::intersect_circle_with_splat(
    grid, parent_a, parent_b, splat_center, splat_normal, splat_radius);

  TEST_ASSERT(actual_points.size() == num_expected_points, "Unexpected number of intersection points");

  for (const Point_3& p : actual_points) {
    TEST_ASSERT(is_on_splat(p, splat_center, splat_normal, splat_radius), "Point is not on splat");

    FT distance_from_p_a = CGAL::sqrt(CGAL::squared_distance(p, parent_a));
    FT distance_from_p_b = CGAL::sqrt(CGAL::squared_distance(p, parent_b));
    
    TEST_ASSERT(CGAL::abs(distance_from_p_a - d) <= 1e-6, "Point is not distance d from parent points");
    TEST_ASSERT(CGAL::abs(distance_from_p_b - d) <= 1e-6, "Point is not distance d from parent points");
  }

  // They must be different points
  if (actual_points.size() == 2) {
    TEST_ASSERT(CGAL::squared_distance(actual_points[0], actual_points[1]) > 1e-6, "Intersection points are not different");
  }

  return true;
}


bool test_intersect_circle_with_splat_no_intersection() {
  // all that matters here is the d (circle is the set of points exactly d away from both pa and pb)
  FT d = 1.0;
  TestGrid grid = GridTestAccess::make_grid(d, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

  Point_3 parent_a(0.0, 0.0, 0.0);
  Point_3 parent_b(1.0, 0.0, 0.0);
  Point_3 splat_center(0.5, 0.0, 0.0);
  Vector_3 splat_normal(0.0, 1.0, 0.0);

  // Barely not touching either side
  FT splat_radius = 0.86603 - 0.01;

  std::size_t num_expected_points = 0;
  
  std::vector<Point_3> actual_points = GridTestAccess::intersect_circle_with_splat(
    grid, parent_a, parent_b, splat_center, splat_normal, splat_radius);

  TEST_ASSERT(actual_points.size() == num_expected_points, "Unexpected number of intersection points");

  return true;
}

bool test_intersect_circle_with_splat_single_intersection() {
  // all that matters here is the d (circle is the set of points exactly d away from both pa and pb)
  FT d = 1.0;
  TestGrid grid = GridTestAccess::make_grid(d, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

  Point_3 parent_a(0.0, 0.0, 0.0);
  Point_3 parent_b(1.0, 0.0, 0.0);
  Point_3 splat_center(0.5, 0.0, 0.86);
  Vector_3 splat_normal(0.0, 1.0, 0.0);
  FT splat_radius = 0.1;

  std::size_t num_expected_points = 1;
  
  std::vector<Point_3> actual_points = GridTestAccess::intersect_circle_with_splat(
    grid, parent_a, parent_b, splat_center, splat_normal, splat_radius);

  TEST_ASSERT(actual_points.size() == num_expected_points, "Unexpected number of intersection points");

  for (const Point_3& p : actual_points) {
    TEST_ASSERT(is_on_splat(p, splat_center, splat_normal, splat_radius), "Point is not on splat");

    FT distance_from_p_a = CGAL::sqrt(CGAL::squared_distance(p, parent_a));
    FT distance_from_p_b = CGAL::sqrt(CGAL::squared_distance(p, parent_b));
    
    TEST_ASSERT(CGAL::abs(distance_from_p_a - d) <= 1e-6, "Point is not distance d from parent points");
    TEST_ASSERT(CGAL::abs(distance_from_p_b - d) <= 1e-6, "Point is not distance d from parent points");
  }

  return true;
}


bool test_intersect_circle_with_splat_single_equidistant_result() {
  // all that matters here is the d (circle is the set of points exactly d away from both pa and pb)
  FT d = 0.5 + 1e-6;
  TestGrid grid = GridTestAccess::make_grid(d, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

  Point_3 parent_a(0.0, 0.0, 0.0);
  Point_3 parent_b(1.0, 0.0, 0.0);
  Point_3 splat_center(0.5, 0.0, 0.0);
  Vector_3 splat_normal(0.0, 1.0, 0.0);
  FT splat_radius = 1.0;

  // It's still counted as two distinct intersection points due to the small perturbation in d
  std::size_t num_expected_points = 2;
  
  std::vector<Point_3> actual_points = GridTestAccess::intersect_circle_with_splat(
    grid, parent_a, parent_b, splat_center, splat_normal, splat_radius);

  TEST_ASSERT(actual_points.size() == num_expected_points, "Unexpected number of intersection points");

  for (const Point_3& p : actual_points) {
    TEST_ASSERT(is_on_splat(p, splat_center, splat_normal, splat_radius), "Point is not on splat");

    FT distance_from_p_a = CGAL::sqrt(CGAL::squared_distance(p, parent_a));
    FT distance_from_p_b = CGAL::sqrt(CGAL::squared_distance(p, parent_b));
    
    TEST_ASSERT(CGAL::abs(distance_from_p_a - d) <= 1e-6, "Point is not distance d from parent points");
    TEST_ASSERT(CGAL::abs(distance_from_p_b - d) <= 1e-6, "Point is not distance d from parent points");
  }

  return true;
}
bool test_intersect_circle_with_splat_single_splat_in_same_plane() {
  // all that matters here is the d (circle is the set of points exactly d away from both pa and pb)
  FT d = 0.5 + 1e-6;
  TestGrid grid = GridTestAccess::make_grid(d, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

  Point_3 parent_a(0.0, 0.0, 0.0);
  Point_3 parent_b(1.0, 0.0, 0.0);
  Point_3 splat_center(0.5, 0.0, 0.0);
  Vector_3 splat_normal(1.0, 0.0, 0.0);
  FT splat_radius = 1.0;

  std::size_t num_expected_points = 0;
  
  std::vector<Point_3> actual_points = GridTestAccess::intersect_circle_with_splat(
    grid, parent_a, parent_b, splat_center, splat_normal, splat_radius);

  TEST_ASSERT(actual_points.size() == num_expected_points, "Unexpected number of intersection points");

  return true;
}

void run_intersect_circle_with_splat_tests(TestData& test_data) {
  test(test_intersect_circle_with_splat_simple, "test_intersect_circle_with_splat_simple", test_data);
  test(test_intersect_circle_with_splat_other_normal, "test_intersect_circle_with_splat_other_normal", test_data);
  test(test_intersect_circle_with_splat_no_intersection, "test_intersect_circle_with_splat_no_intersection", test_data);
  test(test_intersect_circle_with_splat_single_intersection, "test_intersect_circle_with_splat_single_intersection", test_data);
  test(test_intersect_circle_with_splat_single_equidistant_result, "test_intersect_circle_with_splat_single_equidistant_result", test_data);
  test(test_intersect_circle_with_splat_single_splat_in_same_plane, "test_intersect_circle_with_splat_single_splat_in_same_plane", test_data);
}