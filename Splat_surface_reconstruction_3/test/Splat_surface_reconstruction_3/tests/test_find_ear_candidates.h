#pragma once
#include "common.h"

inline bool test_find_ear_candidates_hexagon() {
  std::vector<Point> points = {
    Point(1.0, 0.0, 0.0),
    Point(0.5, std::sqrt(3)/2, 0.0),
    Point(-0.5, std::sqrt(3)/2, 0.0),
    Point(-1.0, 0.0, 0.0),
    Point(-0.5, -std::sqrt(3)/2, 0.0),
    Point(0.5, -std::sqrt(3)/2, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  const auto candidates = find_ear_candidates_helper(points, normals);

  std::cout << "Found " << candidates.size() << " candidates:\n";
  for (const auto& c : candidates) {
    std::cout << "(" << c.first << "," << c.second << ")\n";
  }

  const double angleRad = 120 * (CGAL_PI / 180.0);

  const std::vector<std::pair<double, std::size_t>> expected_candidates = {
    std::pair<double, std::size_t> { angleRad, 0 },
    std::pair<double, std::size_t> { angleRad, 1 },
    std::pair<double, std::size_t> { angleRad, 2 },
    std::pair<double, std::size_t> { angleRad, 3 },
    std::pair<double, std::size_t> { angleRad, 4 },
    std::pair<double, std::size_t> { angleRad, 5 },
  };



  TEST_ASSERT(candidates.size() == expected_candidates.size(), "Number of candidates does not match expected");
  for (std::size_t i = 0; i < candidates.size(); ++i) {
    bool found = false;
    for (std::size_t j = 0; j < expected_candidates.size(); ++j) {
      if (approxEq(candidates[i].first, expected_candidates[j].first) &&
          approxEq(candidates[i].second, expected_candidates[j].second)) {
        found = true;
        break;
      }
    }

    std::cout << "Checking candidate (" << candidates[i].first << "," << candidates[i].second << "): "
              << (found ? "found" : "not found") << "\n";
    TEST_ASSERT(found, "Candidate not found in expected candidates");
  }

  return true;
}

inline bool test_find_ear_candidates_step()
{
  using Access = TestAccess;
  using halfedge_descriptor = Access::halfedge_descriptor;
  using vertex_descriptor = Access::vertex_descriptor;
  using Point_2 = Access::Point_2;

  std::vector<Point> points = {
    Point(0.0,  1.0,  0.0),
    Point(0.0,  1.0,  1.0),
    Point(-1.0, 1.0,  1.0),
    Point(-1.0,-1.0,  1.0),
    Point(0.0, -1.0,  1.0),
    Point(0.0, -1.0,  0.0),
    Point(1.0, -1.0,  0.0),
    Point(1.0,  1.0,  0.0)
  };

  FT half_sqrt_2 = 0.7071067;

  std::vector<Vector_3> normals = {
    Vector_3(half_sqrt_2, 0.0, half_sqrt_2),
    Vector_3(half_sqrt_2, 0.0, half_sqrt_2),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(half_sqrt_2, 0.0, half_sqrt_2),
    Vector_3(half_sqrt_2, 0.0, half_sqrt_2),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  const auto candidates = find_ear_candidates_helper(points, normals);

  std::cout << "Found " << candidates.size() << " candidates:\n";
  for (const auto& c : candidates) {
    std::cout << "(" << c.first << "," << c.second << ")\n";
  }

  const double angleRad90 = 90 * (CGAL_PI / 180.0);
  const double angleRad180 = 180 * (CGAL_PI / 180.0);

  const std::vector<std::pair<double, std::size_t>> expected_candidates = {
    //std::pair<double, std::size_t> { angleRad180, 0 },
    //std::pair<double, std::size_t> { angleRad180, 1 },
    std::pair<double, std::size_t> { angleRad90, 2 },
    std::pair<double, std::size_t> { angleRad90, 3 },
    //std::pair<double, std::size_t> { angleRad180, 4 },
    //std::pair<double, std::size_t> { angleRad180, 5 },
    std::pair<double, std::size_t> { angleRad90, 6 },
    std::pair<double, std::size_t> { angleRad90, 7 },
  };



  TEST_ASSERT(candidates.size() == expected_candidates.size(), "Number of candidates does not match expected");
  for (std::size_t i = 0; i < candidates.size(); ++i) {
    bool found = false;
    for (std::size_t j = 0; j < expected_candidates.size(); ++j) {
      if (approxEq(candidates[i].first, expected_candidates[j].first) &&
          approxEq(candidates[i].second, expected_candidates[j].second)) {
        found = true;
        break;
      }
    }

    std::cout << "Checking candidate (" << candidates[i].first << "," << candidates[i].second << "): "
              << (found ? "found" : "not found") << "\n";
    TEST_ASSERT(found, "Candidate not found in expected candidates");
  }

  return true;
}

inline bool test_find_ear_candidates_wrapping_shape() {
  double n = -1.0;
  double p = 1.0;
  std::vector<Point> points = {
    Point(p,p,n),
    Point(p,p,p),
    Point(n,p,p),
    Point(n,n,p),
    Point(p,n,p),
    Point(p,n,n),
    Point(n,n,n),
    Point(n,p,n),
  };
  std::vector<Vector_3> normals = {
    normalize(Vector_3(1.0, 0.0, -1.0)),
    normalize(Vector_3(1.0, 0.0, 1.0)),
    normalize(Vector_3(0.0, 0.0, 1.0)),
    normalize(Vector_3(0.0, 0.0, 1.0)),
    normalize(Vector_3(1.0, 0.0, 1.0)),
    normalize(Vector_3(1.0, 0.0, -1.0)),
    normalize(Vector_3(0.0, 0.0, -1.0)),
    normalize(Vector_3(0.0, 0.0, -1.0))
  };


  const auto candidates = find_ear_candidates_helper(points, normals);

  std::cout << "Found " << candidates.size() << " candidates:\n";
  for (const auto& c : candidates) {
    std::cout << "(" << c.first << "," << c.second << ")\n";
  }

  const double angleRad90 = 90 * (CGAL_PI / 180.0);
  const double angleRad180 = 180 * (CGAL_PI / 180.0);

  const std::vector<std::pair<double, std::size_t>> expected_candidates = {
    //std::pair<double, std::size_t> { angleRad180, 0 },
    //std::pair<double, std::size_t> { angleRad180, 1 },
    std::pair<double, std::size_t> { angleRad90, 2 },
    std::pair<double, std::size_t> { angleRad90, 3 },
    //std::pair<double, std::size_t> { angleRad180, 4 },
    //std::pair<double, std::size_t> { angleRad180, 5 },
    std::pair<double, std::size_t> { angleRad90, 6 },
    std::pair<double, std::size_t> { angleRad90, 7 },
  };



  TEST_ASSERT(candidates.size() == expected_candidates.size(), "Candidate size mismatch");
  for (std::size_t i = 0; i < candidates.size(); ++i) {
    bool found = false;
    for (std::size_t j = 0; j < expected_candidates.size(); ++j) {
      if (approxEq(candidates[i].first, expected_candidates[j].first) &&
          approxEq(candidates[i].second, expected_candidates[j].second)) {
        found = true;
        break;
      }
    }

    std::cout << "Checking candidate (" << candidates[i].first << "," << candidates[i].second << "): "
              << (found ? "found" : "not found") << "\n";
    TEST_ASSERT(found, "Candidate not found in expected candidates");
  }

  return true;
}


inline void run_find_ear_candidates_tests(TestData& test_data) {
  test(test_find_ear_candidates_hexagon, "test_find_ear_candidates_hexagon", test_data);
  test(test_find_ear_candidates_step, "test_find_ear_candidates_step", test_data);
  test(test_find_ear_candidates_wrapping_shape, "test_find_ear_candidates_wrapping_shape", test_data);
}
