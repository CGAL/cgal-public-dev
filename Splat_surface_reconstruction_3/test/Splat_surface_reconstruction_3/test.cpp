// Main test runner: each function under test has its tests in tests/test_<function>.h
#include "tests/test_commit_ear_step.h"
#include "tests/test_to_grid_coordinates.h"
#include "tests/test_flat_index.h"
#include "tests/test_cell_center.h"
#include "tests/test_splat_intersects_box.h"
#include "tests/test_build_cell_point_associations.h"
#include "tests/test_compute_block_normals.h"
#include "tests/test_ear_clip_and_add_faces.h"
#include "tests/test_find_ear_candidates.h"
#include "tests/test_circle_splat_intersection.h"

int main() {
  TestData test_data;

  run_commit_ear_step_tests(test_data);
  run_to_grid_coordinates_tests(test_data);
  run_flat_index_tests(test_data);
  run_cell_center_tests(test_data);
  run_splat_intersects_box_tests(test_data);
  run_build_cell_point_associations_tests(test_data);
  run_compute_block_normals_tests(test_data);
  run_ear_clip_and_add_faces_tests(test_data);
  run_find_ear_candidates_tests(test_data);
  run_circle_splat_intersection_tests(test_data);

  std::cout << "\n========== Test Summary ==========\n";
  std::cout << "Passed " << test_data.passed_tests << " out of " << test_data.total_tests << " tests\n";

  if (!test_data.failed_test_names.empty()) {
    std::cout << "\nFailed tests:\n";
    for (const auto& name : test_data.failed_test_names) {
      std::cout << "  " << name << "\n";
    }
  }

  std::cout << std::endl;

  return test_data.failed_tests == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
