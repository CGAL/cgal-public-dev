// Main test runner: each function under test has its tests in tests/test_<function>.h
#include "tests/test_commit_ear_step.h"
#include "tests/test_to_grid_coordinates.h"
#include "tests/test_flat_index.h"
#include "tests/test_cell_center.h"
#include "tests/test_build_cell_point_associations.h"
#include "tests/test_compute_block_normals.h"
#include "tests/test_ear_clip_and_add_faces.h"
#include "tests/test_find_ear_candidates.h"
#include "tests/test_circle_splat_intersection.h"
#include "tests/test_nearby_point_ids.h"
#include "tests/test_has_vertex_near.h"
#include "tests/test_compute_priority.h"
#include "tests/test_initialize_grid.h"
#include "tests/test_intersect_2d.h"
#include "tests/test_candidate_queue.h"
#include "tests/test_intersect_circle_with_splat.h"
#include "tests/test_estimate_individual_splat_sizes.h"
#include "tests/test_generate_candidates_from_parents.h"



int main() {
  TestData test_data;

  run_commit_ear_step_tests(test_data);
  run_to_grid_coordinates_tests(test_data);
  run_flat_index_tests(test_data);
  run_cell_center_tests(test_data);
  run_build_cell_point_associations_tests(test_data);
  run_compute_block_normals_tests(test_data);
  run_ear_clip_and_add_faces_tests(test_data);
  run_find_ear_candidates_tests(test_data);
  run_circle_splat_intersection_tests(test_data);
  run_nearby_point_ids_tests(test_data);
  run_has_vertex_near_tests(test_data);
  run_compute_priority_tests(test_data);
  run_initialize_grid_tests(test_data);
  run_intersect_2d_tests(test_data);
  run_candidate_queue_tests(test_data);
  run_intersect_circle_with_splat_tests(test_data);
  run_estimate_individual_splat_sizes_tests(test_data);
  run_generate_candidates_from_parents_tests(test_data);

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
