#pragma once
#include "common.h"

Splat_surface_reconstruction get_reconstruction() {
  Polyhedron mesh;

  double edge_length = 0.5;

  // const FT padding = 1e-8;
  // Build the grid and insert points + normals.
  CGAL::Box_grid<Kernel> grid{FT(edge_length),
                              FT(-1.0), FT(-1.0), FT(-1.0),
                              FT(3.0), FT(3.0), FT(3.0)};

  // Construct the reconstruction after initializing the mesh.
  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, std::vector<Vector_3>{});
  return reconstruction;
}

Candidate create_candidate() {
  Candidate cand(Point_3(0,0,0), CGAL::NULL_VECTOR, 0, vertex_descriptor(), vertex_descriptor()); // Dummy initialization
  return cand;
}

Candidate create_candidate(int priority, Index splat_id) {
  Candidate cand = create_candidate();
  cand.priority = priority;
  cand.splat_id = splat_id;
  return cand;
}

bool test_candidate_queue_initially_empty() {
  auto reconstruction = get_reconstruction();

  TEST_ASSERT(
    TestAccess::candidate_queues_empty(reconstruction),
    "Queue should initially be empty");

  return true;
}

bool test_candidate_queue_pop_empty() {
  auto reconstruction = get_reconstruction();
  Candidate popped = create_candidate();

  TEST_ASSERT(
    !TestAccess::pop_candidate_from_queue(reconstruction, popped),
    "Popping an empty queue should fail");

  TEST_ASSERT(
    TestAccess::candidate_queues_empty(reconstruction),
    "Queue should remain empty after a failed pop");

  return true;
}

bool test_candidate_queue_single_element() {
  auto reconstruction = get_reconstruction();

  constexpr int priority = 2;
  constexpr Index splat_id = 42;

  TestAccess::push_candidate_in_queue(
    reconstruction, create_candidate(priority, splat_id));

  Candidate popped = create_candidate();

  TEST_ASSERT(
    TestAccess::pop_candidate_from_queue(reconstruction, popped),
    "Failed to pop the only candidate");

  TEST_ASSERT(popped.priority == priority,
              "Incorrect priority returned");
  TEST_ASSERT(popped.splat_id == splat_id,
              "Incorrect splat_id returned");
  TEST_ASSERT(TestAccess::candidate_queues_empty(reconstruction),
              "Queue should be empty after popping its only candidate");

  return true;
}

// Test all four valid priorities and their ordering.
bool test_candidate_queue_priority_order() {
  auto reconstruction = get_reconstruction();

  TestAccess::push_candidate_in_queue(reconstruction, create_candidate(0, 10));
  TestAccess::push_candidate_in_queue(reconstruction, create_candidate(3, 11));
  TestAccess::push_candidate_in_queue(reconstruction, create_candidate(1, 12));
  TestAccess::push_candidate_in_queue(reconstruction, create_candidate(2, 13));

  const int expected_priorities[] = {3, 2, 1, 0};
  const Index expected_ids[] = {11, 13, 12, 10};

  for (int i = 0; i < 4; ++i) {
    Candidate popped = create_candidate();

    TEST_ASSERT(
      TestAccess::pop_candidate_from_queue(reconstruction, popped),
      "Failed to pop a candidate");

    TEST_ASSERT(
      popped.priority == expected_priorities[i],
      "Candidates were not popped in descending priority order");

    TEST_ASSERT(
      popped.splat_id == expected_ids[i],
      "Candidate identity does not match its priority");
  }

  TEST_ASSERT(
    TestAccess::candidate_queues_empty(reconstruction),
    "Queue should be empty after popping all candidates");

  return true;
}

// Ensure candidates sharing a priority are all preserved.
bool test_candidate_queue_duplicate_priorities() {
  auto reconstruction = get_reconstruction();

  TestAccess::push_candidate_in_queue(reconstruction, create_candidate(2, 101));
  TestAccess::push_candidate_in_queue(reconstruction, create_candidate(2, 102));
  TestAccess::push_candidate_in_queue(reconstruction, create_candidate(2, 103));

  bool seen[3] = {false, false, false};

  for (int i = 0; i < 3; ++i) {
    Candidate popped = create_candidate();

    TEST_ASSERT(
      TestAccess::pop_candidate_from_queue(reconstruction, popped),
      "Failed to pop a candidate with a duplicate priority");

    TEST_ASSERT(popped.priority == 2,
                "Incorrect priority for duplicate-priority candidate");

    TEST_ASSERT(
      popped.splat_id >= 101 && popped.splat_id <= 103,
      "Unexpected splat_id returned");

    const int index = static_cast<int>(popped.splat_id) - 101;

    TEST_ASSERT(!seen[index],
                "A candidate was returned more than once");

    seen[index] = true;
  }

  TEST_ASSERT(seen[0] && seen[1] && seen[2],
              "Not all unique candidates were returned");

  TEST_ASSERT(TestAccess::candidate_queues_empty(reconstruction),
              "Queue should be empty after popping all candidates");

  return true;
}

// Exercise repeated push/pop operations.
bool test_candidate_queue_interleaved_operations() {
  auto reconstruction = get_reconstruction();
  Candidate popped = create_candidate();

  TestAccess::push_candidate_in_queue(reconstruction, create_candidate(1, 201));
  TestAccess::push_candidate_in_queue(reconstruction, create_candidate(0, 202));

  TEST_ASSERT(
    TestAccess::pop_candidate_from_queue(reconstruction, popped),
    "Failed to pop first candidate");
  TEST_ASSERT(popped.splat_id == 201,
              "Expected priority-1 candidate first");

  TestAccess::push_candidate_in_queue(reconstruction, create_candidate(3, 203));
  TestAccess::push_candidate_in_queue(reconstruction, create_candidate(2, 204));

  const int expected_priorities[] = {3, 2, 0};
  const Index expected_ids[] = {203, 204, 202};

  for (int i = 0; i < 3; ++i) {
    TEST_ASSERT(
      TestAccess::pop_candidate_from_queue(reconstruction, popped),
      "Failed to pop candidate after interleaved operations");

    TEST_ASSERT(popped.priority == expected_priorities[i],
                "Incorrect priority after interleaved operations");
    TEST_ASSERT(popped.splat_id == expected_ids[i],
                "Incorrect candidate identity after interleaved operations");
  }

  TEST_ASSERT(TestAccess::candidate_queues_empty(reconstruction),
              "Queue should be empty at the end");

  return true;
}

// Verify that the queue can be reused after becoming empty.
bool test_candidate_queue_reuse() {
  auto reconstruction = get_reconstruction();
  Candidate popped = create_candidate();

  TestAccess::push_candidate_in_queue(reconstruction, create_candidate(3, 301));

  TEST_ASSERT(
    TestAccess::pop_candidate_from_queue(reconstruction, popped),
    "Failed to pop initial candidate");
  TEST_ASSERT(popped.splat_id == 301,
              "Incorrect initial candidate");

  TEST_ASSERT(TestAccess::candidate_queues_empty(reconstruction),
              "Queue should be empty before reuse");

  TestAccess::push_candidate_in_queue(reconstruction, create_candidate(0, 302));

  TEST_ASSERT(
    TestAccess::pop_candidate_from_queue(reconstruction, popped),
    "Failed to pop candidate after reuse");
  TEST_ASSERT(popped.priority == 0 && popped.splat_id == 302,
              "Incorrect candidate returned after reuse");

  TEST_ASSERT(TestAccess::candidate_queues_empty(reconstruction),
              "Queue should be empty after reuse");

  return true;
}

void run_candidate_queue_tests(TestData& test_data) {
  test(test_candidate_queue_initially_empty,
       "test_candidate_queue_initially_empty", test_data);
  test(test_candidate_queue_pop_empty,
       "test_candidate_queue_pop_empty", test_data);
  test(test_candidate_queue_single_element,
       "test_candidate_queue_single_element", test_data);
  test(test_candidate_queue_priority_order,
       "test_candidate_queue_priority_order", test_data);
  test(test_candidate_queue_duplicate_priorities,
       "test_candidate_queue_duplicate_priorities", test_data);
  test(test_candidate_queue_interleaved_operations,
       "test_candidate_queue_interleaved_operations", test_data);
  test(test_candidate_queue_reuse,
       "test_candidate_queue_reuse", test_data);
}