#pragma once
#include "common.h"

// Creates a candidate with parents (0,0,0) and (1,0,0),
// and one additional edge q0-q1. Checks the projection result.
inline bool run_projection_check_edge_case(
    const Point_3& q0,
    const Point_3& q1,
    const Vector_3& normal,
    bool expected,
    const char* message)
{
  FT d = FT(1.0);

  CGAL::Box_grid<Kernel> grid{
    d, -3.5, -3.5, -3.5, 3.5, 3.5, 3.5
  };

  std::vector<Point_3> points = {
    Point_3(0.0, 0.0, 0.0), // parent 0
    Point_3(1.0, 0.0, 0.0), // parent 1
    q0,
    q1
  };

  std::vector<Vector_3> normals(
    points.size(), Vector_3(0.0, 0.0, 1.0));

  Polyhedron mesh;
  std::vector<vertex_descriptor> point_descriptors;

  for (const Point_3& p : points)
    point_descriptors.push_back(mesh.add_vertex(p));

  add_edge_to_mesh(mesh, point_descriptors[0], point_descriptors[1]);
  add_edge_to_mesh(mesh, point_descriptors[2], point_descriptors[3]);

  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, normals);

  for (const vertex_descriptor& vd : point_descriptors)
    TestAccess::insert_mesh_vertex(reconstruction, vd);

  Candidate cand(
    Point_3(1.0, 1.0, 0.0),
    CGAL::NULL_VECTOR,
    0,
    point_descriptors[0],
    point_descriptors[1]);

  bool result =
    TestAccess::projection_check(reconstruction, cand, normal);

  TEST_ASSERT(result == expected, message);
  return true;
}

// ------------------------------------------------------------
// Basic projection geometry
// ------------------------------------------------------------

inline bool test_projection_check_no_obstruction()
{
  return run_projection_check_edge_case(
    Point_3(-2.0, -2.0, 0.0),
    Point_3(-2.5, -2.0, 0.0),
    Vector_3(0.0, 0.0, 1.0),
    true,
    "A non-intersecting edge should not obstruct the candidate");
}

inline bool test_projection_check_edge_in_the_way()
{
  return run_projection_check_edge_case(
    Point_3(1.0, 0.0, 0.0),
    Point_3(0.0, 1.0, 0.0),
    Vector_3(0.0, 0.0, 1.0),
    false,
    "An edge crossing a candidate-parent edge in projection "
    "should obstruct the candidate");
}

inline bool test_projection_check_intersection_only_after_projection()
{
  return run_projection_check_edge_case(
    Point_3(1.0, 0.0, 0.5),
    Point_3(0.0, 1.0, 0.5),
    Vector_3(0.0, 0.0, 1.0),
    false,
    "An eligible edge should obstruct the candidate when its "
    "projection intersects, even if the 3D edges do not");
}

inline bool test_projection_check_projected_near_miss()
{
  return run_projection_check_edge_case(
    Point_3(0.2, 0.8, 0.0),
    Point_3(0.2, 1.2, 0.0),
    Vector_3(0.0, 0.0, 1.0),
    true,
    "A projected edge that misses the candidate-parent edges "
    "should not obstruct the candidate");
}

inline bool test_projection_check_collinear_overlap()
{
  return run_projection_check_edge_case(
    Point_3(0.2, 0.2, 0.0),
    Point_3(0.8, 0.8, 0.0),
    Vector_3(0.0, 0.0, 1.0),
    false,
    "Collinear overlapping projected edges should be rejected");
}

// ------------------------------------------------------------
// Cylinder inclusion: radius
// ------------------------------------------------------------

inline bool test_projection_check_both_vertices_outside_radius()
{
  // Both endpoints are outside radius 2d from candidate (1,1,0).
  // Their connecting edge crosses the projected candidate-parent
  // diagonal, but neither endpoint should be selected as nearby.
  return run_projection_check_edge_case(
    Point_3(-2.0, 2.5, 0.0),
    Point_3(2.5, -2.0, 0.0),
    Vector_3(0.0, 0.0, 1.0),
    true,
    "An edge whose endpoints are both outside the cylinder radius "
    "should not obstruct the candidate");
}

inline bool test_projection_check_one_vertex_inside_radius()
{
  // q0 is inside the radius; q1 is outside.
  // The edge crosses the projected candidate-parent diagonal.
  // Because q0 is selected, its incident edge should be checked.
  return run_projection_check_edge_case(
    Point_3(0.0, 1.0, 0.0),
    Point_3(3.0, 0.0, 0.0),
    Vector_3(0.0, 0.0, 1.0),
    false,
    "An edge incident to a vertex inside the cylinder should be "
    "checked even if its other endpoint is outside");
}

// ------------------------------------------------------------
// Cylinder inclusion: height along the normal
// ------------------------------------------------------------

inline bool test_projection_check_both_vertices_behind_cylinder()
{
  // Both vertices are beyond the back cap, at z = -1.01.
  return run_projection_check_edge_case(
    Point_3(1.0, 0.0, -1.01),
    Point_3(0.0, 1.0, -1.01),
    Vector_3(0.0, 0.0, 1.0),
    true,
    "An edge with both endpoints behind the cylinder should "
    "not obstruct the candidate");
}

inline bool test_projection_check_both_vertices_in_front_of_cylinder()
{
  // Both vertices are beyond the front cap, at z = 1.01.
  return run_projection_check_edge_case(
    Point_3(1.0, 0.0, 1.01),
    Point_3(0.0, 1.0, 1.01),
    Vector_3(0.0, 0.0, 1.0),
    true,
    "An edge with both endpoints in front of the cylinder should "
    "not obstruct the candidate");
}

inline bool test_projection_check_one_vertex_outside_height()
{
  // q0 is at z = 0 and is inside the cylinder.
  // q1 is beyond the front cap.
  return run_projection_check_edge_case(
    Point_3(1.0, 0.0, 0.0),
    Point_3(0.0, 1.0, 1.01),
    Vector_3(0.0, 0.0, 1.0),
    false,
    "An edge incident to a vertex inside the cylinder should be "
    "checked even if its other endpoint is outside its height");
}

// ------------------------------------------------------------
// Normal direction
// ------------------------------------------------------------

inline bool test_projection_check_reversed_normal()
{
  // Both endpoints are in front of the candidate for +Z.
  // With the normal reversed, they are behind the candidate.
  // They remain outside the cylinder's height range.
  return run_projection_check_edge_case(
    Point_3(1.0, 0.0, 1.01),
    Point_3(0.0, 1.0, 1.01),
    Vector_3(0.0, 0.0, -1.0),
    true,
    "Vertices outside the cylinder should remain excluded "
    "when the normal is reversed");
}
inline bool test_projection_check_reversed_normal_inside()
{
  // Both endpoints are in front of the candidate for +Z.
  // With the normal reversed, they are behind the candidate.
  // They remain outside the cylinder's height range.
  return run_projection_check_edge_case(
    Point_3(1.0, 0.0, 0.99),
    Point_3(0.0, 1.0, 0.99),
    Vector_3(0.0, 0.0, -1.0),
    false,
    "Vertices outside the cylinder should remain included "
    "when the normal is reversed");
}

void run_projection_check_tests(TestData& test_data)
{
  test(test_projection_check_no_obstruction,
       "test_projection_check_no_obstruction", test_data);

  test(test_projection_check_edge_in_the_way,
       "test_projection_check_edge_in_the_way", test_data);

  test(test_projection_check_intersection_only_after_projection,
       "test_projection_check_intersection_only_after_projection", test_data);

  test(test_projection_check_projected_near_miss,
       "test_projection_check_projected_near_miss", test_data);

  test(test_projection_check_collinear_overlap,
       "test_projection_check_collinear_overlap", test_data);

  test(test_projection_check_both_vertices_outside_radius,
       "test_projection_check_both_vertices_outside_radius", test_data);

  test(test_projection_check_one_vertex_inside_radius,
       "test_projection_check_one_vertex_inside_radius", test_data);

  test(test_projection_check_both_vertices_behind_cylinder,
       "test_projection_check_both_vertices_behind_cylinder", test_data);

  test(test_projection_check_both_vertices_in_front_of_cylinder,
       "test_projection_check_both_vertices_in_front_of_cylinder", test_data);

  test(test_projection_check_one_vertex_outside_height,
       "test_projection_check_one_vertex_outside_height", test_data);

  test(test_projection_check_reversed_normal,
       "test_projection_check_reversed_normal", test_data);
  
  test(test_projection_check_reversed_normal_inside,
       "test_projection_check_reversed_normal_inside", test_data);
}