#pragma once
#include "common.h"

inline void add_edge(Polyhedron& mesh, vertex_descriptor v0, vertex_descriptor v1) {
  auto h = mesh.add_edge(v0, v1);
  auto h_opposite = opposite(h, mesh);

  // Set the targets explicitly.
  set_target(h, v1, mesh);
  set_target(h_opposite, v0, mesh);

  // Make the two halfedges form a boundary cycle.
  set_next(h, h_opposite, mesh);
  set_next(h_opposite, h, mesh);

  // Register an incident halfedge for each vertex.
  set_halfedge(v0, h_opposite, mesh);
  set_halfedge(v1, h, mesh);
}

inline bool test_compute_priority_3() {
  std::vector<Point> points = {
    Point(0.0, 0.0, 0.0),
    Point(1.0, 0.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  Polyhedron mesh;

  double edge_length = 0.5;

  // const FT padding = 1e-8;
  // Build the grid and insert points + normals.
  CGAL::Box_grid<Kernel> grid{FT(edge_length),
                              FT(-1.0), FT(-1.0), FT(-1.0),
                              FT(3.0), FT(3.0), FT(3.0)};
  
  grid.build(points, normals);
  grid.compute_block_normals();
  
  std::vector<vertex_descriptor> mesh_verts;
  for (const Point& p : points)
    mesh_verts.push_back(mesh.add_vertex(p));

  // Construct the reconstruction after initializing the mesh.
  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, normals);

  Candidate candidate {
    Point_3(0.5, 1.0, 0.0),
    Vector_3(0.0, 0.0, 1.0),
    0,
    mesh_verts[0],
    mesh_verts[1]
  };

  TEST_ASSERT(TestAccess::compute_priority(reconstruction, candidate) == 3, 
              "Priority should be 3 for the given candidate");

  return true;
}

inline bool test_compute_priority_2() {
  std::vector<Point> points = {
    Point(0.0, 0.0, 0.0),
    Point(1.0, 0.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  Polyhedron mesh;

  double edge_length = 0.5;

  // const FT padding = 1e-8;
  // Build the grid and insert points + normals.
  CGAL::Box_grid<Kernel> grid{FT(edge_length),
                              FT(-1.0), FT(-1.0), FT(-1.0),
                              FT(3.0), FT(3.0), FT(3.0)};
  
  grid.build(points, normals);
  grid.compute_block_normals();
  
  std::vector<vertex_descriptor> mesh_verts;
  for (const Point& p : points)
    mesh_verts.push_back(mesh.add_vertex(p));


  add_edge(mesh, mesh_verts[0], mesh_verts[1]);

  // Construct the reconstruction after initializing the mesh.
  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, normals);

  Candidate candidate {
    Point_3(0.5, 1.0, 0.0),
    Vector_3(0.0, 0.0, 1.0),
    0,
    mesh_verts[0],
    mesh_verts[1]
  };
  
  TEST_ASSERT(TestAccess::compute_priority(reconstruction, candidate) == 2, 
              "Priority should be 2 for the given candidate");

  return true;
}

inline bool test_compute_priority_1() {
  // two separate islands
  std::vector<Point> points = {
    Point(0.0, 0.0, 0.0),
    Point(1.0, 0.0, 0.0),
    Point(1.0, 1.0, 0.0),

    Point(2.0, 0.0, 0.0),
    Point(3.0, 0.0, 0.0),
    Point(2.0, 1.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),

    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  Polyhedron mesh;

  double edge_length = 0.5;

  // const FT padding = 1e-8;
  // Build the grid and insert points + normals.
  CGAL::Box_grid<Kernel> grid{FT(edge_length),
                              FT(-1.0), FT(-1.0), FT(-1.0),
                              FT(3.0), FT(3.0), FT(3.0)};
  
  grid.build(points, normals);
  grid.compute_block_normals();
  
  std::vector<vertex_descriptor> mesh_verts_face_0 {
    mesh.add_vertex(points[0]),
    mesh.add_vertex(points[1]),
    mesh.add_vertex(points[2])
  };
  
  std::vector<vertex_descriptor> mesh_verts_face_1 {
    mesh.add_vertex(points[3]),
    mesh.add_vertex(points[4]),
    mesh.add_vertex(points[5])
  };

  const auto face0 = mesh.add_face(mesh_verts_face_0);
  const auto face1 = mesh.add_face(mesh_verts_face_1);

  // Construct the reconstruction after initializing the mesh.
  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, normals);

  Candidate candidate {
    Point_3(0.5, 1.0, 0.0),
    Vector_3(0.0, 0.0, 1.0),
    0,
    mesh_verts_face_0[1],
    mesh_verts_face_1[0]
  };
  
  TEST_ASSERT(TestAccess::compute_priority(reconstruction, candidate) == 1, 
              "Priority should be 1 for the given candidate");

  return true;
}

inline bool test_compute_priority_0() {
  // single island
  std::vector<Point> points = {
    Point(0.0, 0.0, 0.0),
    Point(1.0, 0.0, 0.0),
    Point(1.0, 1.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  Polyhedron mesh;

  double edge_length = 0.5;

  // const FT padding = 1e-8;
  // Build the grid and insert points + normals.
  CGAL::Box_grid<Kernel> grid{FT(edge_length),
                              FT(-1.0), FT(-1.0), FT(-1.0),
                              FT(3.0), FT(3.0), FT(3.0)};
  
  grid.build(points, normals);
  grid.compute_block_normals();
  
  std::vector<vertex_descriptor> mesh_verts {
    mesh.add_vertex(points[0]),
    mesh.add_vertex(points[1]),
    mesh.add_vertex(points[2])
  };

  const auto face = mesh.add_face(mesh_verts);

  // Construct the reconstruction after initializing the mesh.
  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, normals);

  Candidate candidate {
    Point_3(0.5, 1.0, 0.0),
    Vector_3(0.0, 0.0, 1.0),
    0,
    mesh_verts[0],
    mesh_verts[1]
  };
  
  TEST_ASSERT(TestAccess::compute_priority(reconstruction, candidate) == 0, 
              "Priority should be 0 for the given candidate");

  return true;
}


void run_compute_priority_tests(TestData& test_data) {
  test(test_compute_priority_3, "test_compute_priority_3", test_data);
  test(test_compute_priority_2, "test_compute_priority_2", test_data);
  test(test_compute_priority_1, "test_compute_priority_1", test_data);
  test(test_compute_priority_0, "test_compute_priority_0", test_data);
}