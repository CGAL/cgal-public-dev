#pragma once
#include "common.h"


inline bool test_has_vertex_near_one_vertex() {
  std::vector<Point> points = {
    Point(0.0, 0.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0)
  };

  Polyhedron mesh;

  double average_spacing = CGAL::compute_average_spacing<CGAL::Parallel_if_available_tag>(points, 2);
  std::cout << "Average spacing: " << average_spacing << std::endl;

  double edge_length = 0.5;

  const auto bbox = CGAL::bounding_box(points.begin(), points.end()); // recompute bbox after centering and scaling

  // const FT padding = 1e-8;
  // Build the grid and insert points + normals.
  CGAL::Box_grid<Kernel> grid{FT(edge_length),
                              FT(-1.0), FT(-1.0), FT(-1.0),
                              FT(3.0), FT(3.0), FT(3.0)}; // initialize grid with cell size equal to average spacing and bounding box [-1,1]^3
  
  grid.build(points, normals);
  grid.compute_block_normals();
  
  std::vector<vertex_descriptor> mesh_verts;
  for (const Point& p : points)
    mesh_verts.push_back(mesh.add_vertex(p));

  // Construct the reconstruction after initializing the mesh.
  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, normals);

  FT d = FT(edge_length);

  for (const auto& v : mesh_verts) {
    TestAccess::insert_mesh_vertex(reconstruction, v);
  }

  TEST_ASSERT(
      TestAccess::has_vertex_near(reconstruction, Point(d*0.99, 0.0, 0.0)),
      "The only vertex should be near the query point");
  TEST_ASSERT(
      !TestAccess::has_vertex_near(reconstruction, Point(d*1.01, 0.0, 0.0)),
      "No other vertices should be near the query point");
  TEST_ASSERT(
      !TestAccess::has_vertex_near(reconstruction, Point(d, 0.0, 0.0)),
      "No other vertices should be near the query point");

  return true;
}

void run_has_vertex_near_tests(TestData& test_data) {
    test(test_has_vertex_near_one_vertex, "test_has_vertex_near_one_vertex", test_data);
}