#pragma once
#include "common.h"

inline bool test_commit_ear_step() {
  // Implement the test for commit_ear_step here.
  return true;
}

Vector_3 normalize(const Vector_3& v)
{
  return v / std::sqrt(v.squared_length());
}

std::vector<std::pair<double, std::size_t>> find_ear_candidates_helper(std::vector<Point> points, std::vector<Vector_3> normals) {
  

  Polyhedron mesh;

  double average_spacing = CGAL::compute_average_spacing<CGAL::Parallel_if_available_tag>(points, 2);
  std::cout << "Average spacing: " << average_spacing << std::endl;

  double edge_length = average_spacing;

  const auto bbox = CGAL::bounding_box(points.begin(), points.end()); // recompute bbox after centering and scaling


  const FT padding = FT(1.1) * edge_length;

  // const FT padding = 1e-8;
  // Build the grid and insert points + normals.
  CGAL::Box_grid<Kernel> grid{FT(edge_length),
                              FT(bbox.xmin()-padding), FT(bbox.ymin()-padding), FT(bbox.zmin()-padding),
                              FT(bbox.xmax()+padding), FT(bbox.ymax()+padding), FT(bbox.zmax()+padding)}; // initialize grid with cell size equal to average spacing and bounding box [-1,1]^3
  
  grid.build(points, normals);
  grid.compute_block_normals();
  
  std::vector<vertex_descriptor> mesh_verts;
  for (const Point& p : points)
    mesh_verts.push_back(mesh.add_vertex(p));

  const auto face = mesh.add_face(mesh_verts);
  assert(face != Polyhedron::null_face());

  // Construct the reconstruction after initializing the mesh.
  Splat_surface_reconstruction reconstruction =
    TestAccess::make_reconstruction(grid, mesh, normals);

  // Extract the halfedge cycle and vertices in matching cyclic order.
  std::vector<halfedge_descriptor> cycle;
  std::vector<vertex_descriptor> verts;

  auto h = halfedge(mesh_verts[0], mesh_verts[1], mesh).first;
  const auto start = h;

  do {
    cycle.push_back(h);
    verts.push_back(source(h, mesh));
    h = next(h, mesh);
  } while (h != start);

  // The hexagon lies in the xy-plane, centered at the origin.
  const Point origin = points[0];
  Vector_3 n = Access::compute_cycle_normal(reconstruction, verts);
  std::vector<Vector_3> frame = grid.compute_local_tangent_frame(n);
  Vector_3 u = frame[0];
  Vector_3 v = frame[1];
  CGAL_assertion(u != CGAL::NULL_VECTOR && v != CGAL::NULL_VECTOR);
  if (CGAL::cross_product(u, v) * n < FT(0))
    v = -v;

  // Project each vertex into the local 2D coordinate system.
  std::vector<Point_2> P;
  P.reserve(verts.size());

  for (Point p : points)
  {
    P.push_back(TestAccess::project(reconstruction, p, origin, u, v));
  }

  // Associate each cycle vertex with its corresponding input index.
  std::vector<int> idx;
  idx.reserve(verts.size());

  for (const auto vertex : verts) {
    idx.push_back(static_cast<int>(
      get(CGAL::vertex_index, mesh, vertex)));
    std::cout << "Vertex " << vertex << " has index " << idx.back() << std::endl;
  }

  const auto candidates = Access::find_ear_candidates(
    reconstruction, cycle, verts, P, n, origin, u, v, idx);

  return candidates;
}


inline void run_commit_ear_step_tests(TestData& test_data) {
  test(test_commit_ear_step, "test_commit_ear_step", test_data);
}
