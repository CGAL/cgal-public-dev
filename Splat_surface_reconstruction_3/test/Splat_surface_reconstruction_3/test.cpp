#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/splat_surface_reconstruction.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/IO/polygon_mesh_io.h>

#include <vector>
#include <fstream>
#include <string>
#include <filesystem>
#include <iostream>
#include <stdexcept>

#define TEST_ASSERT(expr, message) \
  do { \
    if (!(expr)) { \
      throw std::runtime_error(message); \
    } \
  } while (false)

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector_3;
typedef std::pair<Point, Vector_3> Pwn;
typedef CGAL::Surface_mesh<Point> Polyhedron;
typedef Kernel::FT FT;

typedef CGAL::Splat_surface_reconstruction_3<std::vector<Point>, std::vector<Vector_3>, Polyhedron> Splat_surface_reconstruction;

bool approxEq(double a, double b, double epsilon = 1e-6) {
  return std::abs(a - b) < epsilon;
}

namespace CGAL {

template <typename PointRange, typename NormalRange, typename TriangleMesh>
struct Splat_surface_reconstruction_3_test_access {
  using Reconstruction =
    CGAL::Splat_surface_reconstruction_3<PointRange, NormalRange, TriangleMesh>;

  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

  using Kernel = typename CGAL::Kernel_traits<
    typename PointRange::value_type>::Kernel;

  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;

  static Reconstruction make_reconstruction(
    const typename Reconstruction::Grid& grid,
    TriangleMesh& mesh,
    std::vector<Vector_3> normals)
  {
    return Reconstruction(grid, mesh, normals, 0);
  }

  static bool commit_ear_step(
    Reconstruction& reconstruction,
    std::vector<halfedge_descriptor>& cycle,
    std::size_t ear_pos)
  {
    return reconstruction.commit_ear_step(cycle, ear_pos);
  }


  static std::vector<std::pair<double, std::size_t>> find_ear_candidates(
    Reconstruction& reconstruction,
    const std::vector<halfedge_descriptor>& cycle,
    const std::vector<vertex_descriptor>& verts,
    const std::vector<Point_2>& P,
    const Vector_3& n,
    const Point_3& origin,
    const Vector_3& u,
    const Vector_3& v,
    const std::vector<int>& idx)
  {
    return reconstruction.find_ear_candidates(cycle, verts, P, n, origin, u, v, idx, false);
  }

  static bool ear_clip_and_add_faces(Reconstruction& reconstruction, std::vector<halfedge_descriptor> cycle) {
    return reconstruction.ear_clip_and_add_faces(cycle);
  }

  static Vector_3 compute_cycle_normal(Reconstruction& reconstruction, const std::vector<vertex_descriptor>& verts)
  {
    return reconstruction.compute_cycle_normal(verts);
  }

  static Point_2 project(Reconstruction& reconstruction, const Point_3& p, const Point_3& origin, const Vector_3& u, const Vector_3& v)
  {
    return reconstruction.project(p, origin, u, v);
  }

  static Polyhedron& mesh(Reconstruction& reconstruction)
  {
    return reconstruction.mesh_;
  }
};

}

using TestAccess =
  CGAL::Splat_surface_reconstruction_3_test_access<std::vector<Point>, std::vector<Vector_3>, Polyhedron>;
using Access = TestAccess;
using halfedge_descriptor = Access::halfedge_descriptor;
using vertex_descriptor = Access::vertex_descriptor;
using Point_2 = Access::Point_2;

bool test_commit_ear_step() {
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

Polyhedron ear_clip_and_add_faces_helper(std::vector<Point> points, std::vector<Vector_3> normals) {
  
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

  auto h = halfedge(face, mesh);
  const auto start = h;

  do {
    cycle.push_back(h);
    verts.push_back(source(h, mesh));
    
    h = next(h, mesh);
  } while (h != start);

  const auto null_f = boost::graph_traits<Polyhedron>::null_face();

  for (halfedge_descriptor h : cycle) {
    set_face(h, null_f, mesh);
  }

  remove_face(face, mesh);

  // Verify that there are no faces yet
  assert(mesh.number_of_faces() == 0); // No faces should exist

  TestAccess::ear_clip_and_add_faces(reconstruction, cycle);

  return TestAccess::mesh(reconstruction);
}

bool test_ear_clip_and_add_faces_single_tri() {
  std::vector<Point> points = {
    Point(1.0, 1.0, 0.0),
    Point(-1.0, -1.0, 0.0),
    Point(1.0, -1.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  Polyhedron mesh = ear_clip_and_add_faces_helper(points, normals);

  assert(num_vertices(mesh) == 3);
  assert(num_edges(mesh) == 3);
  assert(mesh.number_of_faces() == 1);

  for (auto f : faces(mesh)) {
    std::size_t face_degree = 0;
    for (auto h : halfedges_around_face(halfedge(f, mesh), mesh))
      ++face_degree;

    assert(face_degree == 3);
  }

  assert(CGAL::is_valid_polygon_mesh(mesh));

  return true;
}

bool test_ear_clip_and_add_faces_square() {
  std::vector<Point> points = {
    Point(-1.0, 1.0, 0.0),
    Point(-1.0, -1.0, 0.0),
    Point(1.0, -1.0, 0.0),
    Point(1.0, 1.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  Polyhedron mesh = ear_clip_and_add_faces_helper(points, normals);

  std::cout << mesh << std::endl;

  assert(num_vertices(mesh) == 4);
  assert(num_edges(mesh) == 5);
  assert(mesh.number_of_faces() == 2);

  for (auto f : faces(mesh)) {
    std::size_t face_degree = 0;
    for (auto h : halfedges_around_face(halfedge(f, mesh), mesh))
      ++face_degree;

    assert(face_degree == 3);
  }

  assert(CGAL::is_valid_polygon_mesh(mesh));

  return true;
}


bool test_ear_clip_and_add_faces_wrapping_shape() {

  // This test checks the ear clipping and face addition for a shape that wraps around, 
  // and thus has opposing normals
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

  Polyhedron mesh = ear_clip_and_add_faces_helper(points, normals);

  std::cout << mesh << std::endl;

  assert(num_vertices(mesh) == 4);
  assert(num_edges(mesh) == 5);
  assert(mesh.number_of_faces() == 2);

  for (auto f : faces(mesh)) {
    std::size_t face_degree = 0;
    for (auto h : halfedges_around_face(halfedge(f, mesh), mesh))
      ++face_degree;

    assert(face_degree == 3);
  }

  assert(CGAL::is_valid_polygon_mesh(mesh));

  return true;
}

bool test_find_ear_candidates_hexagon() {
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



  assert(candidates.size() == expected_candidates.size());
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
    assert(found);
  }

  return true;
}

bool test_find_ear_candidates_step()
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



  assert(candidates.size() == expected_candidates.size());
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
    assert(found);
  }

  return true;
}

bool test_find_ear_candidates_wrapping_shape() {
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
    assert(found);
  }

  return true;
}

bool test_circle_splat_intersection_exact()
{
  using Grid = CGAL::Box_grid<Kernel>;

  // Create a dummy grid.
  Grid grid(
      FT(1.0),
      FT(-2.0), FT(-2.0), FT(-2.0),
      FT(2.0), FT(2.0), FT(2.0));

  const Point parent_a(-0.5, 0.0, 0.0);
  const Point parent_b( 0.5, 0.0, 0.0);

  const Point splat_center(0.0, 0.0, 0.0);

  // Plane y = 0
  const Vector_3 splat_normal(0.0, 1.0, 0.0);

  const FT splat_radius = FT(10.0);

  auto hits =
      grid.intersect_circle_with_splat(
          parent_a,
          parent_b,
          splat_center,
          splat_normal,
          splat_radius);

  std::cout << "\n========== Circle-Splat Test ==========\n";

  std::cout << "Expected intersections:\n";
  std::cout << "(0,0,+sqrt(3)/2)\n";
  std::cout << "(0,0,-sqrt(3)/2)\n";

  std::cout << "\nFound "
            << hits.size()
            << " intersections:\n";

  for (std::size_t i = 0; i < hits.size(); ++i) {
    std::cout
        << hits[i].x() << " "
        << hits[i].y() << " "
        << hits[i].z() << "\n";
  }

  if (hits.size() != 2) {
    std::cerr << "FAILED: expected exactly 2 points.\n";
    return false;
  }

  const double expected =
      std::sqrt(3.0) / 2.0;

  const double z0 =
      CGAL::to_double(hits[0].z());

  const double z1 =
      CGAL::to_double(hits[1].z());

  const double eps = 1e-6;

  bool pass =
      std::abs(std::abs(z0) - expected) < eps &&
      std::abs(std::abs(z1) - expected) < eps;



  if (pass) {
    std::cout << "PASSED\n";
  } else {
    std::cout << "FAILED\n";
    std::cout << "Expected |z| = "
              << expected
              << std::endl;
  }

  std::cout << "\n========== No-Intersection Test ==========\n";

  const Point far_splat_center(0.0, 2.0, 0.0);
  const Vector_3 far_splat_normal(0.0, 1.0, 0.0);

  auto no_hits =
      grid.intersect_circle_with_splat(
          parent_a,
          parent_b,
          far_splat_center,
          far_splat_normal,
          splat_radius);

  std::cout << "Found "
            << no_hits.size()
            << " intersections.\n";

  if (!no_hits.empty()) {
    std::cerr
        << "FAILED: expected zero intersections."
        << std::endl;

    for (std::size_t i = 0; i < no_hits.size(); ++i) {
      std::cout
          << no_hits[i].x() << " "
          << no_hits[i].y() << " "
          << no_hits[i].z() << "\n";
    }

    return false;
  }

  std::cout << "PASSED\n";

  return pass;
}

struct TestData {
  uint total_tests = 0;
  uint passed_tests = 0;
  uint failed_tests = 0;
  std::vector<std::string> failed_test_names;
};

template <typename Function>
void test(Function&& function, const char* test_name, TestData& test_data) {
  test_data.total_tests++;
  try {
    std::forward<Function>(function)();
    test_data.passed_tests++;
  }
  catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    test_data.failed_tests++;
    test_data.failed_test_names.push_back(test_name);
  }
  catch (...) {
    std::cout << "Unknown exception caught in test " << test_name << '\n';
    test_data.failed_tests++;
    test_data.failed_test_names.push_back(test_name);
  }
}

int main(int argc, char* argv[]) {
  test_circle_splat_intersection_exact();

  TestData test_data;

  test(test_find_ear_candidates_hexagon, "test_find_ear_candidates_hexagon", test_data);

  test(test_find_ear_candidates_step, "test_find_ear_candidates_step", test_data);

  test(test_find_ear_candidates_wrapping_shape, "test_find_ear_candidates_wrapping_shape", test_data);

  test(test_ear_clip_and_add_faces_square, "test_ear_clip_and_add_faces_square", test_data);
  test(test_ear_clip_and_add_faces_single_tri, "test_ear_clip_and_add_faces_single_tri", test_data);

  std::cout << "\n========== Test Summary ==========\n";
  std::cout << "Passed " << test_data.passed_tests << " out of " << test_data.total_tests << " tests\n";

  if (!test_data.failed_test_names.empty()) {
    std::cout << "\nFailed tests:\n";
    for (const auto& name : test_data.failed_test_names) {
      std::cout << "  " << name << "\n";
    }
  }

  std::cout << std::endl;

  return EXIT_SUCCESS;
}