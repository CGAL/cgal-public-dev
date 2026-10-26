#pragma once
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
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef std::pair<Point, Vector_3> Pwn;
typedef CGAL::Surface_mesh<Point> Polyhedron;
typedef Kernel::FT FT;

typedef CGAL::Splat_surface_reconstruction_3<std::vector<Point>, std::vector<Vector_3>, Polyhedron> Splat_surface_reconstruction;

inline bool approxEq(double a, double b, double epsilon = 1e-6) {
  return std::abs(a - b) < epsilon;
}

inline bool approxEq(const Point_3& p1, const Point_3& p2, double epsilon = 1e-6) {
  return approxEq(p1.x(), p2.x(), epsilon) &&
         approxEq(p1.y(), p2.y(), epsilon) &&
         approxEq(p1.z(), p2.z(), epsilon);
}

inline bool approxEq(const Vector_3& v1, const Vector_3& v2, double epsilon = 1e-6) {
  return approxEq(v1.x(), v2.x(), epsilon) &&
         approxEq(v1.y(), v2.y(), epsilon) &&
         approxEq(v1.z(), v2.z(), epsilon);
}

inline bool approxEq(const std::vector<Vector_3>& v1, const std::vector<Vector_3>& v2, double epsilon = 1e-6) {
  if (v1.size() != v2.size()) return false;
  for (std::size_t i = 0; i < v1.size(); ++i) {
    if (!approxEq(v1[i], v2[i], epsilon)) return false;
  }
  return true;
}

// Get a random integer between min and max (inclusive)
inline int rand_int(int min, int max) {
  return min + std::rand() % (max - min + 1);
}

// Get a random float in [min, max]
inline float rand_float_inclusive(float min, float max) {
  return min + (max - min) *
    (static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));
}

// Get a random float in [min, max)
inline float rand_float_exclusive(float min, float max) {
  return min + (max - min) *
    (static_cast<float>(std::rand()) / (static_cast<float>(RAND_MAX) + 1.0f));
}

inline Point_3 rand_point_3(float min_x, float max_x, float min_y, float max_y, float min_z, float max_z) {
  return Point_3(
    rand_float_inclusive(min_x, max_x),
    rand_float_inclusive(min_y, max_y),
    rand_float_inclusive(min_z, max_z)
  );
}

inline std::vector<Point_3> rand_points_3(std::size_t count, float min_x, float max_x, float min_y, float max_y, float min_z, float max_z) {
  std::vector<Point_3> points;
  points.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    points.push_back(rand_point_3(min_x, max_x, min_y, max_y, min_z, max_z));
  }
  return points;
}



struct TestData {
  uint total_tests = 0;
  uint passed_tests = 0;
  uint failed_tests = 0;
  std::vector<std::string> failed_test_names;
};

template <typename Function>
inline void test(Function&& function, const char* test_name, TestData& test_data) {
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










namespace CGAL {

template <typename Kernel_>
struct Grid_test_access {
  using TestGrid = CGAL::Box_grid<Kernel_>;

  using Kernel   = Kernel_;
  using FT       = typename Kernel::FT;
  using Point_3  = typename Kernel::Point_3;
  using Point_2  = typename Kernel::Point_2;
  using Vector_3 = typename Kernel::Vector_3;
  using Index    = std::size_t;
  using Cell = typename Box_grid<Kernel>::Cell;

  static TestGrid make_grid(FT box_size, FT min_x, FT min_y, FT min_z, FT max_x, FT max_y, FT max_z) {
    return TestGrid(box_size, min_x, min_y, min_z, max_x, max_y, max_z);
  }

  static TestGrid make_grid_from_points(FT box_size, const std::vector<Point_3>& points) {
    auto bbox = CGAL::bounding_box(points.begin(), points.end());
    FT min_x = bbox.xmin();
    FT min_y = bbox.ymin();
    FT min_z = bbox.zmin();
    FT max_x = bbox.xmax();
    FT max_y = bbox.ymax();
    FT max_z = bbox.zmax();
    return make_grid(box_size, min_x, min_y, min_z, max_x, max_y, max_z);
  }

  static std::vector<Vector_3> compute_block_normals(TestGrid& grid) {
    return grid.compute_block_normals();
  }

  static bool to_grid_coordinates(TestGrid& grid, const Point_3& p, int& ix, int& iy, int& iz) {
    return grid.to_grid_coordinates(p, ix, iy, iz);
  }
  
  static std::size_t flat_index(TestGrid& grid, int ix, int iy, int iz) {
    return grid.flat_index(ix, iy, iz);
  }

  static Point_3 cell_center(TestGrid& grid, int ix, int iy, int iz) {
    return grid.cell_center(ix, iy, iz);
  }
  
  static const Cell* cell(TestGrid& grid, int ix, int iy, int iz) {
    return grid.cell(ix, iy, iz);
  }

  static void insert(TestGrid& grid, Index point_id, const Point_3& p, const Vector_3& n) {
    grid.insert(point_id, p, n);
  }

  static const std::vector<Index>& point_ids_within_distance_d_of_box(TestGrid& grid, 
    std::size_t ix, std::size_t iy, std::size_t iz) {

    return grid.point_ids_within_distance_d_of_box(ix, iy, iz);
  }

  static const std::vector<Index>& point_ids_in_box(TestGrid& grid, 
    std::size_t ix, std::size_t iy, std::size_t iz) {

    return grid.point_ids_in_box(ix, iy, iz);
  }

  static const std::vector<Index> nearby_point_ids(TestGrid& grid, Point_3& query, FT radius) {
    return grid.nearby_point_ids(query, radius);
  }

  static std::size_t nx(TestGrid& grid) {
    return grid.nx_;
  }
  static std::size_t ny(TestGrid& grid) {
    return grid.ny_;
  }
  static std::size_t nz(TestGrid& grid) {
    return grid.nz_;
  }
};


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

using TestGrid = CGAL::Grid_test_access<Kernel>::TestGrid;
using GridTestAccess = CGAL::Grid_test_access<Kernel>;
using Cell = CGAL::Box_grid<Kernel>::Cell;
using Index = CGAL::Box_grid<Kernel>::Index;

inline Polyhedron ear_clip_and_add_faces_helper(std::vector<Point> points, std::vector<Vector_3> normals) {
  
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








