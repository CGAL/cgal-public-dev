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

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector_3;
typedef std::pair<Point, Vector_3> Pwn;
typedef CGAL::Surface_mesh<Point> Polyhedron;
typedef Kernel::FT FT;

typedef CGAL::Splat_surface_reconstruction_3<std::vector<Point>, std::vector<Vector_3>, Polyhedron> Splat_surface_reconstruction;

bool test_circle_splat_intersection_exact()
{
  using Grid = CGAL::Box_grid<Kernel>;

  // Create a dummy grid.
  Grid grid(
      FT(1.0),
      FT(-2.0), FT(2.0),
      FT(-2.0), FT(2.0),
      FT(-2.0), FT(2.0));

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

int main(int argc, char* argv[]) {
  test_circle_splat_intersection_exact();


  return EXIT_SUCCESS;
}