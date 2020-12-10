#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_reconstruction_3.h>
#include <CGAL/IO/OFF_reader.h>

using EPICK   = CGAL::Exact_predicates_inexact_constructions_kernel;
using Kernel  = EPICK;
using Point_3 = typename Kernel::Point_3;
using KSR     = CGAL::Kinetic_shape_reconstruction_3<Kernel>;

struct Polygon_map {

  using key_type   = std::vector<std::size_t>;
  using value_type = std::vector<Point_3>;
  using reference  = value_type;
  using category   = boost::readable_property_map_tag;

  const std::vector<Point_3>& points;
  Polygon_map(
    const std::vector<Point_3>& vertices) :
  points(vertices)
  { }

  friend reference get(const Polygon_map& map, const key_type& face) {
    reference polygon;
    polygon.reserve(face.size());
    std::transform(
      face.begin(), face.end(),
      std::back_inserter(polygon),
      [&](const std::size_t vertex_index) -> Point_3 {
        return map.points[vertex_index];
      });
    return polygon;
  }
};

const bool run_test(
  const std::string input_filename,
  const std::size_t num_iters,
  std::size_t& num_tests) {

  ++num_tests;
  std::ifstream input_file(input_filename);
  std::vector<Point_3> input_vertices;
  std::vector< std::vector<std::size_t> > input_faces;
  assert(CGAL::read_OFF(input_file, input_vertices, input_faces));

  std::cout << std::endl;
  std::cout << "--INPUT FILE: " << input_filename << std::endl;
  const Polygon_map polygon_map(input_vertices);
  for (unsigned int k = 1; k <= 6; ++k) {
    std::cout << std::endl << "--INPUT K: " << k << std::endl;
    for (std::size_t iter = 0; iter < num_iters; ++iter) {
      std::cout << std::endl << "--ITERATION #" << iter + 1 << " BEGIN!" << std::endl;
      KSR ksr(false, false);
      assert(ksr.partition(input_faces, polygon_map, k));
      ksr.clear();
      std::cout << std::endl << "--ITERATION #" << iter + 1 << " END!" << std::endl;
    }
  }

  std::cout << std::endl << "--INPUT K: " << 100 << std::endl;
  for (std::size_t iter = 0; iter < num_iters; ++iter) {
    std::cout << std::endl << "--ITERATION #" << iter + 1 << " BEGIN!" << std::endl;
    KSR ksr(false, false);
    assert(ksr.partition(input_faces, polygon_map, 100));
    ksr.clear();
    std::cout << std::endl << "--ITERATION #" << iter + 1 << " END!" << std::endl;
  }

  return true;
}

int main (const int argc, const char** argv) {

  std::size_t num_tests = 0;
  const std::size_t num_iters = 3;

  // Stress tests 0.
  assert(run_test("data/stress-test-0/test-1-polygon-a.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-1-polygon-b.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-1-polygon-c.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-1-polygon-d.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-2-polygons-ab.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-2-polygons-ac.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-2-polygons-ad.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-2-polygons-bc.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-2-polygons-bd.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-2-polygons-cd.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-3-polygons-abc.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-3-polygons-abd.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-3-polygons-acd.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-3-polygons-bcd.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-4-polygons-abcd.off", num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-6-polygons.off", num_iters, num_tests));

  // Stress tests 1.
  assert(run_test("data/stress-test-1/test-1-rnd-polygons-1-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-2-rnd-polygons-1-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-3-rnd-polygons-1-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-4-rnd-polygons-1-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-5-rnd-polygons-2-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-6-rnd-polygons-2-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-7-rnd-polygons-2-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-8-rnd-polygons-3-4.off", num_iters, num_tests));

  // Stress tests 2.
  assert(run_test("data/stress-test-2/test-1-rnd-polygons-1-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-2/test-2-rnd-polygons-1-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-2/test-3-rnd-polygons-1-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-2/test-4-rnd-polygons-1-3.off", num_iters, num_tests));
  assert(run_test("data/stress-test-2/test-5-rnd-polygons-2-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-2/test-6-rnd-polygons-3-4.off", num_iters, num_tests));

  // Stress tests 3.
  assert(run_test("data/stress-test-3/test-1-rnd-polygons-2-3.off", num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-2-rnd-polygons-2-3.off", num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-3-rnd-polygons-2-3.off", num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-4-rnd-polygons-2-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-5-rnd-polygons-1-3.off", num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-6-rnd-polygons-2-3.off", num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-7-rnd-polygons-2-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-8-rnd-polygons-2-10.off", num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-9-rnd-polygons-4-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-10-rnd-polygons-5-4.off", num_iters, num_tests));

  // Stress tests 4.
  assert(run_test("data/stress-test-4/test-1-rnd-polygons-2-6.off", num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-2-rnd-polygons-3-8.off", num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-3-rnd-polygons-4-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-4-rnd-polygons-4-6.off", num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-5-rnd-polygons-6-4.off", num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-6-rnd-polygons-5-6.off", num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-7-rnd-polygons-7-6.off", num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-8-rnd-polygons-7-8.off", num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-9-rnd-polygons-12-4.off", num_iters, num_tests));

  // Real data tests.
  assert(run_test("data/real-data-test/building_b_15squares_15planes.off", num_iters, num_tests));

  std::cout << std::endl << "--OUTPUT STATS:" << std::endl;
  std::cout << "* number of iterations per test: " << num_iters << std::endl;
  std::cout << "* k intersections: [1, 6]" << std::endl;

  std::cout << std::endl << "ALL " << num_tests << " TESTS SUCCESS!" << std::endl;
  return EXIT_SUCCESS;
}
