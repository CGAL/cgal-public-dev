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

// ------------------------------------------------------------
// Read points + normals from file
// ------------------------------------------------------------

bool read_xyz_file(const std::string& filename,
                   std::vector<Point>& points,
                   std::vector<Vector_3>& normals) {

  points.clear();
  normals.clear();

  struct Appender {
    std::vector<Point>& points;
    std::vector<Vector_3>& normals;

    Appender(std::vector<Point>& points, std::vector<Vector_3>& normals)
      : points(points), normals(normals) {}

    void operator()(const Pwn& pwn) {
      points.push_back(pwn.first);
      normals.push_back(pwn.second);
    }
  };

  return CGAL::IO::read_points<Pwn>(
    filename,
    boost::make_function_output_iterator(Appender(points, normals)),
    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()).normal_map(CGAL::Second_of_pair_property_map<Pwn>())
  );
}

// ------------------------------------------------------------
// Center and scale into [-1,1]^3
// ------------------------------------------------------------

void center_and_scale_point_cloud(std::vector<Point>& points) {
  if (points.empty()) {
    return;
  }

  const auto bbox = CGAL::bounding_box(points.begin(), points.end());

  const double center_x = 0.5 * (bbox.xmin() + bbox.xmax());
  const double center_y = 0.5 * (bbox.ymin() + bbox.ymax());
  const double center_z = 0.5 * (bbox.zmin() + bbox.zmax());

  const double extent_x = bbox.xmax() - bbox.xmin();
  const double extent_y = bbox.ymax() - bbox.ymin();
  const double extent_z = bbox.zmax() - bbox.zmin();

  const double max_extent = std::max({extent_x, extent_y, extent_z});
  if (max_extent <= 0.0) {
    return;
  }

  // Scale so the largest side fits exactly in [-1,1].
  const double scale = 2.0 / max_extent;

  for (Point& p : points) {
    const double x = (CGAL::to_double(p.x()) - center_x) * scale;
    const double y = (CGAL::to_double(p.y()) - center_y) * scale;
    const double z = (CGAL::to_double(p.z()) - center_z) * scale;
    p = Point(x, y, z);
  }
}

// ------------------------------------------------------------
// Compute normals if needed
// ------------------------------------------------------------

void compute_normals_if_missing(std::vector<Point>& points,
                                std::vector<Vector_3>& normals,
                                std::size_t k = 6) {
  if (points.empty()) {
    return;
  }

  const bool has_normals = !normals.empty() && normals[0] != CGAL::NULL_VECTOR;

  if (has_normals) {
    std::cout<<"Input file already contains normals, skipping normal estimation." << std::endl;
    return;
  }

  std::cout << "Estimating normals..." << std::endl;
  normals.assign(points.size(), CGAL::NULL_VECTOR);

  std::vector<std::size_t> indices(points.size());
  std::iota(indices.begin(), indices.end(), 0);

  CGAL::jet_estimate_normals<CGAL::Sequential_tag>(
    indices, k,
    CGAL::parameters::point_map(CGAL::make_random_access_property_map(points)).normal_map(CGAL::make_random_access_property_map(normals)));

  CGAL::mst_orient_normals(
    indices, k,
    CGAL::parameters::point_map(CGAL::make_random_access_property_map(points)).normal_map(CGAL::make_random_access_property_map(normals)));
}

int main(int argc, char* argv[]) {
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/kitten.xyz");

  std::vector<Point> points;
  std::vector<Vector_3> normals;
  std::cout << "Filename: " << filename << std::endl;

  if (!read_xyz_file(filename, points, normals)) {
    std::cerr << "Error: cannot read input file!" << std::endl;
    return EXIT_FAILURE;
  }

  // Normalize geometry first.
  std::cout<<"Centering and scaling point cloud between [-1,1]^3 ..." << std::endl;
  center_and_scale_point_cloud(points);

  // Then compute normals if the file did not provide them.
  compute_normals_if_missing(points, normals, 6);

  // Average spacing should be computed in the normalized coordinate system.
  double average_spacing = CGAL::compute_average_spacing<CGAL::Parallel_if_available_tag>(points, 6);

  Polyhedron output_mesh;
  const auto bbox = CGAL::bounding_box(points.begin(), points.end()); // recompute bbox after centering and scaling

  // Build the grid and insert points + normals.
  CGAL::Box_grid<Kernel> grid{FT(average_spacing), 
                              FT(bbox.xmin()), FT(bbox.xmax()), 
                              FT(bbox.ymin()), FT(bbox.ymax()), 
                              FT(bbox.zmin()), FT(bbox.zmax())}; // initialize grid with cell size equal to average spacing and bounding box [-1,1]^3
  grid.build(points, normals); // insert points and normals into the grid
  
  std::vector<Vector_3> block_normals = grid.compute_block_normals(); // compute block normals by averaging point normals in each cell
  std::cout<<"Computed " << block_normals.size() << " block normals." << std::endl;

  // if (CGAL::splat_surface_reconstruction(points, normals, output_mesh, average_spacing)) {
  //   std::string fname = std::filesystem::path(filename).stem().string() + ".off";
  //   CGAL::IO::write_polygon_mesh(fname, output_mesh);
  // }

  // else {
  //   return EXIT_FAILURE;
  // }

  return EXIT_SUCCESS;
}