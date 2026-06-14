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

bool read_off_file(const std::string& filename,
                   std::vector<Point>& points,
                   std::vector<Vector_3>& normals)
{
  points.clear();
  normals.clear();

  CGAL::Surface_mesh<Point> mesh;

  if (!CGAL::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Error: cannot read OFF file " << filename << std::endl;
    return false;
  }

  auto normal_map_opt =
    mesh.property_map<CGAL::Surface_mesh<Point>::Vertex_index, Vector_3>("v:normal");

  const bool has_normals = normal_map_opt.has_value();

  for (auto vd : mesh.vertices()) {
    points.push_back(mesh.point(vd));

    if (has_normals) {
      normals.push_back((*normal_map_opt)[vd]);
    } else {
      normals.push_back(CGAL::NULL_VECTOR);
    }
  }

  return true;
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

  const double max_extent = (std::max)({extent_x, extent_y, extent_z});
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

template <typename Mesh>
bool write_mesh_graph_ply(const Mesh& mesh, const std::string& filename)
{
  using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;

  std::ofstream out(filename);
  if (!out) {
    std::cerr << "Error: cannot open " << filename << " for writing.\n";
    return false;
  }

  std::unordered_map<vertex_descriptor, int> vindex;
  int idx = 0;

  for (auto vd : vertices(mesh)) {
    vindex[vd] = idx++;
  }

  std::vector<std::pair<int,int>> edges_list;
  for (auto e : edges(mesh)) {
    halfedge_descriptor h = halfedge(e, mesh);
    vertex_descriptor s = source(h, mesh);
    vertex_descriptor t = target(h, mesh);

    int is = vindex[s];
    int it = vindex[t];
    if (is > it) std::swap(is, it);
    edges_list.emplace_back(is, it);
  }

  out << "ply\n";
  out << "format ascii 1.0\n";
  out << "element vertex " << num_vertices(mesh) << "\n";
  out << "property float x\n";
  out << "property float y\n";
  out << "property float z\n";
  out << "element edge " << edges_list.size() << "\n";
  out << "property int vertex1\n";
  out << "property int vertex2\n";
  out << "end_header\n";

  for (auto vd : vertices(mesh)) {
    const auto& p = get(CGAL::vertex_point, mesh, vd);
    out << CGAL::to_double(p.x()) << " "
        << CGAL::to_double(p.y()) << " "
        << CGAL::to_double(p.z()) << "\n";
  }

  for (const auto& e : edges_list) {
    out << e.first << " " << e.second << "\n";
  }

  return true;
}

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

  // const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/kitten.xyz");
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/plane.off");

  std::vector<Point> points;
  std::vector<Vector_3> normals;
  std::cout << "Filename: " << filename << std::endl;

  // if extension is .off
  if (std::filesystem::path(filename).extension() == ".off") {
    if (!read_off_file(filename, points, normals)) {
      std::cerr << "Error: cannot read input file!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else if (std::filesystem::path(filename).extension() == ".xyz") {
    if (!read_xyz_file(filename, points, normals)) {
      std::cerr << "Error: cannot read input file!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else {
    std::cerr << "Error: unsupported file format!" << std::endl;
    return EXIT_FAILURE;
  }

  // Normalize geometry first.
  std::cout<<"Centering and scaling point cloud between [-1,1]^3 ..." << std::endl;
  center_and_scale_point_cloud(points);

  // Then compute normals if the file did not provide them.
  compute_normals_if_missing(points, normals, 6);

  // Average spacing should be computed in the normalized coordinate system.
  double average_spacing = 0.5*CGAL::compute_average_spacing<CGAL::Parallel_if_available_tag>(points, 6);

  Polyhedron output_mesh;
  const auto bbox = CGAL::bounding_box(points.begin(), points.end()); // recompute bbox after centering and scaling

  // Build the grid and insert points + normals.
  CGAL::Box_grid<Kernel> grid{FT(average_spacing),
                              FT(bbox.xmin()-1e-3), FT(bbox.xmax()+1e-3),
                              FT(bbox.ymin()-1e-3), FT(bbox.ymax()+1e-3),
                              FT(bbox.zmin()-1e-3), FT(bbox.zmax()+1e-3)}; // initialize grid with cell size equal to average spacing and bounding box [-1,1]^3
  grid.build(points, normals); // insert points and normals into the grid

  std::vector<Vector_3> block_normals = grid.compute_block_normals(); // compute block normals by averaging point normals in each cell
  std::cout<<"Computed " << block_normals.size() << " block normals." << std::endl;

  std::vector<FT> splat_sizes = grid.estimate_individual_splat_sizes(); // estimate individual splat sizes based on local point distribution
  std::cout<<"Estimated individual splat sizes for " << splat_sizes.size() << " points." << std::endl;

  grid.write_point_cloud_ply("debug_points.ply");
  grid.write_grid_vertices_ply("debug_grid_vertices.ply");
  grid.write_cell_centers_and_normals_ply("debug_cell_normals.ply", 0.2);


  CGAL::Splat_surface_reconstruction_3<std::vector<Point>, std::vector<Vector_3>, Polyhedron> reconstruction(grid, output_mesh);
  reconstruction.run();

  write_mesh_graph_ply(output_mesh, "graph.ply");
  std::cout << "Mesh graph written to graph.ply" << std::endl;

  CGAL::IO::write_polygon_mesh("mesh.ply", output_mesh);
  std::cout << "Reconstructed mesh written to mesh.ply" << std::endl;

  // Splat_surface_reconstruction(points, normals, output_mesh, average_spacing);

  // if (CGAL::splat_surface_reconstruction(points, normals, output_mesh, average_spacing)) {
  //   std::string fname = std::filesystem::path(filename).stem().string() + ".off";
  //   CGAL::IO::write_polygon_mesh(fname, output_mesh);
  // }

  // else {
  //   return EXIT_FAILURE;
  // }

  return EXIT_SUCCESS;
}