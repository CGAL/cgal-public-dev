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

double center_and_scale_point_cloud(std::vector<Point>& points) {
  if (points.empty()) {
    return 1.0;
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
    return 1.0;
  }

  // Scale so the largest side fits exactly in [-1,1].
  const double scale = 2.0 / max_extent;

  for (Point& p : points) {
    const double x = (CGAL::to_double(p.x()) - center_x) * scale;
    const double y = (CGAL::to_double(p.y()) - center_y) * scale;
    const double z = (CGAL::to_double(p.z()) - center_z) * scale;
    p = Point(x, y, z);
  }
  return scale;
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

    for (std::size_t i = 0; i < normals.size(); ++i) {
      CGAL_assertion(normals[i] != CGAL::NULL_VECTOR);
      const double len2 = CGAL::to_double(normals[i].squared_length());
      CGAL_assertion(len2 > 0.0);
      normals[i] = normals[i] / std::sqrt(len2);
    }
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

int main(int argc, char* argv[]) {
  const std::string filename = (argc>1)?argv[1]:CGAL::data_file_path("meshes/plane.off");
  // const std::string filename = CGAL::data_file_path("meshes/plane.off");
  // const std::string filename = (argc>1)?argv[1]:CGAL::data_file_path("points_3/fold.xyz");
  // const std::string filename = CGAL::data_file_path("meshes/fold.off");

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
  double sc = center_and_scale_point_cloud(points);

  // Then compute normals if the file did not provide them.
  compute_normals_if_missing(points, normals, 6);

  // Average spacing should be computed in the normalized coordinate system.
  double scale = (argc > 2) ? std::stod(argv[2]) : 1.0;
  double average_spacing = scale * CGAL::compute_average_spacing<CGAL::Parallel_if_available_tag>(points, 2);

  // double average_spacing = sc * std::stod(argv[2]);

  Polyhedron output_mesh;
  const auto bbox = CGAL::bounding_box(points.begin(), points.end()); // recompute bbox after centering and scaling

  const FT padding = FT(1.1) * average_spacing;
  // const FT padding = 1e-8;
  // Build the grid and insert points + normals.
  CGAL::Box_grid<Kernel> grid{FT(average_spacing),
                              FT(bbox.xmin()-padding), FT(bbox.xmax()+padding),
                              FT(bbox.ymin()-padding), FT(bbox.ymax()+padding),
                              FT(bbox.zmin()-padding), FT(bbox.zmax()+padding)}; // initialize grid with cell size equal to average spacing and bounding box [-1,1]^3
  grid.build(points, normals); // insert points and normals into the grid

  std::vector<Vector_3> block_normals = grid.compute_block_normals(); // compute block normals by averaging point normals in each cell
  std::cout<<"Computed " << block_normals.size() << " block normals." << std::endl;

  std::vector<FT> splat_sizes = grid.estimate_individual_splat_sizes(); // estimate individual splat sizes based on local point distribution
  std::cout<<"Estimated individual splat sizes for " << splat_sizes.size() << " points." << std::endl;

  grid.fill_empty_block_normals_from_large_splats();

  grid.write_point_cloud_ply("debug_points.ply");
  grid.write_grid_vertices_ply("debug_grid_vertices.ply");
  grid.write_cell_centers_and_normals_ply("debug_cell_normals.ply", 0.2);

  CGAL::Splat_surface_reconstruction_3<std::vector<Point>, std::vector<Vector_3>, Polyhedron> reconstruction(grid, output_mesh);
  reconstruction.run();
  
  CGAL::IO::write_polygon_mesh("mesh.ply", output_mesh);
  std::cout << "Reconstructed mesh written to mesh.ply" << std::endl;

  // std::ofstream out("rejected_candidates.ply");
  // if (!out) {
  //   std::cerr << "Error: cannot open output file!" << std::endl;
  //   return EXIT_FAILURE;
  // }
  // out << "ply\n";
  // out << "format ascii 1.0\n";
  // out << "element vertex " << reconstruction.rejected_candidates_.size() << "\n";
  // out << "property float x\n";
  // out << "property float y\n";
  // out << "property float z\n";
  // // out << "element edge " << 2*reconstruction.rejected_candidates_.size() << "\n";
  // // out << "property int vertex1\n";
  // // out << "property int vertex2\n";
  // out << "end_header\n";

  // // for(const auto& vd : vertices(output_mesh))
  // // {
  // //   const auto& p = get(CGAL::vertex_point, output_mesh, vd);
  // //     out << p.x() << " "
  // //         << p.y() << " "
  // //         << p.z() << "\n";
  // // }
  // for (const auto& cand : reconstruction.rejected_candidates_)
  // {
  //   const auto& p = cand.p;
  //     out << p.x() << " "
  //         << p.y() << " "
  //         << p.z() << "\n";
  // }
  // // Write edges connecting rejected candidates to their parents
  // // for(int i = 0; i < reconstruction.rejected_candidates_.size(); ++i)
  // // {
  // //   // omit the v in the output, just write the index of the vertex
  // //   out << reconstruction.rejected_candidates_[i].parent0 << " " << i + vertices(output_mesh).size() << "\n";
  // //   out << reconstruction.rejected_candidates_[i].parent1 << " " << i + vertices(output_mesh).size() << "\n";
  // //   // out << reconstruction.rejected_candidates_[i].parent1 << " " << reconstruction.rejected_candidates_[i].parent0 << "\n";
  // // }

  return EXIT_SUCCESS;
}