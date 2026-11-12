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
typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

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
                   std::vector<Vector_3>& normals) {
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

  if (!has_normals) {
    std::cerr << "Warning: OFF file " << filename << " does not contain vertex normals." << std::endl;
  }

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

int main(int argc, char* argv[])
{
  std::string filename;
  std::string input_type;
  std::string output_filename = "mesh.obj";
  double d = 0.05;
  double feature_angle_minimum = 30.0;

  using Reconstruction =
    CGAL::Splat_surface_reconstruction_3<
      std::vector<Point>,
      std::vector<Vector_3>,
      Polyhedron>;

  Reconstruction::Initialization_method init_method =
    Reconstruction::Initialization_method::FEATURE_AWARE;

  // Parse arguments.
  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];

    auto next_value = [&]() -> std::string
    {
      if (i + 1 >= argc)
      {
        std::cerr << "Missing value for " << arg << std::endl;
        std::exit(EXIT_FAILURE);
      }
      return argv[++i];
    };

    if (arg == "--input")
      filename = next_value();

    else if (arg == "--type")
      input_type = next_value();

    else if (arg == "--output")
      output_filename = next_value();

    else if (arg == "--d")
      d = std::stod(next_value());

    else if (arg == "--feature-angle")
      feature_angle_minimum = std::stod(next_value());

    else if (arg == "--init")
    {
      const std::string method = next_value();

      if (method == "feature_aware")
        init_method = Reconstruction::Initialization_method::FEATURE_AWARE;
      else if (method == "default")
        init_method = Reconstruction::Initialization_method::DEFAULT;
      else
      {
        std::cerr << "Unknown initialization method: "
                  << method << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      std::cerr << "Unknown argument: " << arg << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Validate arguments.
  if (filename.empty() ||
      (input_type != "mesh" && input_type != "point-cloud"))
  {
    std::cerr
      << "Usage: " << argv[0]
      << " --input <file> --type <mesh|point-cloud>"
      << " [--d <value>]"
      << " [--init <default|feature_aware>]"
      << " [--feature-angle <degrees>]"
      << " [--output <file>]"
      << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input: " << filename << std::endl;
  std::cout << "Input type: " << input_type << std::endl;

  Polyhedron output_mesh;
  double average_spacing = 0.0;

  if (input_type == "point-cloud")
  {
    // Read point cloud.
    std::vector<Point> points;
    std::vector<Vector_3> normals;

    const std::string extension =
      std::filesystem::path(filename).extension().string();

    if (extension == ".off")
    {
      if (!read_off_file(filename, points, normals))
      {
        std::cerr << "Error: cannot read point cloud!"
                  << std::endl;
        return EXIT_FAILURE;
      }
    }
    else if (extension == ".xyz")
    {
      if (!read_xyz_file(filename, points, normals))
      {
        std::cerr << "Error: cannot read point cloud!"
                  << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      std::cerr << "Unsupported point-cloud format: "
                << extension << std::endl;
      return EXIT_FAILURE;
    }

    // Normalize geometry.
    std::cout << "Centering and scaling point cloud..."
              << std::endl;

    center_and_scale_point_cloud(points);

    average_spacing =
      CGAL::compute_average_spacing<
        CGAL::Parallel_if_available_tag>(points, 2);

    // Compute normals if missing.
    compute_normals_if_missing(points, normals, 6);

    std::cout << "Average spacing: " << average_spacing << std::endl;
    std::cout << "d: " << d << std::endl;

    // Construct point-cloud reconstruction.
    Reconstruction reconstruction(
      points,
      normals,
      d,
      output_mesh,
      init_method,
      FT(feature_angle_minimum)
    );

    reconstruction.run();
  }
  else // input_type == "mesh"
  {
    // Read input mesh.
    Polyhedron input_mesh;

    if (!CGAL::IO::read_polygon_mesh(filename, input_mesh) ||
        input_mesh.is_empty() ||
        !CGAL::is_triangle_mesh(input_mesh))
    {
      std::cerr << "Error: cannot read a valid triangle mesh!"
                << std::endl;
      return EXIT_FAILURE;
    }

    // Collect input mesh vertices.
    std::vector<Point> points;
    points.reserve(num_vertices(input_mesh));

    for (vertex_descriptor v : vertices(input_mesh))
      points.push_back(input_mesh.point(v));

    // Normalize geometry using the same normalization as
    // the point-cloud implementation.
    std::cout << "Centering and scaling input mesh..."
              << std::endl;

    center_and_scale_point_cloud(points);

    // Apply the same transformation to the input mesh.
    auto vertex_it = vertices(input_mesh).begin();

    for (std::size_t i = 0; i < points.size(); ++i, ++vertex_it)
      input_mesh.point(*vertex_it) = points[i];

    // Estimate spacing from the mesh vertices.
    average_spacing =
      CGAL::compute_average_spacing<
        CGAL::Parallel_if_available_tag>(points, 2);

    std::cout << "Average spacing: " << average_spacing << std::endl;
    std::cout << "d: " << d << std::endl;

    // Construct mesh-input reconstruction.
    Reconstruction reconstruction(
      input_mesh,
      d,
      output_mesh,
      init_method,
      FT(feature_angle_minimum)
    );

    reconstruction.run();
  }

  // Write output.
  if (!CGAL::IO::write_polygon_mesh(output_filename, output_mesh))
  {
    std::cerr << "Error: could not write output mesh!"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Reconstructed mesh written to "
            << output_filename << std::endl;

  return EXIT_SUCCESS;
}
