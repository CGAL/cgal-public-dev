#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/refine_mesh_at_isolevel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

// Define types
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
// using Kernel = CGAL::Simple_cartesian<double>;
using Point_3 = Kernel::Point_3;
using Triangle_mesh = CGAL::Surface_mesh<Point_3>;
using VertexDescriptor = boost::graph_traits<Triangle_mesh>::vertex_descriptor;
using FaceDescriptor = boost::graph_traits<Triangle_mesh>::face_descriptor;
using EdgeDescriptor = boost::graph_traits<Triangle_mesh>::edge_descriptor;
using HalfedgeDescriptor = boost::graph_traits<Triangle_mesh>::halfedge_descriptor;
using VertexDistanceMap = Triangle_mesh::Property_map<VertexDescriptor, double>;

using HeatMethodDirect = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh,
    CGAL::Heat_method_3::Direct>;
using HeatMethodIDT = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh,
    CGAL::Heat_method_3::Intrinsic_Delaunay>;
using HeatMethodIM = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh,
    CGAL::Heat_method_3::Direct,
    CGAL::Heat_method_3::mollification_scheme_constant>;
using Traits = CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh>;
using Surface_mesh_shortest_path = CGAL::Surface_mesh_shortest_path<Traits>;
typedef typename Traits::Vector_3 Vector_3;

// Function to calculate the mean of a vector
template <typename T>
T calculate_mean(const std::vector<T> &data)
{
  if (data.empty())
    return T(0);
  return std::accumulate(data.begin(), data.end(), T(0)) / data.size();
}

// Function to calculate the standard deviation of a vector
template <typename T>
T calculate_standard_deviation(const std::vector<T> &data)
{
  if (data.size() < 2)
    return T(0);
  T mean_value = calculate_mean(data);
  T squared_sum = std::accumulate(data.begin(), data.end(), T(0), [mean_value](T sum, T value) {
    T diff = value - mean_value;
    return sum + diff * diff;
  });
  return std::sqrt(squared_sum / (data.size() - 1));
}

// Function to copy a mesh including property maps
template <typename MeshType>
MeshType copy_mesh(MeshType &original_mesh)
{
  // Create a new mesh to store the copy
  MeshType copied_mesh;

  // Use an unordered_map to keep track of elements.
  std::unordered_map<VertexDescriptor, VertexDescriptor> v2v;
  std::unordered_map<HalfedgeDescriptor, HalfedgeDescriptor> h2h;
  std::unordered_map<FaceDescriptor, FaceDescriptor> f2f;
  // Copy the face graph from the original mesh to the new mesh
  CGAL::copy_face_graph(original_mesh,
      copied_mesh,
      CGAL::parameters::vertex_to_vertex_output_iterator(std::inserter(v2v, v2v.end()))
          .halfedge_to_halfedge_output_iterator(std::inserter(h2h, h2h.end()))
          .face_to_face_output_iterator(std::inserter(f2f, f2f.end())));

  // Copy vertex properties
  auto original_vertex_point_map = get(CGAL::vertex_point, original_mesh);
  auto copied_vertex_point_map = get(CGAL::vertex_point, copied_mesh);
  for (const auto &v : vertices(original_mesh)) {
    put(copied_vertex_point_map, v2v[v], get(original_vertex_point_map, v));
  }

  // Copy vertex properties with key type double
  auto original_v_property_maps = original_mesh.template properties<VertexDescriptor>();
  for (const std::string &prop_name : original_v_property_maps) {
    if (original_mesh.template property_type<VertexDescriptor>(prop_name) != typeid(double))
      continue;
    auto original_v_map = *original_mesh.template property_map<VertexDescriptor, double>(prop_name);
    auto copied_v_map = copied_mesh
                            .template add_property_map<VertexDescriptor, double>(
                                prop_name, std::numeric_limits<double>::infinity())
                            .first;

    for (const auto &v : vertices(original_mesh)) {
      put(copied_v_map, v2v[v], get(original_v_map, v));
    }
  }

  return copied_mesh;
}

// Function to calculate relative error
double calculate_relative_error(double estimated_value, double true_value) {
    double absolute_error = std::abs(estimated_value - true_value);
    // FIXME: handle zero true_value
    double relative_error =
        true_value != 0 ? absolute_error / std::abs(true_value) : absolute_error;
    return relative_error;
}

int main(int argc, char *argv[])
{
  // Input filename
  const std::string filename =
      (argc > 1) ? argv[1]
                 : CGAL::data_file_path(
                       "/home/islam/code/cgal-public-dev/Data/data/meshes/larger_sphere.off");

  Triangle_mesh mesh;
  if (!CGAL::IO::read_polygon_mesh(filename, mesh) || CGAL::is_empty(mesh) ||
      !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Error: Invalid or empty input file or not a triangle mesh." << std::endl;
    return EXIT_FAILURE;
  }

  // Property maps for distance values
  auto [vertex_distance_direct, vertex_distance_idt, vertex_distance_im, vertex_distance_true] =
      std::tuple{mesh.add_property_map<VertexDescriptor, double>(
                         "v:distance_direct", std::numeric_limits<double>::infinity())
                     .first,
          mesh.add_property_map<VertexDescriptor, double>(
                  "v:distance_idt", std::numeric_limits<double>::infinity())
              .first,
          mesh.add_property_map<VertexDescriptor, double>(
                  "v:distance_im", std::numeric_limits<double>::infinity())
              .first,
          mesh.add_property_map<VertexDescriptor, double>(
                  "v:distance_true", std::numeric_limits<double>::infinity())
              .first};

  // Initialize heat methods
  VertexDescriptor source_vertex = *vertices(mesh).first;
  HeatMethodDirect heat_method_direct(mesh);
  heat_method_direct.add_source(source_vertex);
  heat_method_direct.estimate_geodesic_distances(vertex_distance_direct);

  HeatMethodIDT heat_method_idt(mesh);
  heat_method_idt.add_source(source_vertex);
  heat_method_idt.estimate_geodesic_distances(vertex_distance_idt);

  HeatMethodIM heat_method_im(mesh);
  heat_method_im.add_source(source_vertex);
  heat_method_im.estimate_geodesic_distances(vertex_distance_im);

  // gound truth
  Surface_mesh_shortest_path shortest_paths(mesh);
  bool has_degenerate_faces = CGAL::Heat_method_3::internal::has_degenerate_faces(mesh, Traits());
  if (!has_degenerate_faces) {
    Surface_mesh_shortest_path shortest_paths(mesh);
    shortest_paths.add_source_point(source_vertex);

    for (auto v : vertices(mesh)) {
      double true_distance = shortest_paths.shortest_distance_to_source_points(v).first;
      put(vertex_distance_true, v, true_distance);
    }
  }
  else {
    vertex_distance_true = vertex_distance_idt;
  }

  // Compute errors
  size_t vertices_cnt = num_vertices(mesh);
  std::vector<double> dir_errors(vertices_cnt);
  std::vector<double> im_errors(vertices_cnt);
  std::vector<double> idt_errors(vertices_cnt);

  // Compute min and max distance
  double min_distance = std::numeric_limits<double>::max();
  double max_distance = std::numeric_limits<double>::lowest();
  for (const auto vertex : vertices(mesh)) {
    size_t index = vertex.idx();
    dir_errors[index] = calculate_relative_error(get(vertex_distance_direct, vertex), get(vertex_distance_true, vertex));
    idt_errors[index] = calculate_relative_error(get(vertex_distance_idt, vertex), get(vertex_distance_true, vertex));
    im_errors[index] = calculate_relative_error(get(vertex_distance_im, vertex), get(vertex_distance_true, vertex));
    double dist1 = get(vertex_distance_direct, vertex);
    double dist2 = get(vertex_distance_idt, vertex);
    double dist3 = get(vertex_distance_im, vertex);
    min_distance = std::min(std::min(min_distance, dist1), std::min(dist2, dist3));
    max_distance = std::max(std::max(max_distance, dist1), std::max(dist2, dist3));
  }

  // Calculate and display mean and standard deviation
  double mean_direct_error = calculate_mean(dir_errors);
  double stddev_direct_error = calculate_standard_deviation(dir_errors);

  double mean_idt_error = calculate_mean(idt_errors);
  double stddev_idt_error = calculate_standard_deviation(idt_errors);

  double mean_im_error = calculate_mean(im_errors);
  double stddev_im_error = calculate_standard_deviation(im_errors);

  std::cout << "Direct Method Errors: Mean = " << mean_direct_error
            << ", Std Dev = " << stddev_direct_error << std::endl;
  std::cout << "IDT Method Errors: Mean = " << mean_idt_error << ", Std Dev = " << stddev_idt_error
            << std::endl;
  std::cout << "IM Method Errors: Mean = " << mean_im_error << ", Std Dev = " << stddev_im_error
            << std::endl;

  std::cout << "Distances: Max = " << max_distance << ", min = " << min_distance << std::endl;
  // Property map for colors
  typedef CGAL::Color Color;
  // auto vertex_color_map =
  //     mesh.add_property_map<VertexDescriptor, Color>("v:color", Color(255, 255, 255)).first;
  auto ecm = mesh.add_property_map<EdgeDescriptor, bool>("e:is_constrained", 0).first;
  // default isovalues for cutting the mesh
  std::vector<double> isovalues = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  for (auto &val : isovalues) {
    val *= (max_distance - min_distance) / isovalues.size();
  }

  // Create a map of distance names to filenames
  std::map<std::string, std::string> distance_to_filename = {
      {"vertex_distance_direct", "output_direct.ply"},
      {"vertex_distance_idt", "output_idt.ply"},
      {"vertex_distance_im", "output_im.ply"},
      {"vertex_distance_true", "output_true.ply"}};

  std::array<std::pair<std::string, std::string>, 4> distances_lst = {
      std::make_pair("v:distance_direct", "vertex_distance_direct"),
      std::make_pair("v:distance_idt", "vertex_distance_idt"),
      std::make_pair("v:distance_im", "vertex_distance_im"),
      std::make_pair("v:distance_true", "vertex_distance_true")};
  for (std::size_t i = 0, n = has_degenerate_faces ? 3 : 4; i < n; ++i) {
    auto tm = copy_mesh(mesh);
    auto ecm = tm.add_property_map<EdgeDescriptor, bool>("e:is_constrained", 0).first;

    auto vertex_color_map =
        tm.add_property_map<VertexDescriptor, Color>("v:color", Color(255, 255, 255)).first;

    // auto vertex_distance = dis.first;
    auto vertex_distance = *tm.template property_map<VertexDescriptor, double>(distances_lst[i].first);
    const std::string &distance_name = distances_lst[i].second;

    // refine the tm along isovalues
    for (double isovalue : isovalues)
      CGAL::Polygon_mesh_processing::refine_mesh_at_isolevel(
          tm, vertex_distance, isovalue, CGAL::parameters::edge_is_constrained_map(ecm));

    // Initialize a new distance property map for the refined tm
    auto new_vertex_distance = tm.add_property_map<VertexDescriptor, double>(
                                     "v:distance", std::numeric_limits<double>::infinity())
                                   .first;

    // Interpolate distance values for new vertices
    for (FaceDescriptor fd : faces(tm)) {
      std::vector<VertexDescriptor> vertices_of_face;
      for (auto v : vertices_around_face(halfedge(fd, tm), tm)) {
        vertices_of_face.push_back(v);
      }

      // Assign distance to new vertices based on average of surrounding vertices
      for (VertexDescriptor v : vertices_of_face) {
        double distance = 0;
        if (get(vertex_distance, v) != std::numeric_limits<double>::infinity()) {
          distance = get(vertex_distance, v);
        }
        else {
          // Interpolation: average of neighboring vertex distances
          std::vector<double> distances;
          for (VertexDescriptor neighbor : vertices_around_face(halfedge(v, tm), tm)) {
            if (get(vertex_distance, neighbor) != std::numeric_limits<double>::infinity()) {
              distances.push_back(get(vertex_distance, neighbor));
            }
          }
          if (!distances.empty()) {
            distance = std::accumulate(distances.begin(), distances.end(), 0.0) / distances.size();
          }
        }
        put(new_vertex_distance, v, distance);
      }
    }

    // Map distances to colors
    for (VertexDescriptor vd : vertices(tm)) {
      double dist = get(new_vertex_distance, vd);
      double normalized_dist = (dist - min_distance) / (max_distance - min_distance);
      unsigned char r = static_cast<unsigned char>(255 * normalized_dist);
      unsigned char g = 0;
      unsigned char b = static_cast<unsigned char>(
          255 * (1 - normalized_dist));  // Simple gradient from red to blue
      put(vertex_color_map, vd, Color(r, g, b));
    }

    // Save the tm with distance information to a file
    auto output_filename = distance_to_filename.at(distance_name);
    if (!CGAL::IO::write_polygon_mesh(output_filename, tm, CGAL::parameters::vertex_color_map(vertex_color_map))) {
      std::cerr << "Failed to write output file: " << output_filename << std::endl;
      return EXIT_FAILURE;
    }

    // // Iterate over all edges in the tm
    // std::vector<EdgeDescriptor> edges_to_remove;
    // for (const auto &e : edges(tm)) {
    //   // Check if the edge was marked as a constrained edge
    //   if (get(ecm, e)) {
    //     edges_to_remove.push_back(e);
    //   }
    // }

    // // Remove the marked edges
    // for (const auto &e : edges_to_remove) {
    //   remove_edge(e, tm);
    // }

    // Collect vertices
    std::vector<VertexDescriptor> vs;
    for (const auto &v : vertices(tm)) {
      vs.push_back(v);
    }

    // Collect edges
    std::vector<std::pair<int, int>> es;
    for (const auto &e : edges(tm)) {
      auto src = source(e, tm);
      auto tr = target(e, tm);
      if (get(ecm, e)) {
        es.emplace_back(src, tr);
      }
    }

    // Open file for writing
    std::ofstream file(std::string("edges_") + output_filename);
    if (!file) {
      std::cerr << "Error: Cannot open file " << filename << " for writing." << std::endl;
      return -1;
    }

    // Write PLY header
    file << "ply\n";
    file << "format ascii 1.0\n";
    file << "element vertex " << vs.size() << "\n";
    file << "property float x\n";
    file << "property float y\n";
    file << "property float z\n";
    file << "element edge " << es.size() << "\n";
    file << "property int vertex1\n";
    file << "property int vertex2\n";
    file << "end_header\n";

    // Write vertices (assuming CGAL's default point type for the mesh)
    auto point_map = get(CGAL::vertex_point, tm);
    for (const auto &v : vs) {
      const auto &point = get(point_map, v);
      file << point.x() << " " << point.y() << " " << point.z() << "\n";
    }

    // Write edges
    for (const auto &edge : es) {
      file << edge.first << " " << edge.second << "\n";
    }

    file.close();
  }

  return EXIT_SUCCESS;
}
