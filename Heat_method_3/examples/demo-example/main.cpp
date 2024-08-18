#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include "CGAL/boost/graph/graph_traits_Surface_mesh.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <limits>


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

// Function to calculate the mean of a vector
template <typename T>
T calculate_mean(const std::vector<T>& data) {
    if (data.empty()) return T(0);
    return std::accumulate(data.begin(), data.end(), T(0)) / data.size();
}

// Function to calculate the standard deviation of a vector
template <typename T>
T calculate_standard_deviation(const std::vector<T>& data) {
    if (data.size() < 2) return T(0);
    T mean_value = calculate_mean(data);
    T squared_sum = std::accumulate(data.begin(), data.end(), T(0), 
        [mean_value](T sum, T value) {
            T diff = value - mean_value;
            return sum + diff * diff;
        });
    return std::sqrt(squared_sum / (data.size() - 1));
}

int main(int argc, char* argv[]) {
    // Input filename
    const std::string filename = (argc > 1) ? argv[1] 
                                            : CGAL::data_file_path("/home/islam/code/cgal-public-dev/Data/data/meshes/larger_sphere.off");
    
    Triangle_mesh mesh;
    if (!CGAL::IO::read_polygon_mesh(filename, mesh) || CGAL::is_empty(mesh) || !CGAL::is_triangle_mesh(mesh)) {
        std::cerr << "Error: Invalid or empty input file or not a triangle mesh." << std::endl;
        return EXIT_FAILURE;
    }

    // Property maps for distance values
    auto [vertex_distance_direct, vertex_distance_idt, vertex_distance_im, vertex_distance_true] = 
        std::tuple{
            mesh.add_property_map<VertexDescriptor, double>("v:distance_direct", 0).first,
            mesh.add_property_map<VertexDescriptor, double>("v:distance_idt", 0).first,
            mesh.add_property_map<VertexDescriptor, double>("v:distance_im", 0).first,
            mesh.add_property_map<VertexDescriptor, double>("v:distance_true", 0).first
        };

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
    shortest_paths.add_source_point(source_vertex);

    for (auto v : vertices(mesh)) {
      double true_distance = shortest_paths.shortest_distance_to_source_points(v).first;
      put(vertex_distance_true, v, true_distance);
    }
    // Compute errors
    size_t vertices_cnt = num_vertices(mesh);
    std::vector<double> dir_errors(vertices_cnt);
    std::vector<double> im_errors(vertices_cnt);
    std::vector<double> idt_errors(vertices_cnt);

    for (const auto vertex : vertices(mesh)) {
        size_t index = vertex.idx();
        dir_errors[index] = get(vertex_distance_direct, vertex) - get(vertex_distance_true, vertex);
        idt_errors[index] = get(vertex_distance_idt, vertex) - get(vertex_distance_true, vertex);
        im_errors[index] = get(vertex_distance_im, vertex) - get(vertex_distance_true, vertex);
    }

    // Calculate and display mean and standard deviation
    double mean_direct_error = calculate_mean(dir_errors);
    double stddev_direct_error = calculate_standard_deviation(dir_errors);

    double mean_idt_error = calculate_mean(idt_errors);
    double stddev_idt_error = calculate_standard_deviation(idt_errors);

    double mean_im_error = calculate_mean(im_errors);
    double stddev_im_error = calculate_standard_deviation(im_errors);

    std::cout << "Direct Method Errors: Mean = " << mean_direct_error << ", Std Dev = " << stddev_direct_error << std::endl;
    std::cout << "IDT Method Errors: Mean = " << mean_idt_error
              << ", Std Dev = " << stddev_idt_error << std::endl;
    std::cout << "IM Method Errors: Mean = " << mean_im_error << ", Std Dev = " << stddev_im_error << std::endl;

    return EXIT_SUCCESS;
}
