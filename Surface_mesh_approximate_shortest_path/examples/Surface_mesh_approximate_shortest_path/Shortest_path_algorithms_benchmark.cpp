//#ifndef CGAL_EIGEN3_ENABLED
//#define CGAL_EIGEN3_ENABLED
//#endif

//#include <Eigen/Dense.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximate_shortest_path.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <iostream>
#include <chrono>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef Kernel::Point_3 Point_3;

typedef boost::graph_traits<Surface_mesh> Graph_traits;

typedef Graph_traits::vertex_descriptor vertex_descriptor;
typedef Graph_traits::edge_descriptor edge_descriptor;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;
typedef Graph_traits::face_descriptor face_descriptor;

typedef CGAL::Surface_mesh_approximate_shortest_path_traits<Kernel, Surface_mesh>   Traits;

typedef Traits::Visibility_heuristic    Visibility_heuristic;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Never_skip_condition                 Skip_condition;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Static_speed_limiter<Kernel>         Enqueue_policy;
typedef CGAL::Surface_mesh_approximate_shortest_path<Traits, Visibility_heuristic, Skip_condition, Enqueue_policy>  Surface_mesh_approximate_shortest_path;

// Surface_mesh_shortest_path
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Surface_mesh>   Shortest_path_traits;
typedef CGAL::Surface_mesh_shortest_path<Shortest_path_traits>          Surface_mesh_shortest_path;

// Heat Method
typedef Surface_mesh::Property_map<vertex_descriptor, double>           Vertex_property_map;

typedef VTKWriter<Kernel, Surface_mesh> Writer;

std::vector<vertex_descriptor>
GetSourceVertices(int number_of_sources)
{
    std::vector<vertex_descriptor> source_vertices(number_of_sources);

    for (int i = 0; i < number_of_sources; i++)
    {
        source_vertices[i] = vertex_descriptor(364*i + 23);
    }

    return source_vertices;
}

std::vector<vertex_descriptor>
GetTargetVertices(int number_of_targets)
{
    std::vector<vertex_descriptor> target_vertices(number_of_targets);

    for (int i = 0; i < number_of_targets; i++)
    {
        target_vertices[i] = vertex_descriptor(1004*i + 692);
    }

    return target_vertices;
}

double
GeodesicDistanceWithApproximateShortestPath(
        Surface_mesh& mesh, vertex_descriptor source_vertex, vertex_descriptor target)
{
    Surface_mesh_approximate_shortest_path shortest_path(mesh);

    shortest_path.add_source(source_vertex);
    shortest_path.add_target(target);
    shortest_path.propagate_geodesic_source();

    double distance = shortest_path.get_geodesic_distance_to_targets()[0];

    return distance;
}

double
GeodesicDistanceWithShortestPath(
        Surface_mesh& mesh, vertex_descriptor source_vertex, vertex_descriptor target)
{
    Surface_mesh_shortest_path shortest_path(mesh);
    shortest_path.add_source_point(source_vertex);

    double distance = shortest_path.shortest_distance_to_source_points(target).first;

    return distance;
}

double
GeodesicDistanceWithHeatMethod(
        Surface_mesh& mesh, vertex_descriptor source_vertex, vertex_descriptor target)
{
    Vertex_property_map vertex_distance_map = mesh.add_property_map<vertex_descriptor,double>().first;

    CGAL::Heat_method_3::estimate_geodesic_distances(mesh, vertex_distance_map, source_vertex);

    double distance = get(vertex_distance_map, target);

    vertex_distance_map.reset();

    return distance;
}

void setZero(std::vector<std::pair<double, double>>& vector_of_pairs)
{
    for (std::pair<double, double> pair : vector_of_pairs)
    {
        pair.first = 0.;
        pair.second = 0.;
    }
}

void LoadMesh(int mesh_number, Surface_mesh& mesh)
{
    std::string filename;
    if (mesh_number == 0)
    {
        std::cout << "====== ARMADILLO ======" << std::endl;
        filename.assign(CGAL::data_file_path("meshes/armadillo.off"));
    }
    else if (mesh_number == 1)
    {
        std::cout << "====== TURBINE ======" << std::endl;
        filename = CGAL::data_file_path("meshes/turbine.off");
    }
    else if (mesh_number == 2)
    {
        std::cout << "====== ELEPHANT ======" << std::endl;
        filename = CGAL::data_file_path("meshes/elephant.off");
    }
    else if (mesh_number == 3)
    {
        std::cout << "====== CHEESE ======" << std::endl;
        filename = CGAL::data_file_path("meshes/cheese.off");
    }

    std::cout << filename << std::endl;

    Surface_mesh tmesh;
    if(!CGAL::IO::read_polygon_mesh(filename, tmesh) ||
        !CGAL::is_triangle_mesh(tmesh))
    {
        std::cerr << "Invalid input file." << std::endl;
    }

    mesh = tmesh;
}

typedef std::pair<double, double> geodesic_result;

int main(int argc, char** argv)
{
    std::chrono::steady_clock::time_point begin, end;

    for (int i = 0; i < 4; i++)
    {
        Surface_mesh mesh;
        LoadMesh(i, mesh);
        int n = 7;

        std::vector<geodesic_result> approx_shopa_results, shopa_results, heat_results;
        double max_shopa = 0.;

        std::vector<vertex_descriptor> source_vertices = GetSourceVertices(n);
        std::cout << "Source vertices: " << std::endl;
        for (auto s : source_vertices)
        {
            std::cout << s << "\t";
        }
        std::cout << std::endl << std::endl;

        std::cout << "Target vertices: " << std::endl;
        std::vector<vertex_descriptor> target_vertices = GetTargetVertices(n);
        for (auto s : target_vertices)
        {
            std::cout << s << "\t";
        }
        std::cout << std::endl << std::endl;

        for (vertex_descriptor source : source_vertices)
        {
            for (vertex_descriptor target : target_vertices)
            {
                std::cout << source << " --> " << target << std::endl;

                // get distances and measure time
                begin = std::chrono::steady_clock::now();
                double appr_dist = GeodesicDistanceWithApproximateShortestPath(mesh, source, target);
                end = std::chrono::steady_clock::now();
                double appr_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                approx_shopa_results.push_back({appr_dist, appr_time});

                begin = std::chrono::steady_clock::now();
                double shopa_dist = GeodesicDistanceWithShortestPath(mesh, source, target);
                end = std::chrono::steady_clock::now();
                double shopa_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                shopa_results.push_back({shopa_dist, shopa_time});
                if (shopa_dist > max_shopa) { max_shopa = shopa_dist; }

                begin = std::chrono::steady_clock::now();
                double heat_dist = GeodesicDistanceWithHeatMethod(mesh, source, target);
                end = std::chrono::steady_clock::now();
                double heat_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                heat_results.push_back({heat_dist, heat_time});

                std::cout << "    Approximate Shortest Path: " << appr_dist << " (" << appr_time << " micros)" << std::endl
                          << "    Shortest Path: " << shopa_dist  << " (" << shopa_time << " micros)"<< std::endl
                          << "    Heat Method: " << heat_dist << " (" << heat_time << " micros)" << std::endl << std::endl;
            }
        }

        std::cout << "longest geodesic: " << max_shopa << std::endl;

        // post-processing: bin results and compute averages
        int bins = 1;
        std::vector<int> bin_sizes(bins);
        std::fill(bin_sizes.begin(), bin_sizes.end(), 0.);

        std::vector<std::pair<double, double>> appr_shopa_bins(bins);
        setZero(appr_shopa_bins);

        std::vector<std::pair<double, double>> shopa_bins(bins);
        setZero(shopa_bins);

        std::vector<std::pair<double, double>> heat_bins(bins);
        setZero(heat_bins);

        for (int i = 0; i < shopa_results.size()-1; i++)
        {
            std::cout << shopa_results[i].first << ", " << shopa_results[i].second;
            int bin_idx = std::floor((shopa_results[i].first - 1e-4) * bins / max_shopa);
            std::cout << " is in bin " << bin_idx << std::endl;
            bin_sizes[bin_idx] ++;

            shopa_bins[bin_idx].first += shopa_results[i].first;
            shopa_bins[bin_idx].second += shopa_results[i].second;

            appr_shopa_bins[bin_idx].first += std::abs(approx_shopa_results[i].first - shopa_results[i].first);
            appr_shopa_bins[bin_idx].second += approx_shopa_results[i].second;

            heat_bins[bin_idx].first += std::abs(heat_results[i].first - shopa_results[i].first);
            heat_bins[bin_idx].second += heat_results[i].second;
        }

        // output bins
        std::cout << std::endl << "++++ average results ++++" << std::endl;
        for (int i = 0; i < bins; i++)
        {
            std::cout << "bin " << i << " [" << i*max_shopa/bins << ", " << (i+1)*max_shopa/bins << "] with size " << bin_sizes[i] << std::endl;

            std::cout << "Shortest Paths: " << std::endl << "\t avg distance: " << shopa_bins[i].first/bin_sizes[i]
                      << std::endl << "\t avg time: " << shopa_bins[i].second/bin_sizes[i] << " micros" << std::endl;
            std::cout << "Approximate Shortest Paths: " << std::endl << "\t avg (abs) error to Shortest Paths: " << appr_shopa_bins[i].first/bin_sizes[i]
                      << std::endl << " \t avg time: " << appr_shopa_bins[i].second/bin_sizes[i]<< " micros" << std::endl;
            std::cout << "Heat method: " << std::endl << "\t avg (abs) error to Shortest Paths: " << heat_bins[i].first/bin_sizes[i]
                      << std::endl << "\t avg time: " << heat_bins[i].second/bin_sizes[i] << " micros" << std::endl;
            std::cout << std::endl << std::endl;
        }
    }

    return 0;
}
