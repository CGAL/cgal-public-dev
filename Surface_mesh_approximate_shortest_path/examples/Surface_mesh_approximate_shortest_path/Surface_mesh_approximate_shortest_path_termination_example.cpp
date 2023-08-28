#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximate_shortest_path.h>
#include <iostream>
#include <chrono>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh_approximate_shortest_path_traits<Kernel, Surface_mesh>   Traits;

typedef Traits::Visibility_heuristic    Visibility_heuristic;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Never_skip_condition        Skip_condition;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Embedding_space_distance_limiter<Kernel>         Enqueue_policy;
typedef CGAL::Surface_mesh_approximate_shortest_path<Traits, Visibility_heuristic, Skip_condition, Enqueue_policy>  Surface_mesh_approximate_shortest_path;

typedef boost::graph_traits<Surface_mesh> Graph_traits;

typedef Graph_traits::vertex_descriptor vertex_descriptor;

typedef VTKWriter<Kernel, Surface_mesh> Writer;

int main(int argc, char** argv)
{
    const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/armadillo.off");

    Surface_mesh tmesh;
    if(!CGAL::IO::read_polygon_mesh(filename, tmesh) ||
        !CGAL::is_triangle_mesh(tmesh))
    {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "number of mesh faces: " << tmesh.number_of_faces() << std::endl;

    // Shortest Paths
    Surface_mesh_approximate_shortest_path shortest_path(tmesh, true, true);

    // add source point
    Point_3 source_point(tmesh.point(vertex_descriptor(3143)));
    shortest_path.add_source(source_point);

    Point_3 target_point(tmesh.point(vertex_descriptor(3139)));
    shortest_path.add_target(target_point);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    shortest_path.propagate_geodesic_source();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double appr_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    double dist = shortest_path.get_geodesic_distance_to_targets()[0];

    int total_updates = shortest_path.number_of_total_face_updates();
    int iter_count = shortest_path.get_iter_count();

    // output
    std::cout << "The geodesic distance is " << dist << " (" << appr_time << " micros)" << std::endl;
    std::cout << "The number of total face value updates is " << total_updates
              << " over " << iter_count << " iterations." << std::endl;

    // write to file
    std::vector<double> distances = shortest_path.get_geodesic_distances();
    Writer writer{};
    writer(tmesh, distances, "armadillo_partial_geodesic.vtk");

    std::cout << "writing geodesic indicator map" << std::endl;
    std::vector<bool> indicator_map = shortest_path.extract_shortest_path_face_indicator_map(0);
    writer(tmesh, indicator_map, "armadillo_geodesic_path.vtk");

    for (int i = 0; i < iter_count; i++)
    {
        std::string filename = "face_updates_after" + std::to_string(i) +"th_iteration.vtk";
        std::vector<int> updates = shortest_path.get_cumulative_update_counts_after_iterations(i);
        writer(tmesh, updates, filename);
    }

    return 0;
}
