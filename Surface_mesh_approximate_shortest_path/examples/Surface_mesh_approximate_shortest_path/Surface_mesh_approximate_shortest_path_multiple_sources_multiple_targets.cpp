#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximate_shortest_path.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh_approximate_shortest_path_traits<Kernel, Surface_mesh>   Traits;

typedef Traits::Visibility_heuristic    Visibility_heuristic;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Never_skip_condition                 Skip_condition;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Static_speed_limiter<Kernel>         Enqueue_policy;
typedef CGAL::Surface_mesh_approximate_shortest_path<Traits, Visibility_heuristic, Skip_condition, Enqueue_policy>  Surface_mesh_approximate_shortest_path;

typedef boost::graph_traits<Surface_mesh> Graph_traits;

typedef Graph_traits::vertex_descriptor vertex_descriptor;

typedef VTKWriter<Kernel, Surface_mesh> Writer;

int main(int argc, char** argv)
{
    const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/turbine.off");

    Surface_mesh tmesh;
    if(!CGAL::IO::read_polygon_mesh(filename, tmesh) ||
        !CGAL::is_triangle_mesh(tmesh))
    {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "number of mesh faces: " << tmesh.number_of_faces() << std::endl;

    // Shortest Paths
    Surface_mesh_approximate_shortest_path shortest_path(tmesh);

    // add source point
    std::vector<vertex_descriptor> source_vertices = {vertex_descriptor(3),
                                                      vertex_descriptor(972),
                                                      vertex_descriptor(8002),
                                                      vertex_descriptor(3720),
                                                      vertex_descriptor(5204),
                                                      vertex_descriptor(6926),
                                                      vertex_descriptor(2285),
                                                      };
    for (auto source : source_vertices)
    {
        shortest_path.add_source(source);
    }

    // add target points
    std::vector<vertex_descriptor> target_vertices = {vertex_descriptor(5720),
                                                      vertex_descriptor(7134),
                                                      vertex_descriptor(149),
                                                      vertex_descriptor(2801),
                                                      vertex_descriptor(8072),
                                                      vertex_descriptor(4798),
                                                      vertex_descriptor(1193)};
    for (auto target : target_vertices)
    {
        shortest_path.add_target(target);
    }

    // compute geodesic distances
    shortest_path.propagate_geodesic_source();

    std::cout << "The following geodesic distances were computed:" << std::endl;
    std::vector<double> distances_to_targets = shortest_path.get_geodesic_distance_to_targets();
    for (int i = 0; i < target_vertices.size(); i ++)
    {
        std::cout << "\t target " << i << ": " << distances_to_targets[i] << std::endl;
    }

    // write to file
    std::vector<double> distances = shortest_path.get_geodesic_distances();
    Writer writer{};
    writer(tmesh, distances, "turbine_full_geodesics.vtk");

    return 0;
}
