#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximate_shortest_path.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh_approximate_shortest_path_traits<Kernel, Surface_mesh>   Traits;

typedef Traits::Visibility_heuristic    Visibility_heuristic;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Never_skip_condition        Skip_condition;
//typedef CGAL::Surface_mesh_approximate_shortest_path_3::Always_enqueue_in_A<Kernel>         Enqueue_policy;
//typedef CGAL::Surface_mesh_approximate_shortest_path_3::Static_speed_limiter<Kernel>         Enqueue_policy;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Embedding_space_distance_limiter<Kernel>         Enqueue_policy;
typedef CGAL::Surface_mesh_approximate_shortest_path<Traits, Visibility_heuristic, Skip_condition, Enqueue_policy>  Surface_mesh_approximate_shortest_path;

typedef VTKWriter<Kernel, Surface_mesh> Writer;

int main(int argc, char** argv)
{
    const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

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
    Point_3 source_point(0.1796, 0.4966, 0.0785);
    shortest_path.add_source(source_point);
    shortest_path.propagate_geodesic_source();

    // write to file
    std::vector<double> distances = shortest_path.get_geodesic_distances();
    Writer writer{};
    writer(tmesh, distances, "elephant_full_geodesics.vtk");

    return 0;
}
