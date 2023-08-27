#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximate_shortest_path.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef Kernel::Point_3 Point_3;
typedef Kernel::FT FT;

typedef CGAL::Surface_mesh_approximate_shortest_path_traits<Kernel, Surface_mesh>   Traits;

typedef Traits::Visibility_heuristic    Visibility_heuristic;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Never_skip_condition        Skip_condition;
//typedef CGAL::Surface_mesh_approximate_shortest_path_3::Always_enqueue_in_A<Kernel>         Enqueue_policy;
//typedef CGAL::Surface_mesh_approximate_shortest_path_3::Static_speed_limiter<Kernel>         Enqueue_policy;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Embedding_space_distance_limiter<Kernel>         Enqueue_policy;
typedef CGAL::Surface_mesh_approximate_shortest_path<Traits, Visibility_heuristic, Skip_condition, Enqueue_policy>  Surface_mesh_approximate_shortest_path;

typedef boost::graph_traits<Surface_mesh> Graph_traits;

typedef Graph_traits::vertex_descriptor vertex_descriptor;
typedef Graph_traits::edge_descriptor edge_descriptor;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;
typedef Graph_traits::face_descriptor face_descriptor;

typedef VTKWriter<Kernel, Surface_mesh> Writer;

typedef Surface_mesh_approximate_shortest_path::Face_values Face_values;


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

    // Shortest Paths
    Surface_mesh_approximate_shortest_path shortest_path(tmesh, true, true);

    // add source point
    Point_3 source_point(tmesh.point(vertex_descriptor(3143)));
    shortest_path.add_source(source_point);

    Point_3 target_point(tmesh.point(vertex_descriptor(3139)));
    shortest_path.add_target(target_point);

    shortest_path.propagate_geodesic_source();

    // output
    std::cout << "Geodesic Distances to Target Points:" << std::endl;
    std::vector<FT> distances_to_target = shortest_path.get_geodesic_distance_to_targets();
    for (FT d : distances_to_target)
    {
        std::cout << d << std::endl;
    }

    // write to file
    std::vector<double> distances = shortest_path.get_geodesic_distances();
    Writer writer{};
    writer(tmesh, distances, "armadillo_full_geodesics.vtk");

    std::cout << "writing indicator map" << std::endl;
    std::vector<bool> indicator_map = shortest_path.extract_shortest_path_face_indicator_map(0);
    writer(tmesh, indicator_map, "armadillo_geodesic_path.vtk");

    return 0;
}
