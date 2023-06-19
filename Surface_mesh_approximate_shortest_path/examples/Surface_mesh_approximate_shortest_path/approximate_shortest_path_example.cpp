#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximate_shortest_path.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Triangle_mesh;
typedef Kernel::Point_2 Point_2;

typedef CGAL::Surface_mesh_approximate_shortest_path_traits<Kernel, Triangle_mesh> Traits;
typedef CGAL::Surface_mesh_approximate_shortest_path<Traits> Surface_mesh_approximate_shortest_path;

typedef boost::graph_traits<Triangle_mesh> Graph_traits;

typedef Graph_traits::edge_descriptor edge_descriptor;
typedef Graph_traits::halfedge_iterator halfedge_iterator;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;
typedef Graph_traits::face_descriptor face_descriptor;

typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;

typedef Surface_mesh_approximate_shortest_path::Face_values Face_values;

typedef CGAL::Surface_mesh_approximate_shortest_path_3::Compute_squared_edge_length<Kernel, Triangle_mesh> Compute_squared_edge_length;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Unfold_triangle_3_along_halfedge<Kernel, Triangle_mesh> unfold_object;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Reconstruct_source_point_in_triangle_tangent_space<Kernel, Triangle_mesh> reconstructed_source;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Construct_triangle_centroid_2<Kernel> Construct_triangle_centroid;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Construct_heuristic_point_2<Kernel, Triangle_mesh> Construct_heuristic_point_2;

int main(int argc, char** argv)
{
    const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/triangle.off");

    Triangle_mesh tmesh;
    if(!CGAL::IO::read_polygon_mesh(filename, tmesh) ||
        !CGAL::is_triangle_mesh(tmesh))
    {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "This is the mesh we consider (the reference triangle)" << std::endl;
    std::cout << tmesh;

    // Let's play around with property maps
    Triangle_mesh::Property_map<edge_descriptor, Kernel::FT> squared_edge_lengths;
    bool created;
    boost::tie(squared_edge_lengths, created) = tmesh.add_property_map<edge_descriptor, Kernel::FT>("edges");
    assert(created);

    // get face property map
    Triangle_mesh::Property_map<face_descriptor, Face_values> face_values_map;
    bool face_created;
    boost::tie(face_values_map, face_created) = tmesh.add_property_map<face_descriptor, Face_values>("faces");
    assert(face_created);

    // Shortest Paths
    Surface_mesh_approximate_shortest_path shortest_path(tmesh);
    //unfold_object unfolded_triangle = shortest_path.unfold_triangle_3_along_halfedge_object();
    halfedge_descriptor h(0);
    std::cout << "for halfedge " << h.idx() << std::endl;
    //auto tri = unfolded_triangle(tmesh, h, squared_edge_lengths);

    // construct centroid of the triangle
    //Construct_triangle_centroid C = shortest_path.construct_centroid_object();
    //Point_2 center = C(tri.B, tri.P);

    // construct heuristic point
    Point_2 Q = shortest_path.construct_heuristic_point_object()(tmesh, h, squared_edge_lengths);

    // reconstruct source
    Point_2 S = shortest_path.reconstruct_source_point_in_triangle_tangent_space_object()(tmesh, h, squared_edge_lengths, face_values_map);

    // do some intersection tests
    Point_2 source(2., -0.5);
    Point_2 heu(0.5, 0.5);
    Point_2 A(0., 0.);
    Point_2 B(1., 0.);

    auto inter = shortest_path.edge_intersection_test_object()(source, heu, A, B);
    std::cout << "intersection test yields:" << std::endl
              << "intersection: " << inter.first << "\t right: " << inter.second << std::endl;

    return 0;
}
