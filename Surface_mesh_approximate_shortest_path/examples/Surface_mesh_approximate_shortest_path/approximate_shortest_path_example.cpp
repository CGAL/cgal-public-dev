#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximate_shortest_path.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Triangle_mesh;
typedef Surface_mesh_approximate_shortest_path_traits<Kernel, Triangle_mesh> Traits;
typedef Surface_mesh_approximate_shortest_path<Traits> Surface_mesh_approximate_shortest_path_type;
typedef boost::graph_traits<Triangle_mesh> Graph_traits;

typedef Graph_traits::halfedge_iterator halfedge_iterator;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;

typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;

int main(int argc, char** argv)
{
    const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

    Triangle_mesh tmesh;
    if(!CGAL::IO::read_polygon_mesh(filename, tmesh) ||
        !CGAL::is_triangle_mesh(tmesh))
    {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }

    // unfold a triangle in the tangent plane of the triangle with the opposite halfedge
    halfedge_iterator h_it = halfedges(tmesh).begin();
    int h_idx = 10;
    std::advance(h_it, h_idx);

    Surface_mesh_approximate_shortest_path_type shortest(tmesh);

    /*
    // pick up a random face
    const unsigned int randSeed = argc > 2 ? boost::lexical_cast<unsigned int>(argv[2]) : 7915421;
    CGAL::Random rand(randSeed);
    const int target_face_index = rand.get_int(0, static_cast<int>(num_faces(tmesh)));
    face_iterator face_it = faces(tmesh).first;
    std::advance(face_it,target_face_index);
    // ... and define a barycentric coordinates inside the face
    Traits::Barycentric_coordinates face_location = {{0.25, 0.5, 0.25}};

    // construct a shortest path query object and add a source point
    Surface_mesh_shortest_path shortest_paths(tmesh);
    shortest_paths.add_source_point(*face_it, face_location);

    // For all vertices in the tmesh, compute the points of
    // the shortest path to the source point and write them
    // into a file readable using the CGAL Polyhedron demo
    std::ofstream output("shortest_paths_with_id.polylines.txt");
    vertex_iterator vit, vit_end;
    for ( boost::tie(vit, vit_end) = vertices(tmesh);
         vit != vit_end; ++vit)
    {
        std::vector<Traits::Point_3> points;
        shortest_paths.shortest_path_points_to_source_points(*vit, std::back_inserter(points));

        // print the points
        output << points.size() << " ";
        for (std::size_t i = 0; i < points.size(); ++i)
            output << " " << points[i];
        output << std::endl;
    }
    */
    return 0;
}
