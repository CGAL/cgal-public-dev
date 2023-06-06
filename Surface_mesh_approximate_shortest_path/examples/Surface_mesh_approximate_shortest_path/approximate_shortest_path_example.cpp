#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximate_shortest_path.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Triangle_mesh;

typedef CGAL::Surface_mesh_approximate_shortest_path_traits<Kernel, Triangle_mesh> Traits;
typedef CGAL::Surface_mesh_approximate_shortest_path<Traits> Surface_mesh_approximate_shortest_path;

typedef boost::graph_traits<Triangle_mesh> Graph_traits;

typedef Graph_traits::halfedge_iterator halfedge_iterator;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;

typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;

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

    // unfold a triangle in the tangent plane of the triangle with the opposite halfedge
    halfedge_iterator h_it = halfedges(tmesh).begin();
    int h_idx = 10;
    std::advance(h_it, h_idx);

    Surface_mesh_approximate_shortest_path shortest_path(tmesh);
    //auto unfolded_triangle =
    //std::cout << "unfolded triangle" << std::endl;
    //std::cout << "B: " << unfolded_triangle.B[0] << "," << unfolded_triangle.B[1] << std::endl;
    //std::cout << "P: " << unfolded_triangle.P[0] << "," << unfolded_triangle.P[1] << std::endl;

    return 0;
}
