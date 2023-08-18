#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximate_shortest_path.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh_approximate_shortest_path_traits<Kernel, Surface_mesh>   Traits;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Never_skip_condition        Skip_condition;
//typedef CGAL::Surface_mesh_approximate_shortest_path_3::Always_enqueue_in_A         Enqueue_policy;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Static_speed_limiter         Enqueue_policy;

typedef CGAL::Surface_mesh_approximate_shortest_path<Traits, Skip_condition, Enqueue_policy>  Surface_mesh_approximate_shortest_path;

typedef boost::graph_traits<Surface_mesh> Graph_traits;

typedef Graph_traits::vertex_descriptor vertex_descriptor;
typedef Graph_traits::edge_descriptor edge_descriptor;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;
typedef Graph_traits::face_descriptor face_descriptor;

typedef Surface_mesh_approximate_shortest_path::Face_values Face_values;

void
WriteVTK(const char* filename, Surface_mesh& mesh, std::vector<double> face_data)
{
    std::ofstream out(filename);

    // header for (legacy) vtk file
    out << "# vtk DataFile Version 2.0\n";
    out << "description: data on regular grid\n";
    out << "ASCII\n";

    // write geometry/topology
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << mesh.num_vertices() << " float\n";

    for (vertex_descriptor vd : vertices(mesh))
    {
        Point_3 v = mesh.point(vd);
        out << v.x() << " "
            << v.y() << " "
            << v.z() << std::endl;
    }

    // 4*_num_faces is the "size of the cell list"
    // see https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    out << "CELLS " << mesh.num_faces() << " " << 4 * mesh.num_faces() << std::endl;
    for (face_descriptor fd : faces(mesh))
    {
        out << "3 ";
        halfedge_descriptor h0 = mesh.halfedge(fd);
        for (halfedge_descriptor hd : CGAL::halfedges_around_face(h0, mesh))
        {
            out << mesh.source(hd).idx() << " ";
        }
        out << std::endl;
    }

    // write cell types (5 = VTK_TRIANGLE)
    out << "CELL_TYPES " << mesh.num_faces() << std::endl;
    for (int face_num = 0; face_num < mesh.num_faces(); ++face_num)
    {
        out << "5" << std::endl;
    }

    out << "CELL_DATA " << mesh.num_faces() << std::endl;
    out << "SCALARS " << "cell_scalars " << "float " << "1" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;

    for (int i = 0; i < mesh.num_faces(); ++i)
    {
        out << face_data[i] << std::endl;
    }

    out.close();
}

std::vector<double>
ExtractDistanceData(Surface_mesh& mesh, Surface_mesh_approximate_shortest_path& shopa)
{
    std::vector<double> distances(mesh.num_faces());
    for (face_descriptor fd : faces(mesh))
    {
        distances[fd.idx()] = shopa.get_face_values(fd).d;
    }
    return distances;
}

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

    // Shortest Paths
    Surface_mesh_approximate_shortest_path shortest_path(tmesh);

    // add source point
    Point_3 source_point(0.1796, 0.4966, 0.0785);
    shortest_path.propagate_geodesic_source(source_point);

    // write to file
    std::vector<double> distances = ExtractDistanceData(tmesh, shortest_path);
    WriteVTK("turbine_geodesics.vtk", tmesh, distances);

    return 0;
}
