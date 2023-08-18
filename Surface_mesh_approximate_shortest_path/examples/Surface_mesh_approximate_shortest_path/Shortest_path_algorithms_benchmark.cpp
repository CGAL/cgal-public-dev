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

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef Kernel::Point_3 Point_3;

typedef boost::graph_traits<Surface_mesh> Graph_traits;

typedef Graph_traits::vertex_descriptor vertex_descriptor;
typedef Graph_traits::edge_descriptor edge_descriptor;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;
typedef Graph_traits::face_descriptor face_descriptor;

// Surface_mesh_approximate_shortest_path
typedef CGAL::Surface_mesh_approximate_shortest_path_traits<Kernel, Surface_mesh>   Traits;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Never_skip_condition        Skip_condition;
//typedef CGAL::Surface_mesh_approximate_shortest_path_3::Always_enqueue_in_A         Enqueue_policy;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Static_speed_limiter         Enqueue_policy;

typedef CGAL::Surface_mesh_approximate_shortest_path<Traits, Skip_condition, Enqueue_policy>  Surface_mesh_approximate_shortest_path;

// Surface_mesh_shortest_path
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Surface_mesh>   Shortest_path_traits;
typedef CGAL::Surface_mesh_shortest_path<Shortest_path_traits>          Surface_mesh_shortest_path;

// Heat Method
typedef Surface_mesh::Property_map<vertex_descriptor, double>           Vertex_property_map;

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

std::vector<vertex_descriptor>
GetTargetVertices(int number_of_targets)
{
    std::vector<vertex_descriptor> target_vertices(number_of_targets);

    for (int i = 0; i < number_of_targets; i++)
    {
        target_vertices.push_back(vertex_descriptor(100*i));
    }
}

std::vector<double>
GeodesicDistanceWithApproximateShortestPath(
        Surface_mesh& mesh, vertex_descriptor source_vertex)
{
    Surface_mesh_approximate_shortest_path shortest_path(mesh);

    shortest_path.propagate_geodesic_source(mesh.point(source_vertex));

    // write to file
    std::vector<double> distances = ExtractDistanceData(mesh, shortest_path);
    //WriteVTK("armadillo_geodesics.vtk", mesh, distances);

    return distances;
}

std::vector<double>
GeodesicDistanceWithShortestPath(
        Surface_mesh& mesh, vertex_descriptor source_vertex)
{
    Surface_mesh_shortest_path shortest_path(mesh);
    shortest_path.add_source_point(source_vertex);

    std::vector<vertex_descriptor> target_vertices;
    target_vertices = GetTargetVertices(10);

    std::vector<double> distances;
    for (vertex_descriptor target_vertex : target_vertices)
    {
        double dist = shortest_path.shortest_distance_to_source_points(target_vertex).first;
        distances.push_back(dist);
    }

    return distances;
}

std::vector<double>
GeodesicDistanceWithHeatMethod(
        Surface_mesh& mesh, vertex_descriptor source_vertex)
{
    Vertex_property_map vertex_distance_map = mesh.add_property_map<vertex_descriptor,double>().first;

    CGAL::Heat_method_3::estimate_geodesic_distances(mesh, vertex_distance_map, source_vertex);

    std::vector<vertex_descriptor> target_vertices = GetTargetVertices(10);
    std::vector<double> distances;
    for(vertex_descriptor target : target_vertices)
    {
        distances.push_back( get(vertex_distance_map, target) );
    }

    return distances;
}

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
    Surface_mesh_approximate_shortest_path shortest_path(tmesh);

    // add source point
    Point_3 source_point(0.1796, 0.4966, 0.0785);
    shortest_path.propagate_geodesic_source(source_point);

    // write to file
    std::vector<double> distances = ExtractDistanceData(tmesh, shortest_path);
    WriteVTK("armadillo_geodesics.vtk", tmesh, distances);

    return 0;
}
