#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Triangle_mesh;
typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
typedef Triangle_mesh::Property_map<vertex_descriptor, double> Vertex_distance_map;
typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh, CGAL::Heat_method_3::Direct> Heat_method_direct;
typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh, CGAL::Heat_method_3::Intrinsic_Delaunay> Heat_method_idt;
typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh, CGAL::Heat_method_3::Intrinsic_mollification_constant> Heat_method_im;

int main(int argc, char* argv[])
{
    // Input filename
    const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("/home/islam/code/cgal-public-dev/Data/data/meshes/larger_sphere.off");
    Triangle_mesh tm;
    if (!CGAL::IO::read_polygon_mesh(filename, tm) || CGAL::is_empty(tm) || !CGAL::is_triangle_mesh(tm))
    {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }

    // Property map for the distance values to the source set
    Vertex_distance_map vertex_distance = tm.add_property_map<vertex_descriptor, double>("v:distance", 0).first;

    // Add the first vertex as the source set
    vertex_descriptor source = *(vertices(tm).first);

    const auto mode = (argc > 2) ? atoi(argv[2]) : 2;
    switch (mode)
    {
    case 0:
    {
      Heat_method_direct hm(tm);

      hm.add_source(source);

      // Estimate geodesic distances
      hm.estimate_geodesic_distances(vertex_distance);
      break;
    }
    case 1:
    {
      Heat_method_idt hm(tm);

      hm.add_source(source);

      // Estimate geodesic distances
      hm.estimate_geodesic_distances(vertex_distance);
      break;
    }

    default:
    {
      Heat_method_im hm(tm);

      hm.add_source(source);

      // Estimate geodesic distances
      hm.estimate_geodesic_distances(vertex_distance);
      break;
    }
    }
    // Compute min and max distance
    double min_distance = std::numeric_limits<double>::max();
    double max_distance = std::numeric_limits<double>::lowest();

    for (vertex_descriptor vd : vertices(tm))
    {
      double dist = get(vertex_distance, vd);
      if (dist < min_distance)
        min_distance = dist;
      if (dist > max_distance)
        max_distance = dist;
    }

    // Property map for colors
    typedef CGAL::Color Color;
    auto vertex_color_map = tm.add_property_map<vertex_descriptor, Color>("v:color", Color(255, 255, 255)).first;

    // Map distances to colors
    for (vertex_descriptor vd : vertices(tm)) {
        double dist = get(vertex_distance, vd);
        double normalized_dist = (dist - min_distance) / (max_distance - min_distance);
        unsigned char r = static_cast<unsigned char>(255 * normalized_dist);
        unsigned char g = 0;
        unsigned char b = static_cast<unsigned char>(255 * (1 - normalized_dist));  // Simple gradient from red to blue
        put(vertex_color_map, vd, Color(r, g, b));
    }

    // Save the mesh with distance information to a file
    const std::string output_filename = "output_with_distances.ply";
    if (!CGAL::IO::write_polygon_mesh(output_filename, tm))
    {
        std::cerr << "Failed to write output file." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Geodesic distances saved to " << output_filename << std::endl;
    return 0;
}
