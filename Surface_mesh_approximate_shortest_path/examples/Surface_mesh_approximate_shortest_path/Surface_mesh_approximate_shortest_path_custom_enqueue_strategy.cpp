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

/*!
 * \brief This particular Enqueue policy prefers face whose barycenter has high z coordinates
 * by keeping them A. Thus, an algorithm using the enqueue policy will have a tendency to find
 * vertical geodesics faster.
 */
template<class Kernel>
class Enqueue_policy
{
public:
    typedef typename Kernel::FT         FT;
    typedef typename Kernel::Point_3    Point_3;

public:
    Enqueue_policy() {};

    CGAL::EnqueueResult operator() (FT geodesic_dist, FT geodesic_radius,
                             Point_3 prev_face_target = Point_3(),
                             Point_3 curr_face_target = Point_3(),
                             Point_3 overall_target = Point_3())
    {
        if (prev_face_target.z() < curr_face_target.z())
        {
            return CGAL::ENQUEUE_IN_A;
        }
        else
        {
            return CGAL::ENQUEUE_IN_B;
        }
    }
};

typedef CGAL::Surface_mesh_approximate_shortest_path<Traits, Visibility_heuristic, Skip_condition, Enqueue_policy<Kernel>>  Surface_mesh_approximate_shortest_path;

typedef boost::graph_traits<Surface_mesh> Graph_traits;

typedef Graph_traits::vertex_descriptor vertex_descriptor;

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

    // Shortest Paths
    Surface_mesh_approximate_shortest_path shortest_path(tmesh);

    // add source point
    Point_3 source_point(tmesh.point(vertex_descriptor(294)));
    shortest_path.add_source(source_point);

    shortest_path.propagate_geodesic_source();

    // write to file
    std::string outfilename = "elephant_full_geodesics.vtk";
    std::cout << "writing results to " << outfilename << std::endl;
    std::vector<double> distances = shortest_path.get_geodesic_distances();
    Writer writer{};
    writer(tmesh, distances, outfilename);

    return 0;
}
