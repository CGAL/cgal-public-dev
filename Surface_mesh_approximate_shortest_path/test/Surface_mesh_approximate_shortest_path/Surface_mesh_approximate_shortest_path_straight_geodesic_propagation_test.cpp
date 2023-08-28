#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximate_shortest_path.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT                                          FT;
typedef Kernel::Point_2                                     Point_2;
typedef Kernel::Point_3                                     Point_3;
typedef CGAL::Surface_mesh<Point_3>                         Surface_mesh;

typedef boost::graph_traits<Surface_mesh> Graph_traits;

typedef Graph_traits::vertex_descriptor vertex_descriptor;
typedef Graph_traits::edge_descriptor edge_descriptor;
typedef Graph_traits::halfedge_iterator halfedge_iterator;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;
typedef Graph_traits::face_descriptor face_descriptor;

typedef CGAL::Surface_mesh_approximate_shortest_path_traits<Kernel, Surface_mesh>   Traits;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Never_skip_condition        Skip_condition;
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Static_speed_limiter<Kernel>         Enqueue_policy;

typedef Traits::Visibility_heuristic    Visibility_heuristic;
typedef CGAL::Surface_mesh_approximate_shortest_path<Traits, Visibility_heuristic, Skip_condition, Enqueue_policy>  Surface_mesh_approximate_shortest_path;

typedef CGAL::Face_values<Kernel>   Face_values;

template <typename FT, typename FT2>
void CHECK_EQUAL(const FT& a, const FT2& b)
{
    if (a != b)
        std::cerr << "ERROR: a (" << a << ") is not equal to b (" << b << ").\n";

    assert(a == b);
}

template <typename FT, typename FT2>
void CHECK_CLOSE(const FT& a, const FT2& b, const FT& eps)
{
    if ( !(CGAL::abs(a-b) < eps) )
        std::cerr << "ERROR: difference (" << CGAL::abs(a-b) << ") is larger than eps (" << eps << ").\n";
    assert(CGAL::abs(a-b) < eps);
}

Face_values correct_propagated_face_values(face_descriptor face)
{
    if (face.idx() == 0) {
        FT sigma(0.);
        FT d(2./3.);
        std::array<FT,3> d2verts = {FT(9.), FT(4.), FT(5.)};
        return Face_values(sigma, d, d2verts);
    }
    else if (face.idx() == 1) {
        FT sigma(0.);
        FT d = sqrt(FT(4.+16./9.));
        std::array<FT,3> d2verts = {FT(20.), FT(9.), FT(5.)};
        return Face_values(sigma, d, d2verts);
    }
    else if (face.idx() == 2) {
        FT sigma(0.);
        FT d = sqrt(FT(25./9.+9.));
        std::array<FT,3> d2verts = {FT(29.), FT(20.), FT(5.)};
        return Face_values(sigma, d, d2verts);
    }
    else if (face.idx() == 3) {
        FT sigma(0.);
        FT d = sqrt(FT(196./9.+169./9.));
        std::array<FT,3> d2verts = {FT(100.), FT(20.), FT(29.)};
        return Face_values(sigma, d, d2verts);
    }
    else {
        std::cerr << "face index out of bounds, something went wrong.";
        return EXIT_FAILURE;
    }
}

void test_propagated_face_values(face_descriptor face, Surface_mesh_approximate_shortest_path& shopa)
{
    Face_values values = shopa.get_geodesic_distances()[face.idx()];
    Face_values desired_face_values = correct_propagated_face_values(face);

    CHECK_CLOSE(values.sigma, desired_face_values.sigma, 1e-10);
    //CHECK_CLOSE(values.d, desired_face_values.d, 1e-10);

    //CHECK_CLOSE(values.d2verts[0], desired_face_values.d2verts[0], 1e-10);
    //CHECK_CLOSE(values.d2verts[1], desired_face_values.d2verts[1], 1e-10);
    //CHECK_CLOSE(values.d2verts[2], desired_face_values.d2verts[2], 1e-10);
}

int main()
{
    // straight_geodesic_propagation_mesh
    const std::string filename = "../data/straight_geodesic_propagation_mesh.off";
    //const std::string filename = "../data/traits_test_mesh.off";

    Surface_mesh mesh;
    if (!CGAL::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
    {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "running progatation test for straight geodesic (no turns, only intersections)" << std::endl;
    Surface_mesh_approximate_shortest_path shopa(mesh);
    Point_3 source(FT(0.), FT(0.), FT(0.));
    shopa.add_source(source);
    shopa.propagate_geodesic_source();

    // run tests
    for (face_descriptor face : faces(mesh))
    {
        test_propagated_face_values(face, shopa);
    }
    std::cout << "tests successful" << std::endl << std::endl;
}
