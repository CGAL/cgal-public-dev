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
//typedef CGAL::Surface_mesh_approximate_shortest_path_3::Always_enqueue_in_A         Enqueue_policy;
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

std::pair<FT, FT> correct_geodesic_dists(face_descriptor face)
{
    if (face.idx() == 4) {
        FT sigma(0.);
        FT d (sqrt( CGAL::square(1./3.-0.9) + CGAL::square(5./3.+0.9) ));
        return std::make_pair(sigma, d);
    }
    else if (face.idx() == 5) {
        Face_values vals;
        FT sigma(sqrt( CGAL::square(1.9) + CGAL::square(0.1) ));
        FT d (sigma + sqrt(5.)/3);
        return std::make_pair(sigma, d);
    }
    else if (face.idx() == 6) {
        FT sigma(sqrt( CGAL::square(1.9) + CGAL::square(0.1) ));
        FT d (sigma + sqrt(5.)/3);
        return std::make_pair(sigma, d);
    }
    else {
        std::cerr << "face index out of bounds, something went wrong.";
    }
}

void test_propagated_face_values(face_descriptor face, Surface_mesh_approximate_shortest_path& shopa)
{
    if (face.idx() >= 4 && face.idx() <= 6)
    {
        FT face_values = shopa.get_geodesic_distances()[face.idx()];
        std::pair<FT, FT> desired_geodesic_dists = correct_geodesic_dists(face);

        //CHECK_CLOSE(face_values.sigma, desired_geodesic_dists.first, 1e-8);
        CHECK_CLOSE(face_values, desired_geodesic_dists.second, 1e-8);
    }
}

int main()
{
    // straight_geodesic_propagation_mesh
    const std::string filename = "../data/right_turn_test_mesh.off";

    Surface_mesh mesh;
    if (!CGAL::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
    {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "running right turn test" << std::endl;
    Surface_mesh_approximate_shortest_path shopa(mesh);
    Point_3 source(FT(0.9), FT(-0.9), FT(1./6.));
    shopa.add_source(source);
    shopa.propagate_geodesic_source();

    // run tests
    for (face_descriptor face : faces(mesh))
    {
        test_propagated_face_values(face, shopa);
    }
    std::cout << "tests successful" << std::endl << std::endl;

    /* // check vertices
    for (vertex_descriptor v : vertices(mesh))
    {
        std::cout << "vertex " << v << " at position " << mesh.point(v) << std::endl;
    }

    // check halfedges
    for (halfedge_descriptor h : halfedges(mesh))
    {
        std::cout << "halfedge " << h << " with source vertex " << mesh.source(h)
                  << " and target vertex " << mesh.target(h) << " is in face " << mesh.face(h)
                  << "." << std::endl;

    }

    // check faces
    for (face_descriptor f : faces(mesh))
    {
        halfedge_descriptor h0 = mesh.halfedge(f);
        std::cout << "face " << f << " with vertices";
        for (vertex_descriptor fv : mesh.vertices_around_face(h0))
        {
            std::cout << " " << fv;
        }
        std::cout << "." << std::endl;
    } */
}
