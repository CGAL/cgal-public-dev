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
typedef CGAL::Surface_mesh_approximate_shortest_path_3::Always_enqueue_in_A         Enqueue_policy;

typedef CGAL::Surface_mesh_approximate_shortest_path<Traits, Skip_condition, Enqueue_policy>  Surface_mesh_approximate_shortest_path;

typedef Surface_mesh_approximate_shortest_path::unfold_triangle_3           unfold_triangle_3;
typedef Surface_mesh_approximate_shortest_path::Construct_heuristic_point_2    Construct_heuristic_point_2;

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

int read_test_mesh(Surface_mesh& mesh)
{
    const std::string filename = "../data/traits_test_mesh.off";

    if(!CGAL::IO::read_polygon_mesh(filename, mesh) ||
        !CGAL::is_triangle_mesh(mesh))
    {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }

    return 0;
}

void test_unfolding()
{
    Surface_mesh mesh;
    read_test_mesh(mesh);

    // unfold the test_mesh (traits_test_mesh.off) along select halfedges
    Surface_mesh_approximate_shortest_path appr_shortest_path(mesh);
    auto edge_lengths = appr_shortest_path.Get_edge_length_map();

    // halfedge h1
    halfedge_descriptor h1(1);
    auto untri_h1 = unfold_triangle_3()(mesh, h1, edge_lengths);
    CHECK_CLOSE(untri_h1.B.x(), FT(1.5), FT(1e-10));
    CHECK_CLOSE(untri_h1.B.y(), FT(0.), FT(1e-10));
    CHECK_CLOSE(untri_h1.P.x(), FT(1.), FT(1e-10));
    CHECK_CLOSE(untri_h1.P.y(), FT(2.5), FT(1e-10));

    // halfedge h3
    halfedge_descriptor h3(3);
    auto untri_h3 = unfold_triangle_3()(mesh, h3, edge_lengths);
    CHECK_CLOSE(untri_h3.B.x(), CGAL::sqrt(FT(6.5)), FT(1e-10));
    CHECK_CLOSE(untri_h3.B.y(), FT(0.), FT(1e-10));
    FT Px_h3(1.5/(2.*CGAL::sqrt(6.5)));
    CHECK_CLOSE(untri_h3.P.x(), Px_h3, FT(1e-10));
    CHECK_CLOSE(untri_h3.P.y(), CGAL::sqrt(2.25-CGAL::square(Px_h3)), FT(1e-10));

    // halfedge h1
    halfedge_descriptor h5(5);
    auto untri_h5 = unfold_triangle_3()(mesh, h5, edge_lengths);
    CHECK_CLOSE(untri_h5.B.x(), CGAL::sqrt(FT(7.25)), FT(1e-10));
    CHECK_CLOSE(untri_h5.B.y(), FT(0.), FT(1e-10));
    FT Px_h5(11.5/(2.*CGAL::sqrt(7.25)));
    CHECK_CLOSE(untri_h5.P.x(), Px_h5, FT(1e-10));
    CHECK_CLOSE(untri_h5.P.y(), CGAL::sqrt(6.5-CGAL::square(Px_h5)), FT(1e-10));
}

void test_heuristic_point_construction()
{
    Surface_mesh mesh;
    read_test_mesh(mesh);

    // unfold the test_mesh (traits_test_mesh.off) along select halfedges
    Surface_mesh_approximate_shortest_path appr_shortest_path(mesh);
    auto edge_lengths = appr_shortest_path.Get_edge_length_map();

    // halfedge h1
    halfedge_descriptor h1(1);
    //auto Q_h1 = appr_shortest_path.construct_heuristic_point_object()(mesh, h1, edge_lengths);
    auto Q_h1 = Construct_heuristic_point_2()(mesh, h1, edge_lengths);
    CHECK_CLOSE(Q_h1.x(), FT(0.83333333333333), FT(1e-10));
    CHECK_CLOSE(Q_h1.y(), FT(0.83333333333333), FT(1e-10));

    // halfedge h3
    halfedge_descriptor h3(3);
    auto Q_h3 = Construct_heuristic_point_2()(mesh, h3, edge_lengths);
    CHECK_CLOSE(Q_h3.x(), FT(0.66616074204193), FT(1e-10));
    CHECK_CLOSE(Q_h3.y(), FT(0.91289120453441), FT(1e-10));

    // halfedge h1
    halfedge_descriptor h5(5);
    auto Q_h5 = Construct_heuristic_point_2()(mesh, h5, edge_lengths);
    CHECK_CLOSE(Q_h5.x(), FT(1.83610877748221), FT(1e-10));
    CHECK_CLOSE(Q_h5.y(), FT(0.86438395711515), FT(1e-10));
}

void test_source_reconstruction()
{

}

void test_intersection_tests()
{

}

int main()
{
    const std::string filename = "../data/traits_test_mesh.off";

    Surface_mesh mesh;
    if(!CGAL::IO::read_polygon_mesh(filename, mesh) ||
        !CGAL::is_triangle_mesh(mesh))
    {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "running traits tests" << std::endl;
    test_unfolding();
    test_heuristic_point_construction();
    std::cout << "tests successful" << std::endl << std::endl;
}


