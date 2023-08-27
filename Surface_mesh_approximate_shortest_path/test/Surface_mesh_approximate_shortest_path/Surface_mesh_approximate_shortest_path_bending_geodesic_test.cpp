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
//typedef CGAL::Surface_mesh_approximate_shortest_path_3::Never_skip_condition        Skip_condition;
//typedef CGAL::Surface_mesh_approximate_shortest_path_3::Always_enqueue_in_A         Enqueue_policy;
//typedef CGAL::Surface_mesh_approximate_shortest_path_3::Static_speed_limiter         Enqueue_policy;

class Skip_condition
{
    // Implement a skip condition that skips
public:
    Skip_condition() {};

    CGAL::SkipResult operator() ()
    {
        return CGAL::DO_NOT_SKIP;
    }
};

class Enqueue_policy
{
public:
    Enqueue_policy() {};

    CGAL::EnqueueResult operator() (
        double new_geodesic_dist,
        double geodesic_dist_to_overall_target,
        Point_3& prev_face_target_point,
        Point_3& face_target_point,
        Point_3& overall_target_point)
    {
        FT prev_sq_embedding_dist = CGAL::squared_distance(prev_face_target_point, overall_target_point);
        FT sq_embedding_dist = CGAL::squared_distance(face_target_point, overall_target_point);
        FT embedding_dist = CGAL::sqrt(sq_embedding_dist);

        if (sq_embedding_dist <= prev_sq_embedding_dist) {
            return CGAL::ENQUEUE_IN_A; // keep enqueuing in A => these halfedges will be prioritized
        }
        else if (new_geodesic_dist + embedding_dist <= geodesic_dist_to_overall_target) {
            return CGAL::ENQUEUE_IN_B; // enqueue in B => explore later, after we have found a first geodesic distance
        }
        else {
            return CGAL::DO_NOT_ENQUEUE; // this can no longer be a sensible geodesic path anyways
        }
    }
};

typedef CGAL::Surface_mesh_approximate_shortest_path<Traits, Skip_condition, Enqueue_policy>  Surface_mesh_approximate_shortest_path;

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
    if (face.idx() == 7) {
        FT sigma(0.);
        FT d (sqrt(401.)/6.);
        return std::make_pair(sigma, d);
    }
    else if (face.idx() == 10) {
        Face_values vals;
        FT sigma(sqrt(265.)/6);
        FT d (sigma + sqrt(104.)/6.);
        return std::make_pair(sigma, d);
    }
    else if (face.idx() == 12) {
        FT sigma(sqrt(265.)/6 + 1.);
        FT d (sigma + sqrt(8)/3.);
        return std::make_pair(sigma, d);
    }
    else if (face.idx() == 24) {
        FT sigma(sqrt(265.)/6 + 1. + sqrt(2.));
        FT d (sigma + sqrt(29.)/3);
        return std::make_pair(sigma, d);
    }
    else if (face.idx() == 29) {
        FT sigma(sqrt(265.)/6 + 1. + sqrt(2.) + 1.);
        FT d (sigma + sqrt(26.)/3);
        return std::make_pair(sigma, d);
    }
    else {
        std::cerr << "face index out of bounds, something went wrong.";
    }
}

void test_propagated_face_values(face_descriptor face, Surface_mesh_approximate_shortest_path& shopa)
{
    std::set<int> face_idx_set = {7, 10, 12, 24, 29};
    if (face_idx_set.find(face.idx()) != face_idx_set.end())
    {
        Face_values face_values = shopa.get_face_values(face);
        std::pair<FT, FT> desired_geodesic_dists = correct_geodesic_dists(face);

        CHECK_CLOSE(face_values.sigma, desired_geodesic_dists.first, 1e-7);
        CHECK_CLOSE(face_values.d, desired_geodesic_dists.second, 1e-7);
    }
}

int main()
{
    // straight_geodesic_propagation_mesh
    const std::string filename = "../data/bending_geodesic_test_mesh.off";

    Surface_mesh mesh;
    if (!CGAL::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
    {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "running progatation test with bending geodesics" << std::endl;
    Surface_mesh_approximate_shortest_path shopa(mesh);
    Point_3 source(FT(3./2.), FT(13./4.), FT(0.));
    Point_3 target(3.0, 1.0, 0.0);
    shopa.add_target(target);
    shopa.propagate_geodesic_source(source);

    // run tests
    for (face_descriptor face : faces(mesh))
    {
        test_propagated_face_values(face, shopa);
    }
    std::cout << "tests successful" << std::endl << std::endl;
}
