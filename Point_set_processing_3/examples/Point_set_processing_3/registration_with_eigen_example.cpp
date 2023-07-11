#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#define _USE_MATH_DEFINES
#include <math.h>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Vector_3                                     Vector;
typedef CGAL::Surface_mesh<Point>                            Mesh;
typedef CGAL::Aff_transformation_3<Kernel>                   Transform;
typedef CGAL::Search_traits_3<Kernel>                        TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits>       Neighbor_search;
typedef Neighbor_search::Tree                                Tree;


Mesh merge_meshes(Mesh a, Mesh b) {
    for (auto f : a.faces()) {
        std::array<Mesh::Vertex_index, 3> triangle;
        int i = 0;
        for (Mesh::Vertex_index v : a.vertices_around_face(a.halfedge(f))) {
            triangle[i++] = b.add_vertex(a.point(v));
        }
        b.add_face(triangle[0], triangle[1], triangle[2]);
    }
    return b;
}


std::pair<Transform, Mesh> rigid_registration(const Mesh& source, const Mesh& target,
                                              double w1 = 0.1, double w2 = 1.0, int max_iter = 100) {
    // build k-d tree for nearest neighbor search
    Tree tree(target.points().begin(), target.points().end());
    // resulting transformation is (R, t)
    Eigen::Matrix3d R_all = Eigen::Matrix3d::Identity();
    Eigen::Vector3d t_all = Eigen::Vector3d::Zero();
    int N = source.number_of_vertices();
    auto z = source;
    // solver
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
    cg.setMaxIterations(1000);
    for (int iter = 0; iter < max_iter; ++iter) {
        // build blocks of sparse coefficient matrix A
        Eigen::SparseMatrix<double> A(6 + 3 * N, 6 + 3 * N);
        // column-wise sparsity pattern
        Eigen::VectorXd a = Eigen::VectorXd::Ones(6 + 3 * N);
        a.head(6) *= 6 + 3 * N;
        a.tail(3 * N) *= 6 + 1;
        A.reserve(a);
        // build blocks of result vector b
        Eigen::VectorXd b(6 + 3 * N);
        b.setZero();
        auto start = std::chrono::high_resolution_clock::now();
        for (auto it = source.vertices().begin(); it != source.vertices().end(); ++it) {
            auto i = it->idx();
            Point xp = source.point(*it);
            Eigen::Vector3d x(xp.x(), xp.y(), xp.z());
            // apply transformation to original x
            Eigen::Vector3d x_t = R_all * x + t_all;
            // build a matrix equivalent to cross product
            Eigen::Matrix3d X_t;
            X_t << 0, -x_t(2), x_t(1),
                   x_t(2), 0, -x_t(0),
                   -x_t(1), x_t(0), 0;
            Eigen::Matrix3d XX_t(X_t * X_t);
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    A.coeffRef(j, k) -= XX_t(j, k); // A_rr
                    A.coeffRef(j, 3 + k) += X_t(j, k); // A_tr
                    A.coeffRef(3 + k, j) -= X_t(k, j); // A_rt
                    A.insert(6 + 3 * i + j, k) = X_t(j, k); // A_zir
                    A.insert(k, 6 + 3 * i + j) = -X_t(k, j); // A_rzj
                }
                A.coeffRef(3 + j, 3 + j) += 1; // A_tt
                A.insert(6 + 3 * i + j, 3 + j) = -1; // A_zit
                A.insert(3 + j, 6 + 3 * i + j) = -1; // A_tzi
                A.insert(6 + 3 * i + j, 6 + 3 * i + j) = 1 + w1 / w2; // A_zizj
			}
            b.segment<3>(3) -= x_t;
            // search nearest neighbor
            Point zp = z.point(*it);
            Point pnn = Neighbor_search(tree, zp, 1).begin()->first;
            Eigen::Vector3d PI(pnn.x(), pnn.y(), pnn.z());
            b.segment<3>(6 + 3 * i) = w1 / w2 * PI + x_t;
        }
        auto stop = std::chrono::high_resolution_clock::now();
        // solve system of linear equations
        cg.compute(A);
        Eigen::VectorXd x = cg.solve(b);
        if (cg.info() == Eigen::Success) {
            std::cout << "Convergence achieved!" << std::endl;
        }
        std::cout << "#iterations:     " << cg.iterations() << std::endl;
        std::cout << "estimated error: " << cg.error() << std::endl;
        // update r
        Eigen::Vector3d r(x[0], x[1], x[2]);
        auto roll = Eigen::AngleAxisd(x[0], Eigen::Vector3d::UnitX());
        auto pitch = Eigen::AngleAxisd(x[1], Eigen::Vector3d::UnitY());
        auto yaw = Eigen::AngleAxisd(x[2], Eigen::Vector3d::UnitZ());
        Eigen::Matrix3d R = (roll * pitch * yaw).matrix();
        R_all *= R;
        // update t
        Eigen::Vector3d t(x[3], x[4], x[5]);
        t_all += t;
        // update z
        for (auto it = z.vertices().begin(); it != z.vertices().end(); ++it) {
            auto i = it->idx();
            Point p(x[6 + 3*i], x[6+ 3*i + 1], x[6 + 3*i+2]);
            z.point(*it) = p;
        }
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start).count();
        std::cout << iter << ' ' << duration << '\n' << R_all << '\n' << t_all << std::endl;
    }
    Transform transform(R_all(0, 0), R_all(0, 1), R_all(0, 2), t_all(0),
                        R_all(1, 0), R_all(1, 1), R_all(1, 2), t_all(1),
                        R_all(2, 0), R_all(2, 1), R_all(2, 2), t_all(2), 1);
    return std::make_pair(transform, z);
}

int main(int argc, char* argv[]) {
    const std::string filename1 = argv[1]; //CGAL::data_file_path("meshes/wolf1.ply");
    Mesh mesh1;
    CGAL::IO::read_PLY(filename1, mesh1);
    const std::string filename2 = argv[2]; //CGAL::data_file_path("meshes/wolf2.ply");
    Mesh mesh2;
    CGAL::IO::read_PLY(filename2, mesh2);
    CGAL::draw(merge_meshes(mesh1, mesh2));
    auto result = rigid_registration(mesh1, mesh2, 0.1, 1.0, 10);
    Transform transform = result.first;
    Mesh z = result.second;
    CGAL::Polygon_mesh_processing::transform(transform, mesh1);
    //CGAL::Polygon_mesh_processing::transform(Transform(CGAL::Translation(), Vector(-40, 0, 0)), mesh1);
    CGAL::draw(merge_meshes(mesh1, mesh2));
    return EXIT_SUCCESS;
}
