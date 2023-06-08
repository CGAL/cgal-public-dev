#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/draw_point_set_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
#include <chrono>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>


typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Point_set_3<Point>                             PointSet;
typedef CGAL::Aff_transformation_3<Kernel>                   Transform;
//typedef CGAL::Search_traits_3<Kernel>                        TreeTraits;
//typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits>       Neighbor_search;
//typedef Neighbor_search::Tree                                Tree;


int closest_point(const PointSet& sp, const Eigen::Vector3f& p) {
    int closest = 0;
    float min_dist = std::numeric_limits<float>::max();
    for (auto i = 0; i < sp.size(); ++i) {
        Eigen::Vector3f p_other(sp.point(i)[0], sp.point(i)[1], sp.point(i)[2]);
        auto dist = (p - p_other).norm();
        if (dist < min_dist) {
            min_dist = dist;
            closest = i;
        }
    }
    return closest;
}


Transform register_point_sets(const PointSet& source, const PointSet& target, float w1 = 1.0, float w2 = 1.0, bool point2plane = true) {

    float ratio = w1 / w2;
    // intermediate representation
    auto z = source;
    // build tree
    // Tree tree(target.points().begin(), target.points().end());
    // initialize rigid transformation
    Eigen::Matrix3f R_all = Eigen::Matrix3f::Identity();
    Eigen::Vector3f t_all = Eigen::Vector3f::Zero();
    for (auto time = 0; time < 100; ++time) {
        // matrix A
        Eigen::SparseMatrix<float> A(6 + 3 * z.size(), 6 + 3 * z.size());
        // column-wise sparsity pattern
        Eigen::VectorXf a = Eigen::VectorXf::Ones(6 + 3 * z.size());
        a.head(6) *= 6 + 3 * z.size();
        a.tail(3 * z.size()) *= 9;
        A.reserve(a);
        // vector b
        Eigen::VectorXf b(6 + 3 * z.size());
        b.setZero();
        // fill matrix A and vector b
        auto start = std::chrono::high_resolution_clock::now();
        for (auto i = 0; i < z.size(); ++i) {
            Eigen::Vector3f sp(source.point(i)[0], source.point(i)[1], source.point(i)[2]);
            Eigen::Vector3f p = R_all * sp + t_all;
            //Eigen::Vector3f p = Eigen::Vector3f(z.point(i)[0], z.point(i)[1], z.point(i)[2]);
            Eigen::Matrix3f Z;
            Z << 0, -p[2], p[1], p[2], 0, -p[0], -p[1], p[0], 0;
            Eigen::Matrix3f A_rr = Z * Z;
            Eigen::Matrix3f dyad = ratio * Eigen::Matrix3f::Identity();
            Eigen::Vector3f zp(z.point(i)[0], z.point(i)[1], z.point(i)[2]);
            int closest = closest_point(target, zp);
            if (point2plane) {
                Eigen::Vector3f tn(target.normal(closest)[0], target.normal(closest)[1], target.normal(closest)[2]);
                dyad *= tn * tn.transpose();
            }
            Eigen::Matrix3f A_zz = Eigen::Matrix3f::Identity() + dyad;
            for (auto j = 0; j < 3; ++j) {
                for (auto k = 0; k < 3; ++k) {
                    A.coeffRef(j, k) += A_rr(j, k);
                    A.coeffRef(j, 3 + k) -= Z(j, k);
                    A.coeffRef(3 + k, j) -= Z(j, k);
                    A.insert(6 + 3 * i + j, k) = Z(j, k);
                    A.insert(k, 6 + 3 * i + j) = Z(j, k);
                    A.insert(6 + 3 * i + j, 6 + 3 * i + k) = A_zz(j, k);
                }
                A.coeffRef(3 + j, 3 + j) += 1;
                A.insert(6 + 3 * i + j, 3 + j) = -1;
                A.insert(3 + j, 6 + 3 * i + j) = -1;
            }
            b.segment<3>(3) -= p;
            //Neighbor_search search(tree, z.point(i), 1);
            //Point pnn = search.begin()->first;
            Eigen::Vector3f tp(target.point(closest)[0], target.point(closest)[1], target.point(closest)[2]);
            b.segment<3>(6 + 3 * i) = p + dyad * tp;
        }
        auto stop = std::chrono::high_resolution_clock::now();
        std::cout << time << ' ' << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << std::endl;
        // solve the linear system
        Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> cg;
        cg.setMaxIterations(20);
        cg.compute(A);
        Eigen::VectorXf x = cg.solve(b);
        if (cg.info() == Eigen::Success) {
            std::cout << "Convergence achieved!" << std::endl;
        }
        std::cout << "#iterations:     " << cg.iterations() << std::endl;
        std::cout << "estimated error: " << cg.error() << std::endl;
        Eigen::Vector3f r(x[0], x[1], x[2]);
        auto roll = Eigen::AngleAxisf(x[0], Eigen::Vector3f::UnitX());
        auto pitch = Eigen::AngleAxisf(x[1], Eigen::Vector3f::UnitY());
        auto yaw = Eigen::AngleAxisf(x[2], Eigen::Vector3f::UnitZ());
        Eigen::Matrix3f R = (roll * pitch * yaw).matrix();
        Eigen::Vector3f t(x[3], x[4], x[5]);
        //std::cout << "R = " << R << std::endl;
        //std::cout << "t = " << t << std::endl;
        R_all *= R;
        t_all += t;
        std::cout << "R_all = " << std::endl << R_all << std::endl;
        std::cout << "t_all = " << std::endl << t_all << std::endl;
        // update z
        for (auto i = 0; i < z.size(); ++i) {
            Eigen::Vector3f p(x[6 + 3 * i], x[6 + 3 * i + 1], x[6 + 3 * i + 2]);
            z.point(i) = Point(p[0], p[1], p[2]);
        }
        CGAL::draw(z);
    }
    Transform transform(R_all(0, 0), R_all(0, 1), R_all(0, 2), t_all[0],
        				R_all(1, 0), R_all(1, 1), R_all(1, 2), t_all[1],
        				R_all(2, 0), R_all(2, 1), R_all(2, 2), t_all[2], 1);
    return transform;
}


int main(int argc, char* argv[]) {

    const std::string filename = argv[1];
    PointSet sp1;
    if (!CGAL::IO::read_PLY(filename, sp1)) {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }
    // copy the point set
    PointSet sp2 = sp1;

    // define a transformation
    const float alpha = M_PI / 6;
    const float sin_alpha = sin(alpha);
    const float cos_alpha = cos(alpha);
    const Transform transform(cos_alpha, -sin_alpha, 0, 0,
                              sin_alpha,  cos_alpha, 0, 0,
                              0,          0,         1, 0, 1);
    std::cout << std::fixed << std::setprecision(4) << "transform = " << std::endl << transform << std::endl;
    // transform the points
    for (auto i = 0; i < sp2.size(); ++i) {
        Point& p = sp2.point(i);
        p = transform.transform(p);
    }
    //auto g = std::default_random_engine(42);
    //for (auto i = sp2.size() - 1; i > 0; --i) {
    //    std::uniform_int_distribution<decltype(i)> d(0, i);
    //    std::swap(sp2.point(i), sp2.point(d(g)));
    //}
    CGAL::draw(sp1);
    CGAL::draw(sp2);
    // register the point sets
    auto result = register_point_sets(sp2, sp1);
    return EXIT_SUCCESS;
}