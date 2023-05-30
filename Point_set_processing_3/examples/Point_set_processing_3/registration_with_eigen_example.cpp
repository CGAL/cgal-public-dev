#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/draw_point_set_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>
#include <Eigen/Dense>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Point_set_3<Point>                             PointSet;
typedef CGAL::Aff_transformation_3<Kernel>                   Transform;

std::pair<bool, Transform> register_point_sets(const PointSet& sp1, const PointSet& sp2, float max_iterations) {
    // initialize the transformation
    Transform T;
    for (auto ti = 0; ti < max_iterations; ++ti) {
        // reserve memory for matrix A and vector b
        Eigen::MatrixXf A(sp2.size(), 6);
        Eigen::VectorXf b(sp2.size());
        // fill matrix A and vector b
        for (auto i = 0; i < sp2.size(); ++i) {
            Point p = T(sp1.point(i));
            Eigen::Vector3f s(p[0], p[1], p[2]);
            Eigen::Vector3f d(sp2.point(i)[0], sp2.point(i)[1], sp2.point(i)[2]);
            Eigen::Vector3f n(sp2.normal(i)[0], sp2.normal(i)[1], sp2.normal(i)[2]);
            auto a = s.cross(n);
            A.row(i) << a[0], a[1], a[2], n[0], n[1], n[2];
            b[i] = n.dot(d - s);
        }
        // solve the linear system
        Eigen::VectorXf x = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
        // transform the points
        auto roll = Eigen::AngleAxisf(x[0], Eigen::Vector3f::UnitX());
        auto pitch = Eigen::AngleAxisf(x[1], Eigen::Vector3f::UnitY());
        auto yaw = Eigen::AngleAxisf(x[2], Eigen::Vector3f::UnitZ());
        Eigen::Matrix3f R = (roll * pitch * yaw).matrix();
        Eigen::Vector3f t(x[3], x[4], x[5]);
        T = Transform(R(0, 0), R(0, 1), R(0, 2), t[0],
                      R(1, 0), R(1, 1), R(1, 2), t[1],
                      R(2, 0), R(2, 1), R(2, 2), t[2], 1) * T;
        // return if the transformation is small enough
        std::cout << ti << ": " << x.norm() << std::endl;
        if (x.norm() < 1e-4) {
			return std::pair<bool, Transform>(true, T);
		}
    }
    return std::pair<bool, Transform>(false, T);
}


int main(int argc, char* argv[]) {
    const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("pointsets/hippo1.ply");
    PointSet sp1;
    if (!CGAL::IO::read_PLY(filename, sp1)) {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }
    // copy the point set in point set 2
    PointSet sp2 = sp1;
    // define a transformation
    const float alpha = M_PI / 6;
    const float sin_alpha = sin(alpha);
    const float cos_alpha = cos(alpha);
    const CGAL::Aff_transformation_3<Kernel> transform(cos_alpha, -sin_alpha, 0, 2,
                                                       sin_alpha, cos_alpha, 0, 2,
                                                       0, 0, 1, 2, 1);
    std::cout << std::fixed << std::setprecision(4) << "transform = " << std::endl << transform << std::endl;
    // transform the points of point set 2
    for (auto i = 0; i < sp2.size(); ++i) {
        Point& p = sp2.point(i);
        p = transform.transform(p);
    }
    // register the point sets
    auto result = register_point_sets(sp1, sp2, 100);
    auto success = result.first;
    auto T = result.second;
    if (success) {
		std::cout << "registration succeeded" << std::endl;
        std::cout << std::fixed << std::setprecision(4) << "transform = " << std::endl << result.second << std::endl;
    }
    else {
		std::cout << "registration failed" << std::endl;
	}
    return EXIT_SUCCESS;
}
