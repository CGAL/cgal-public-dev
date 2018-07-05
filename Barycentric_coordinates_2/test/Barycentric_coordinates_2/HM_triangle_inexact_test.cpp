// Author: Keyu CHEN.
// In this test we compute harmonic coordinates for ~2400 strictly interior points
// with respect to a triangle and compare them with those from triangle coordinates.
// In Harmonic_2 class and related sub-class, we used series of inexact types like sqrt() and Eigen sparse solver,
// So there is some unstable inconsistency in our results (up to 0.3 nuw).

// Todo: Fix Harmonic_2 class with exact kernel.

#include <cassert>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT         Scalar;
typedef Kernel::Point_2    Point;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::vector<Point>  Point_vector;

typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Kernel> Triangle_coordinates;

typedef CGAL::Barycentric_coordinates::Harmonic_mesh_2<Kernel> Mesh;
typedef CGAL::Barycentric_coordinates::Harmonic_solver_2<Kernel> Solver;
typedef CGAL::Barycentric_coordinates::Harmonic_interpolator_2<Kernel> Interpolator;

typedef CGAL::Barycentric_coordinates::Harmonic_2<Kernel, Mesh, Interpolator, Solver> Harmonic;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Harmonic, Kernel> Harmonic_coordinates;

typedef boost::optional<Vector_insert_iterator> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    const Point first_vertex  = Point(0.0, 0.0);
    const Point second_vertex = Point(1.0, 0.0);
    const Point third_vertex  = Point(0.0, 1.0);

    Triangle_coordinates triangle_coordinates(first_vertex, second_vertex, third_vertex);

    Point_vector vertices(3);
    vertices[0] = first_vertex; vertices[1] = second_vertex; vertices[2] = third_vertex;

    Harmonic_coordinates harmonic_coordinates(vertices.begin(), vertices.end());

    Coordinate_vector tri_coordinates;
    Coordinate_vector  hm_coordinates;

    const Scalar step  = Scalar(1) / Scalar(10);
    const Scalar scale = Scalar(5);

    int count = 0;
    const Scalar limit = scale*step;

    for(Scalar x = step; x < limit; x += step) {
        for(Scalar y = step; y < limit; y += step) {
            const Point point(x, y);

            const Output_type tri_result = triangle_coordinates(point, tri_coordinates);
            const Output_type  hm_result = harmonic_coordinates(point, hm_coordinates);

            //assert(tri_coordinates[count + 0] - hm_coordinates[count + 0] <= Scalar(0.3) &&
            //       tri_coordinates[count + 1] - hm_coordinates[count + 1] <= Scalar(0.3) &&
            //       tri_coordinates[count + 2] - hm_coordinates[count + 2] <= Scalar(0.3) );

            if( tri_coordinates[count + 0] - hm_coordinates[count + 0] > Scalar(0.3) ||
                tri_coordinates[count + 1] - hm_coordinates[count + 1] > Scalar(0.3) ||
                tri_coordinates[count + 2] - hm_coordinates[count + 2] > Scalar(0.3)  )
            {
                cout << endl << "HM_triangle_inexact_test: FAILED." << endl << endl;
                exit(EXIT_FAILURE);
            }
            count += 3;
        }
    }

    cout << endl << "HM_triangle_inexact_test: PASSED." << endl << endl;

    return EXIT_SUCCESS;
}
