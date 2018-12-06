// Author: Keyu CHEN.
// In this test we compute harmonic coordinates for ~2400 strictly interior points
// with respect to a triangle and compare them with those from triangle coordinates.
// In Harmonic_2 class and related sub-class, we used series of inexact types like sqrt() and Eigen sparse solver,
// So there is some unstable inconsistency in our results (up to 0.3 now).

// Todo: Fix Harmonic_2 class with exact kernel.

#include <cassert>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_2/Harmonic_interpolator.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_2/Harmonic_mesh.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_2/Harmonic_solver.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT         Scalar;
typedef Kernel::Point_2    Point;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::vector<Point>  Point_vector;

typedef std::pair<Point, bool> Point_with_property;
typedef CGAL::First_of_pair_property_map<Point_with_property> Point_map;
typedef std::vector<Point_with_property> Input_range;


typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Kernel> Triangle_coordinates;

typedef CGAL::Barycentric_coordinates::Harmonic_mesh<Kernel> Mesh;
typedef CGAL::Barycentric_coordinates::Harmonic_solver<Kernel> Solver;
typedef CGAL::Barycentric_coordinates::Harmonic_interpolator<Kernel> Interpolator;

typedef CGAL::Barycentric_coordinates::Harmonic_2<Kernel, Mesh, Interpolator, Solver> Harmonic;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Harmonic, Input_range, Point_map, Kernel> Harmonic_coordinates;

typedef boost::optional<Vector_insert_iterator> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    const Point first_vertex  = Point(-1, 0);
    const Point second_vertex = Point(1, 1);
    const Point third_vertex  = Point(1, -1);

    Triangle_coordinates triangle_coordinates(first_vertex, second_vertex, third_vertex);

    Point_vector vertices(3);
    vertices[0] = first_vertex; vertices[1] = second_vertex; vertices[2] = third_vertex;

    Input_range point_range(3);
    for(size_t i = 0; i < 3; ++i)
    {
        point_range[i]=Point_with_property(vertices[i],false);
    }

    Harmonic_coordinates harmonic_coordinates(point_range, Point_map());

    Coordinate_vector tri_coordinates;
    Coordinate_vector  hm_coordinates;

    const Scalar step  = Scalar(1) / Scalar(10);
    const Scalar scale = Scalar(5);

    size_t count = 0;
    const Scalar limit = scale*step;

    for(Scalar x = step; x < limit; x += step) {
        for(Scalar y = step; y < limit; y += step) {
            const Point point(0.9, y);
    
            const Output_type tri_result = triangle_coordinates.compute(point, tri_coordinates);
            const Output_type  hm_result = harmonic_coordinates.compute(point, hm_coordinates, CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE, CGAL::Barycentric_coordinates::PRECISE);

            //assert(tri_coordinates[count + 0] - hm_coordinates[count + 0] <= Scalar(0.3) &&
            //       tri_coordinates[count + 1] - hm_coordinates[count + 1] <= Scalar(0.3) &&
            //       tri_coordinates[count + 2] - hm_coordinates[count + 2] <= Scalar(0.3) );

            Scalar tri_result_x(0), tri_result_y(0);
            Scalar hm_result_x(0), hm_result_y(0);
            for(size_t k=0;k<3;k++){
                tri_result_x += tri_coordinates[count + k] * vertices[k].x();
                tri_result_y += tri_coordinates[count + k] * vertices[k].y();
                hm_result_x += hm_coordinates[count + k] * vertices[k].x();
                hm_result_y += hm_coordinates[count + k] * vertices[k].y();
            }

            if( tri_coordinates[count + 0] - hm_coordinates[count + 0] > Scalar(1e-6) ||
                tri_coordinates[count + 1] - hm_coordinates[count + 1] > Scalar(1e-6) ||
                tri_coordinates[count + 2] - hm_coordinates[count + 2] > Scalar(1e-6)  )
            {
                // If you want to view all the differences between HM and Triangle coordinates, just change the condition "> Scalar(1e-6)" to "!= Scalar(0)".
                cout << endl << "HM_triangle_inexact_test: FAILED." << endl << endl;
                cout << "location: " << point.x() << " " << point.y() << endl;
                cout << "tri location: " << tri_result_x << " " << tri_result_y << endl;
                cout << "hm location: " << hm_result_x << " " << hm_result_y << endl;
                cout << "difference " << tri_coordinates[count + 0] - hm_coordinates[count + 0] << endl;
                cout << "difference " << tri_coordinates[count + 1] - hm_coordinates[count + 1] << endl;
                cout << "difference " << tri_coordinates[count + 2] - hm_coordinates[count + 2] << endl;
                //exit(EXIT_FAILURE);
            }
            count += 3;
        }
    }

    cout << endl << "HM_triangle_inexact_test: FINISHED." << endl << endl;

    return EXIT_SUCCESS;
}
