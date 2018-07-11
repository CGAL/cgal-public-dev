// Author: Keyu CHEN.
// In this test we compute maximum entropy coordinates for the same interior points
// as we used in gbc package (https://github.com/danston/gbc/blob/master/2d/main.cpp).
// The input polygons are the same and the interior query points (mesh-based result in gbc
// package) are loaded from "./MEC_test/mesh_vertices_location_gbc.txt".
// The computing result should be as same as gbc's, which is stored in "./MEC_test/mesh_vertices_coordinates_gbc.txt".

// Note that we can only use inexact kernel now, there are little inconsistency and we set up a max tolerance 1e-5.

// Todo: Fix Maximum_entropy_2 class with exact kernel.

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Maximum_entropy_2.h>
#include <CGAL/Barycentric_coordinates_2/Maximum_entropy_2/Maximum_entropy_solver.h>
#include <CGAL/Barycentric_coordinates_2/Maximum_entropy_2/Maximum_entropy_prior_function_type_one.h>
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

typedef CGAL::Barycentric_coordinates::Maximum_entropy_newton_solver<Kernel> MEC_newton_solver;
typedef CGAL::Barycentric_coordinates::Maximum_entropy_prior_function_type_one<Kernel> MEC1_prior;

typedef CGAL::Barycentric_coordinates::Maximum_entropy_2<Kernel, MEC1_prior, MEC_newton_solver> Maximum_entropy;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Maximum_entropy, Input_range, Point_map, Kernel> Maximum_entropy_coordinates;

typedef boost::optional<Vector_insert_iterator> Output_type;

using std::cout; using std::endl; using std::string;
using std::ifstream; using std::ofstream;

using namespace std;

int main()
{

    Point_vector vertices(4);
    Input_range point_range(4);

    vertices[0] = Point(0.1, 0.1);
    vertices[1] = Point(1.0, 0.0);
    vertices[2] = Point(0.9, 0.9);
    vertices[3] = Point(0.2, 1.0);

    point_range[0]=Point_with_property(vertices[0],false);
    point_range[1]=Point_with_property(vertices[1],false);
    point_range[2]=Point_with_property(vertices[2],false);
    point_range[3]=Point_with_property(vertices[3],false);


    Maximum_entropy_coordinates maximum_entropy_coordinates(point_range, Point_map());

    Coordinate_vector  me_coordinates;

    ifstream gbc_query_points("../MEC_test/mesh_vertices_location_gbc.txt");
    ifstream gbc_coordinates("../MEC_test/mesh_vertices_coordinates_gbc.txt");

    //gbc_query_points.open("../MEC_test/mesh_vertices_location_gbc.txt");
    //gbc_coordinates.open("../MEC_test/mesh_vertices_coordinates_gbc.txt");

    int count = 0;

    string query_point_line;
    string coordinates_line;
    while( getline(gbc_query_points, query_point_line) && getline(gbc_coordinates, coordinates_line)) {

        std::stringstream query_point_in(query_point_line);
        std::stringstream coordinates_in(coordinates_line);

        double x, y;
        query_point_in >> x;
        query_point_in >> y;

        Point point(x, y);

        const Output_type  dh_result = maximum_entropy_coordinates.compute(point, me_coordinates);
        //assert(tri_coordinates[count + 0] - me_coordinates[count + 0] == Scalar(0) &&
        //       tri_coordinates[count + 1] - me_coordinates[count + 1] == Scalar(0) &&
        //       tri_coordinates[count + 2] - me_coordinates[count + 2] == Scalar(0) );

        double c1, c2, c3, c4;
        coordinates_in >> c1;
        coordinates_in >> c2;
        coordinates_in >> c3;
        coordinates_in >> c4;

        if( Scalar(c1) - me_coordinates[count + 0] > Scalar(1e-5) ||
            Scalar(c2) - me_coordinates[count + 1] > Scalar(1e-5) ||
            Scalar(c3) - me_coordinates[count + 2] > Scalar(1e-5) ||
            Scalar(c4) - me_coordinates[count + 3] > Scalar(1e-5) )
        {
            cout << endl << "MEC_compare_with_gbc_test: FAILED." << endl << endl;
            exit(EXIT_FAILURE);
        }
        //if( tri_coordinates[count + 0] - me_coordinates[count + 0] != Scalar(0) ||
        //    tri_coordinates[count + 1] - me_coordinates[count + 1] != Scalar(0) ||
        //    tri_coordinates[count + 2] - me_coordinates[count + 2] != Scalar(0)  )
        //{
        //    cout << endl << "MEC_triangle_test: FAILED." << endl << endl;
        //    exit(EXIT_FAILURE);
        //}
        count += 4;
        }

    cout << endl << "MEC_compare_with_gbc_test: PASSED." << endl << endl;

    return EXIT_SUCCESS;
}
