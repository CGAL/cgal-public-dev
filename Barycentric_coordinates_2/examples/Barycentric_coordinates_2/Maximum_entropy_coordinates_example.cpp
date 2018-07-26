#include <CGAL/Barycentric_coordinates_2/Maximum_entropy_2.h>
#include <CGAL/Barycentric_coordinates_2/Maximum_entropy_2/Maximum_entropy_solver.h>
#include <CGAL/Barycentric_coordinates_2/Maximum_entropy_2/Maximum_entropy_prior_function_type_one.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>
#include <CGAL/property_map.h>

// Some convenient typedefs.
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef std::vector<Scalar> Scalar_vector;
typedef std::vector<Point>  Point_vector;

typedef std::pair<Point, bool> Point_with_property;
typedef CGAL::First_of_pair_property_map<Point_with_property> Point_map;
typedef std::vector<Point_with_property> Input_range;


typedef std::back_insert_iterator<Scalar_vector> Vector_insert_iterator;
typedef boost::optional<Vector_insert_iterator> Output_type;

typedef CGAL::Barycentric_coordinates::Maximum_entropy_newton_solver<Kernel> MEC_newton_solver;
typedef CGAL::Barycentric_coordinates::Maximum_entropy_prior_function_type_one<Kernel> MEC1_prior;

typedef CGAL::Barycentric_coordinates::Maximum_entropy_2<Kernel, MEC1_prior, MEC_newton_solver> Maximum_entropy;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Maximum_entropy, Input_range, Point_map, Kernel> Maximum_entropy_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    // Construct a concave polygon.
    const int number_of_vertices = 6;
    Point_vector vertices(number_of_vertices);

    Input_range point_range(number_of_vertices);

    vertices[0]=Point(0.0,0.0); vertices[1]=Point(-1.0,2.0); vertices[2]=Point(-2.0,-2.0);
    vertices[3]=Point(0.0,-0.5); vertices[4]=Point(3.0,-3.0); vertices[5]=Point(1.0,2.0);

    for(size_t i = 0; i < number_of_vertices; ++i)
    {
        point_range[i]=Point_with_property(vertices[i],false);
    }

    // Create an std::vector to store coordinates.
    Scalar_vector coordinates;

    // Instantiate the class with mean value coordinates for the polygon defined above.
    Maximum_entropy_coordinates maximum_entropy_coordinates(point_range, Point_map());

    maximum_entropy_coordinates.print_information();


    // Instantiate some interior points in the polygon.
    const int number_of_interior_points = 7;
    const Point interior_points[] = { Point(-1.0f , 1.0f ), Point(-1.5f, -1.5f ), Point(0.0f , -0.25f), Point(1.0f , 1.5f),
                                   Point(0.75f , -0.25f), Point(2.0f, -2.0f ), Point(2.75f, -2.75f )
                                 };

    const CGAL::Barycentric_coordinates::Type_of_algorithm type_of_algorithm = CGAL::Barycentric_coordinates::FAST;

    // We also speed up the computation by using the parameter query_point_location = CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE.
    const CGAL::Barycentric_coordinates::Query_point_location query_point_location = CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE;

    for(int i = 0; i < number_of_interior_points; ++i) {
        const Output_type result = maximum_entropy_coordinates.compute(interior_points[i], std::back_inserter(coordinates), query_point_location, type_of_algorithm);

        // Output the coordinates for each point.
        const string status = (result ? "SUCCESS." : "FAILURE.");
        cout << endl << "For the point " << i + 1 << " status of the computation: " << status << endl;

        for(int j = 0; j < number_of_vertices; ++j)
            cout << "Coordinate " << j + 1 << " = " << coordinates[i * number_of_vertices + j] << endl;
    }

    return EXIT_SUCCESS;
}
