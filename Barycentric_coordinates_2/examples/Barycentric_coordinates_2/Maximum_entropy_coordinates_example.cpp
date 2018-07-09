#include <CGAL/Barycentric_coordinates_2/Maximum_entropy_2.h>
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

typedef CGAL::Barycentric_coordinates::Maximum_entropy_newton_solver<Kernel, Point_with_property, Point_map> MEC_newton_solver;
typedef CGAL::Barycentric_coordinates::Maximum_entropy_prior_function_type_one_2<Kernel, Point_with_property, Point_map> MEC1_prior;

typedef CGAL::Barycentric_coordinates::Maximum_entropy_2<Kernel, MEC1_prior, MEC_newton_solver, Point_with_property, Point_map> Maximum_entropy;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Maximum_entropy, Kernel, Point_with_property, Point_map> Maximum_entropy_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    // Construct a star-shaped polygon.
    const int number_of_vertices = 3;
    Point_vector vertices(number_of_vertices);

    Input_range point_range(number_of_vertices);


    //vertices[0] = Point(0.0, 0.0); vertices[1] = Point(0.1, -0.8); vertices[2] = Point(0.3, 0.0); vertices[3] = Point(0.6, -0.5); vertices[4]  = Point(0.6 , 0.1);
    //vertices[5] = Point(1.1, 0.6); vertices[6] = Point(0.3,  0.2); vertices[7] = Point(0.1, 0.8); vertices[8] = Point(0.1,  0.2); vertices[9] = Point(-0.7, 0.0);
    vertices[0]=Point(0.0,0.0); vertices[1]=Point(2.0,0.5);
    vertices[2]=Point(1.0,2.0);

    point_range[0]=Point_with_property(vertices[0],false);
    point_range[1]=Point_with_property(vertices[1],false);
    point_range[2]=Point_with_property(vertices[2],false);

    Point_map point_map;


    // Create an std::vector to store coordinates.
    Scalar_vector coordinates;

    //MEC1_prior prior(vertices);
    //MEC_newton_solver newton_solver(vertices);
    //MEC_newton_solver mec_newton_solver(vertices);
    //MEC1_prior mec1_prior(vertices);

    // Instantiate the class with mean value coordinates for the polygon defined above.
    Maximum_entropy_coordinates maximum_entropy_coordinates(point_range.begin(), point_range.end(), Point_map());
    //Maximum_entropy_coordinates maximum_entropy_coordinates(vertices.begin(), vertices.end());

    // Print some information about the polygon and coordinates.
    maximum_entropy_coordinates.print_information();


    // Instantiate some interior points in the polygon.
    //const int number_of_interior_points = 2;
    //const Point interior_points[] = { Point(0.12, -0.45), Point(0.55, -0.3), Point(0.9 , 0.45),
    //                                  Point(0.15,  0.35), Point(-0.4, 0.04), Point(0.11, 0.11),
    //                                  Point(0.28,  0.12), // the only point in the kernel of the star shaped polygon
    //                                  Point(0.55,  0.11)
    //                                };
    //const Point interior_points[] = { Point(0.5, 0.5), Point(0.01, 0.01) };
    const int number_of_interior_points = 9;
    const Point interior_points[] = { Point(0.5f , 0.5f ), Point(1.0f, 0.5f ), Point(1.0f , 0.75f), Point(1.0f , 1.0f),                     // interior query points
                                   Point(1.0f , 1.25f), Point(1.0f, 1.5f ), Point(0.75f, 1.0f ), Point(1.25f, 1.0f), Point(1.5f, 0.75f),
                                   //Point(1.0f , 0.25f), Point(0.5f, 1.0f ), Point(1.5f , 1.25f), Point(1.0f , 2.0f), Point(2.0f, 0.5f ), // boundary query points
                                   //Point(0.25f, 1.0f ), Point(0.5f, 1.75f), Point(1.5f , 1.75f), Point(1.75f, 1.5f)                      // exterior query points
                                 };

    // Compute mean value coordinates for all the defined interior points.

    // We speed up the computation using the O(n) algorithm called with the parameter CGAL::Barycentric_coordinates::FAST.
    // The default one is CGAL::Barycentric_coordinates::PRECISE.
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

    // If we need only the unnormalized weights for some point (lets take the last one), we can compute them as follows.

    /*
    // Instantiate an std::vector to store weights.
    Scalar_vector weights;

    // Compute mean value weights.
    const int last_point_index = number_of_interior_points - 1;
    const Output_type result = maximum_entropy_coordinates.compute_weights(interior_points[last_point_index], std::back_inserter(weights));

    // Compute their sum.
    Scalar mv_denominator = Scalar(0);
    for(int j = 0; j < number_of_vertices; ++j) mv_denominator += weights[j];

    // Invert this sum.
    const Scalar mv_inverted_denominator = Scalar(1) / mv_denominator;

    // Output the mean value weights.
    const string status = (result ? "SUCCESS." : "FAILURE.");
    cout << endl << "Status of the weights' computation for the point " << last_point_index + 1 << ": " << status << endl;

    for(int j = 0; j < number_of_vertices; ++j)
        cout << "Weight " << j + 1 << " = " << weights[j] << endl;

    // Now, if we normalize the weights, we recover values of the mean value coordinates for the last point computed earlier.
    cout << endl << "After normalization, for the point " << last_point_index + 1 << " mean value coordinates are " << endl;
    for(int j = 0; j < number_of_vertices; ++j)
        cout << "Coordinate " << j + 1 << " = " << weights[j] * mv_inverted_denominator << endl;
    cout << endl;
    */


    return EXIT_SUCCESS;
}
