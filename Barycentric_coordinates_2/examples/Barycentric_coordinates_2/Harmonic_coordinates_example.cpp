#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

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

typedef CGAL::Barycentric_coordinates::Harmonic_mesh_2<Kernel, Point_with_property, Point_map> Mesh;
typedef CGAL::Barycentric_coordinates::Harmonic_solver_2<Kernel, Point_with_property, Point_map> Solver;
typedef CGAL::Barycentric_coordinates::Harmonic_interpolator_2<Kernel> Interpolator;

typedef CGAL::Barycentric_coordinates::Harmonic_2<Kernel, Mesh, Interpolator, Solver, Point_with_property, Point_map> Harmonic;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Harmonic, Kernel, Point_with_property, Point_map> Harmonic_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    // Construct a unit square.
    const int number_of_vertices = 4;
    Point_vector vertices(number_of_vertices);

    Input_range point_range(number_of_vertices);

    vertices[0] = Point(0, 0); vertices[1] = Point(1, 0); vertices[2] = Point(1, 1); vertices[3] = Point(0, 1);
    //vertices[0] = Point(0, 0); vertices[1] = Point(1, 1); vertices[2] = Point(2, 0); vertices[3] = Point(1, -0.5);
    //vertices[4] = Point(2, -1); vertices[5] = Point(0, -2); vertices[6] = Point(-2, -1); vertices[7] = Point(0, -1);
    point_range[0]=Point_with_property(vertices[0],false);
    point_range[1]=Point_with_property(vertices[1],false);
    point_range[2]=Point_with_property(vertices[2],false);
    point_range[3]=Point_with_property(vertices[3],false);

    Point_map point_map;

    // Create an std::vector to store coordinates.
    Scalar_vector coordinates;

    // Instantiate the class with Harmonic coordinates for the unit square defined above.
    Harmonic_coordinates harmonic_coordinates(point_range.begin(), point_range.end(), Point_map());

    int number_of_interior_points = 9;
    const Point interior_points[] = { Point(0.5f , 0.5f ), Point(0.9f, 0.5f ), Point(0.9f , 0.75f), Point(0.9f , 0.9f),
                                   Point(0.1f , 0.3f), Point(0.1f, 0.1f ), Point(0.75f, 0.9f ), Point(0.25f, 0.9f), Point(0.5f, 0.75f)
                                 };
    //const Point interior_points[] = {Point(-1, -1.1), Point(1, 0.5), Point(1, 0), Point(0.5, 0), Point(0.5, -0.2), Point(0.6, -1), Point(1.5, -1)};

    const CGAL::Barycentric_coordinates::Query_point_location query_point_location = CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE;

    const CGAL::Barycentric_coordinates::Type_of_algorithm type_of_algorithm = CGAL::Barycentric_coordinates::PRECISE;

    for(int i = 0; i < number_of_interior_points; ++i) {
        const Output_type result = harmonic_coordinates.compute(interior_points[i], std::back_inserter(coordinates), query_point_location, type_of_algorithm);

        // Output the coordinates for each point.
        const string status = (result ? "SUCCESS." : "FAILURE.");
        cout << endl << "For the point " << i + 1 << " status of the computation: " << status << endl;
        cout << "Location " << interior_points[i].x() << " " << interior_points[i].y() << endl;

        for(int j = 0; j < number_of_vertices; ++j)
            cout << "Coordinate " << j + 1 << " = " << coordinates[i * number_of_vertices + j] << endl;
    }

    return EXIT_SUCCESS;
}
