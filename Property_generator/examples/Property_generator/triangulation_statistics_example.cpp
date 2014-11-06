/*
    An example to show the generation of various statistics of a randomly
    generated Delaunay triangulation.
*/

#include <vector>
#include <iomanip>
#include <iostream>
#include <CGAL/property_generator_algorithms.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/properties/triangulation_2.h>


// #include <CGAL/properties/meta_functions.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_2                                      Point;
typedef CGAL::Delaunay_triangulation_2<Kernel>               Delaunay;
typedef CGAL::Creator_uniform_2<double, Point>               Creator;

// Alias the namespaces to save typing.
namespace Properties    = CGAL::Properties::Triangulation_2;

int foo(Delaunay::Face_handle)
{
    return 1;
}

int main()
{
    // Number of points to generate.
    unsigned           n = 1<<20;
    Delaunay           dt;
    std::vector<Point> points;

    // Random triangulation.
    CGAL::Random_points_in_square_2<Point,Creator> g(0.5);
    CGAL::cpp11::copy_n(g, n, std::back_inserter(points));
    dt.insert(points.begin(), points.end());

    // Store iterator range as temporaries for convenience.
    Delaunay::Finite_faces_iterator begin = dt.finite_faces_begin();
    Delaunay::Finite_faces_iterator end   = dt.finite_faces_end();

    // Initialise the functors that we shall use.
    Properties::Area<Delaunay>         area(dt);
    Properties::Circumradius<Delaunay> circumradius(dt);
    Properties::Aspect_ratio<Delaunay> aspect_ratio(dt);
    Properties::Max_angle<Delaunay>    max_angle(dt);


    std::vector<Delaunay::Face_handle> temp;
    CGAL::max_result(temp.begin(), temp.end(), foo);
    CGAL::max_result(begin, end, foo);

    // std::cout <<
    //     // int (CGAL::internal::do_dereference<decltype(aspect_ratio), std::vector<Delaunay::Face_handle>::iterator>::value)
    //     int (CGAL::internal::is_callable_with<decltype(aspect_ratio), decltype(begin)>::value )
    //     // int( CGAL::internal::is_callable_with)
    //     << std::endl;


    
    // // Display some statistics about the triangulation.
    // std::cout
    //     << "-- Information about the triangulation --"

    //     << "Mean face area"
    //     << STL_extension::mean_result(begin, end, area)
    //     << std::endl << std::left << std::setw(50)

    //     << "Largest face area"
    //     << STL_extension::max_result(begin, end, area)
    //     << std::endl << std::left << std::setw(50)

    //     << "Mean circumradius"
    //     << STL_extension::mean_result(begin, end, circumradius)
    //     << std::endl << std::left << std::setw(50)

    //     << "Maximum aspect ratio"
    //     << STL_extension::max_result(begin, end, aspect_ratio)
    //     << std::endl << std::left << std::setw(50)

    //     << "Correlation between area and aspect ratio"
    //     << STL_extension::pearson(begin, end, area, aspect_ratio)
    //     << std::endl << std::left << std::setw(50)

    //     << "Number of angles larger than 100 degrees"
    //     <<  STL_extension::count_result_in_interval(begin, end, max_angle, 100,360)
    // << std::endl;
}


