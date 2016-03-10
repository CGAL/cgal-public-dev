#include <vector>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <CGAL/iterator.h>
#include <CGAL/algorithm.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/properties/triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_2                                      Point;
typedef CGAL::Delaunay_triangulation_2<Kernel>               Delaunay;
typedef CGAL::Creator_uniform_2<double, Point>               Creator;



// Alias the namespaces to save typing.
namespace Properties = CGAL::Properties::Triangulation_2;

int main()
{
    // Number of points to generate.
    unsigned           n = 1<<10;
    Delaunay           dt;
    std::vector<Point> points;

    // Random triangulation.
    CGAL::Random_points_in_square_2<Point,Creator> g(0.5);
    CGAL::cpp11::copy_n(g, n, std::back_inserter(points));
    dt.insert(points.begin(), points.end());

    // The property functions take handle types, so we wrap the iterators.
    CGAL::No_deref_iterator<Delaunay::Finite_faces_iterator> begin, end;

    begin = make_no_deref_iterator(dt.finite_faces_begin());
    end   = make_no_deref_iterator(dt.finite_faces_end());

    // Initialise the functors that we shall use.
    Properties::Area<Delaunay>         area(dt);
    Properties::Circumradius<Delaunay> circumradius(dt);
    Properties::Aspect_ratio<Delaunay> aspect_ratio(dt);
    Properties::Max_angle<Delaunay>    max_angle(dt);

    // Display some statistics about the triangulation.
    std::cout
        << "-- Information about the triangulation --" 
        << std::endl << std::left << std::setw(50)
        
        << "Mean face area"
        << mean_result(begin, end, area)
        << std::endl << std::left << std::setw(50)

        << "Largest face area"
        << max_result(begin, end, area)
        << std::endl << std::left << std::setw(50)

        << "Mean circumradius"
        << mean_result(begin, end, circumradius)
        << std::endl << std::left << std::setw(50)

        << "Maximum aspect ratio"
        << max_result(begin, end, aspect_ratio)
        << std::endl << std::left << std::setw(50)

	<< "Correlation between area and aspect ratio"
        << CGAL::pearson(begin, end, area, aspect_ratio)
        << std::endl << std::left << std::setw(50)

        << "Number of angles larger than three"
        <<  CGAL::count_result_in_interval(begin, end, max_angle, 3,4)
        << std::endl;
}


