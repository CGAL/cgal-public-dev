/*
    An example to demonstrate the use of CGAL::No_deref_iterator for providing
    iterators having a value type equivalent to a CGAL handle.
*/

#include <vector>
#include <algorithm>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/iterator.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/property_functors/triangulation_2_face_properties.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Creator_uniform_2<double, Kernel::Point_2>     Creator;
typedef CGAL::Delaunay_triangulation_2<Kernel>               Delaunay;

int main()
{
    // Random triangulation.
    CGAL::Delaunay_triangulation_2<Kernel>                     dt;
    CGAL::Random_points_in_square_2<Kernel::Point_2, Creator>  g(1);
    CGAL::cpp11::copy_n(g, 10000, std::back_inserter(dt));

    // Vector to contain face areas.
    std::vector<double> areas( dt.number_of_faces() );

    // Compute area of all faces in a triangulation and store them in a vector.
    std::transform(
        CGAL::make_no_deref_iterator( dt.finite_faces_begin() ),
        CGAL::make_no_deref_iterator( dt.finite_faces_end()   ),
        areas.begin(),
        CGAL::Properties::Triangulation_2_properties::make_area(dt)
    );
}