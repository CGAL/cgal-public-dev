#include <iostream>
#include <iterator>
#include <CGAL/algorithm.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/properties/triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Delaunay_triangulation_2<Kernel>              Delaunay;
typedef CGAL::Creator_uniform_2<double,Kernel::Point_2>     Creator;

// The property functors for Triangulation_2.
namespace Properties = CGAL::Properties::Triangulation_2;

int main()
{ 
    Delaunay dt;
    
    // Random triangulation.
    CGAL::Random_points_in_square_2<Kernel::Point_2,Creator> g(0.5);
    CGAL::cpp11::copy_n( g, 100, std::back_inserter(dt) );
    
    // Create a functor for computing aspect ratio of faces.
    Properties::Aspect_ratio<Delaunay> aspect_ratio(dt);

    // A face from the triangulation.
    Delaunay::Face_handle face = dt.finite_faces_begin();
    
    // Compute and display the aspect ratio of a face.
    std::cout << aspect_ratio(face) << std::endl;
}