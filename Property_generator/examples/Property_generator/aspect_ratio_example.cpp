#include <iostream>
#include <iterator>
#include <CGAL/iterator.h>
#include <CGAL/algorithm.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/properties/triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Creator_uniform_2<double, Kernel::Point_2>    Creator;

using CGAL::Properties::Triangulation_2::make_aspect_ratio;

int main()
{
  // Random triangulation.
  CGAL::Delaunay_triangulation_2<Kernel> dt;
  CGAL::Random_points_in_square_2<Kernel::Point_2, Creator> g(1);
  CGAL::cpp11::copy_n(g, 10000, std::back_inserter(dt));

  // Max aspect ratio.
  std::cout <<  CGAL::max_result(
				 CGAL::make_no_deref_iterator(dt.finite_faces_begin()),
				 CGAL::make_no_deref_iterator(dt.finite_faces_end()),
				 make_aspect_ratio(dt, CGAL::No_finite_test_tag())
				 )
            << std::endl;
}
