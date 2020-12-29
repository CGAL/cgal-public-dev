/*
  This test checks that the closed form covariance matrices are correct
  Author: Pierre Alliez
 */

#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include "include/CGAL/test_utils.h"
#include <CGAL/Default_diagonalize_traits.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;

typedef Kernel::Line_2 Line_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Iso_rectangle_2 Iso_rectangle_2;
typedef CGAL::Default_diagonalize_traits<typename Kernel::FT, 2> Diagonalize_traits_2;

void test_iso_rectangle_2()
{
	Point_2 origin_2(0.0, 0.0);

	std::list<Iso_rectangle_2> Iso_rectangles;
	Iso_rectangles.push_back(Iso_rectangle_2(origin_2, Point_2(1.0, 1.0)));

	Diagonalize_traits_2::Covariance_matrix moment;
		CGAL::compute_moment_2(Iso_rectangles.begin(), Iso_rectangles.end(), 
			moment, origin_2, Kernel(), CGAL::Dimension_tag<2>());

	std::cerr << "- Iso rectangle: ";
	std::cerr << moment[0] << " " << moment[1] << " " << moment[2];

	// reference = 1/3, 1/4, 1/3
	if(moment[0] != 1.0 / 3.0 || 
	   moment[1] != 0.25 ||
	   moment[2] != 1.0 / 3.0)
	{
		std::cerr << " (failure)" << std::endl;
		std::exit(1); // failure
	}
	else
		std::cerr << " (success)" << std::endl;
}

void test_segment_2()
{
	Point_2 origin_2(0.0, 0.0);

	std::list<Segment_2> segments;
	segments.push_back(Segment_2(Point_2(0.0, 1.0), Point_2(1.0, 0.0)));

	Diagonalize_traits_2::Covariance_matrix moment;
	CGAL::compute_moment_2(segments.begin(), segments.end(),
		moment, origin_2, Kernel(), CGAL::Dimension_tag<1>());

	const double SQRT2 = std::sqrt(2.0);

	std::cerr << "- Segment: ";
	std::cerr << moment[0] << " " << moment[1] << " " << moment[2] << std::endl;

	std::cerr << "  expected: ";
	std::cerr << SQRT2 /3.0 << " " << SQRT2 / 6.0 << " " << SQRT2 / 3.0 << std::endl;

	// reference = sqrt(2)/3, sqrt(2)/6, sqrt(2)/3
	const double diff = std::fabs(moment[0] - SQRT2 / 3.0) +
		std::fabs(moment[1] - SQRT2 / 6.0) +
		std::fabs(moment[2] - SQRT2 / 3.0);
	std::cerr << "  diff: " << diff;

	if (diff > 1e-6)
	{
		std::cerr << " (failure)" << std::endl;
		std::exit(1); // failure
	}
	else
		std::cerr << " (success)" << std::endl;
}

int main()
{
  std::cerr << "Closed form tests in 2D" << std::endl;
  test_iso_rectangle_2();
  test_segment_2();

  return EXIT_SUCCESS;
}
