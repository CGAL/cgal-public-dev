/*
  This test checks that a point sample converges 
  to the correct closed form covariance matrix
  Author: Pierre Alliez
 */

#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include "include/CGAL/test_utils.h"
#include <CGAL/Default_diagonalize_traits.h>
#include <CGAL/Random.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;

typedef Kernel::Line_2 Line_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Iso_rectangle_2 Iso_rectangle_2;
typedef CGAL::Default_diagonalize_traits<typename Kernel::FT, 2> Diagonalize_traits_2;
typedef Diagonalize_traits_2::Covariance_matrix Moment;

// compute distance between two 2D moment matrices
// by summing absolute value of term-by-term differences
double distance_moment_2(Moment& ref, Moment& moment)
{
	return std::fabs(moment[0] - ref[0]) +
		   std::fabs(moment[1] - ref[1]) +
		   std::fabs(moment[2] - ref[2]);
}

void test_pca_convergence_iso_rectangle_2()
{
	Point_2 origin_2(0.0, 0.0);

	std::list<Iso_rectangle_2> iso_rectangles;
	iso_rectangles.push_back(Iso_rectangle_2(origin_2, Point_2(1.0, 1.0)));

	Moment ref =
		CGAL::order_2_moment_2(iso_rectangles.begin(), iso_rectangles.end(), origin_2,
			Kernel(), CGAL::Dimension_tag<2>(), Diagonalize_traits_2());

	std::cerr << "Iso rectangle closed form reference: ";
	std::cerr << ref[0] << " " << ref[1] << " " << ref[2] << std::endl;

	// convergence loop
	const int max_nbsamples = 1e8;
	const int init = 100;
	CGAL::Random random(0); // random generator

	double previous_distance = 0.0; // used to verify decreasing of error

	for (int nb_samples = init; nb_samples < max_nbsamples; nb_samples *= 10)
	{
		std::list<Point_2> points;
		for (int i = 0; i < nb_samples; i++)
		{
			const FT x = random.get_double(0.0, 1.0);
			const FT y = random.get_double(0.0, 1.0);
			points.push_back(Point_2(x, y));
		}

		Moment moment =
			CGAL::order_2_moment_2(points.begin(), points.end(), origin_2,
				Kernel(), CGAL::Dimension_tag<0>(), Diagonalize_traits_2());

		moment[0] /= (FT)points.size();
		moment[1] /= (FT)points.size();
		moment[2] /= (FT)points.size();

		double curr_distance = distance_moment_2(ref, moment);

		std::cerr << "Moment matrix for " << points.size() << " random points: ";
		std::cerr << moment[0] << " " << moment[1] << " " << moment[2] << 
			" distance: " << curr_distance << std::endl;

		// verify convergence = decreasing error
		if(previous_distance != 0.0) 
		{
			if (curr_distance > previous_distance)
			{
				std::cerr << "failure: increasing error" << std::endl;
				std::exit(1); // failure
			}
		}

		// update previous distance
		previous_distance = curr_distance; 
	}
}

int main()
{
  std::cerr << "Convergence tests in 2D" << std::endl;
  test_pca_convergence_iso_rectangle_2();

  return EXIT_SUCCESS;
}
