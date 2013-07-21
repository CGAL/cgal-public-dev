// Tester for Random_points_in_triangle_3 getting a triangle_3 as constructor
// argument

#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <CGAL/double.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <CGAL/point_generators_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel 	K;
typedef K::Point_3 						Point_3;
typedef K::Vector_3						Vector_3;
typedef K::Triangle_3						Triangle_3;

using namespace std;

int main()
{
	CGAL::Random rand;
	int D = 3;
	vector<Point_3> output;
	double **coord = new double*[D];
	for (int i = 0; i < D; i++) {
		coord[i] = new double[D];
		for(int j = 0; j < D; j++) {
			coord[i][j] = rand.get_double(50, 60);
		}
	}

	cout << "Coordinates of the edges of the " << D << "-Simplex are:\n";
	for (int i = 0; i < D; i++) {
		for(int j = 0; j < D; j++) {
			cout << coord[i][j] << " ";
		}
		cout << '\n';
	}

	Point_3 *pts = new Point_3[3];
	for (int i = 0; i < 3; i++) {
		pts[i] = Point_3(coord[i][0], coord[i][1], coord[i][2]);
	}

	Triangle_3 tr(pts[0], pts[1], pts[2]);

	CGAL::Random_points_in_triangle_3 <Point_3> g(tr);
	CGAL::cpp11::copy_n( g, 2, std::back_inserter(output));

	cout << "The generated points are: " << '\n';
	for (int i = 0; i < 2; i++) {
		cout << output[i].x() << " " << output[i].y() << " " <<
			output[i].z() << '\n';
	}

	return 0;
}

