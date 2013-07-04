#include <iostream>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <iterator>
#include <CGAL/enum.h>
#include "generate_barycoords.h"
#include "barycoords_to_cartesian.h"
#include <CGAL/Timer.h>

#define D 2 // the d dimension
#define total 10 // number of points that will be generated
#define NUMBER_OF_TESTS 1

#define VERBOSE

using namespace std;

void benchmark() {
	CGAL::Random rand;
	CGAL::Timer t;

	vector<double> output;
	double **coord = new double*[D + 1];
	for (int i = 0; i < D + 1; i++) {
		coord[i] = new double[D];
		for(int j = 0; j < D; j++) {
			coord[i][j] = rand.get_double(50, 60);
		}
	}

#ifdef VERBOSE	
	cout << "Coordinates of the edges of the " << D << "-Simplex are:\n";
	for (int i = 0; i < D + 1; i++) {
		for(int j = 0; j < D; j++) {
			cout << coord[i][j] << " ";
		}
		cout << '\n';
	}
#endif

	Point *pts = new Point[D+1];
	for (int i = 0; i < D+1; i++) {
		pts[i] = Point(D, coord[i], coord[i] + D);
	}

	vector<double>::iterator it;

	double total_time = 0;
	
	for (int k = 0; k < NUMBER_OF_TESTS; k++) {
#ifdef VERBOSE
		cout << "Coordinates of the randomly generated points are:\n";
#endif

		t.start();
		for(int i = 0; i < total; ++i) {
			output.clear();
			output.reserve(D);
			CGAL::internal::barycoords_d<Point *,
				vector<double>::iterator> (D, pts,output.begin(), rand);
			Point p(D, output.begin(), output.begin() + D);

#ifdef VERBOSE
			for (it = output.begin(); it != output.begin() + D; it++) {
				 cout << *it << " ";
			}
			cout << '\n';
#endif

		}
		t.stop();
		double aux = t.time();
		total_time += aux;
	//	cout << aux << '\n';
		t.reset();
	}

	double avg_time = total_time / NUMBER_OF_TESTS;
	cout << "The average time is " << avg_time << '\n';

	for (int i = 0; i < D + 1; i++) {
		delete coord[i];
	}
	delete coord;
	delete[] pts;
}

int main() {
	benchmark();
	return 0;
}
