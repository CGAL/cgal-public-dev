#include <iostream>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <iterator>
#include <CGAL/enum.h>
#include <CGAL/internal/generate_barycoords_3_IFs.h>
#include <CGAL/Timer.h>

#define D 3 // the dimension is set to 2
#define total 1000000 // number of points that will be generated
#define NUMBER_OF_TESTS 30

//#define VERBOSE

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

	Point *pts = new Point[4];
	for (int i = 0; i < 4; i++) {
		pts[i] = Point(coord[i][0], coord[i][1], coord[i][2]);
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
			CGAL::internal::barycoords_2<Point *,
				vector<double>::iterator> (pts,output.begin(), rand);
			Point p(output[0], output[1], output[2]);

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
		t.reset();
	}

	double avg_time = total_time / NUMBER_OF_TESTS;
	cout << "The average time is " << avg_time << '\n';
}

int main() {
	benchmark();
	return 0;
}
