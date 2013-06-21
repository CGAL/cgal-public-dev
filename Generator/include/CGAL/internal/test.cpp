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

#define D 2 // the d dimension
#define total 1000000 // number of points that will be generated

using namespace std;
vector<double> output;

int main() {
	double coord1[] = {50.5, 60.7};
	double coord2[] = {90.0, 45.4};
	double coord3[] = {42.0, 34.5};
	Point pts[] = {Point(D, coord1,coord1+2, 1),Point(D, coord2,coord2+2,1),Point(D,coord3,coord3+2,1)};

	vector<double>::iterator it;
	for(int i = 0; i < total; ++i) {
		output.clear();
		output.reserve(D);
		CGAL::internal::barycoords_d<Point *, vector<double>::iterator> (2, pts,output.begin());
		Point p(2, output.begin(), output.begin() + D);
		cout << "Coordinates for this point are: ";
		for (it = output.begin(); it != output.begin() + D; it++) {
			 cout << *it << " ";
		}
		cout << '\n';
	}
	return 0;
}

