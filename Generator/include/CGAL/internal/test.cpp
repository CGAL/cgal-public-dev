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
#define total 10 // number of points that will be generated

using namespace std;

int main() {
	double coord1[] = {50.5, 60.7};
	double coord2[] = {90.0, 45.4};
	double coord3[] = {42.0, 34.5};
	Point pts[] = {Point(D, coord1,coord1+2, 1),Point(D, coord2,coord2+2,1),Point(D,coord3,coord3+2,1)};

	vector<double> out;
	vector<double>::iterator it;
	for(int i = 0; i < total; ++i) {
		CGAL::internal::barycoords_d<Point *, vector<double>::iterator> (2, pts,out.begin());
		Point p(2, out.begin(), out.end());
		for (it = out.begin(); it != out.end(); it++) {
			cout << *it << " ";
		}
	}

	return 0;
}

