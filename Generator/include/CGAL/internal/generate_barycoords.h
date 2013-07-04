#ifndef CGAL_INTERNAL_GENERATE_BARYCOORDS_H
#define CGAL_INTERNAL_GENERATE_BARYCOORDS_H

#include <iostream>
#include <vector>
#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <iterator>
#include <cstdlib>
#include <time.h>

//#define WITH_STD_VECTOR

using namespace std;
using namespace CGAL;
typedef Cartesian_d<double>          K;
typedef Point_d< K >                 Point;
typedef Vector_d< K >                Vector;

namespace CGAL { namespace internal {

	// general case for dD space
	std::vector<double> generator(int d, CGAL::Random& rand) {
		vector<double> a;
		a.push_back(0.0);
		for(int i = 0; i < d; ++i) {
			a.push_back(rand.get_double(0,1));
		}
		a.push_back(1.0);
		sort(a.begin(),a.end());
		vector<double> b;
		for(int i = 0; i <= d; ++i) {
			b.push_back(a[i+1]-a[i]);
		}
		return b;
	}

	// optimized implementation for 2D space
#ifdef WITH_STD_VECTOR
	std::vector<double> generator_2(CGAL::Random& rand) {
		vector<double> a;
		a.push_back(0.0);
		a.push_back(rand.get_double(0,1));
		double aux = rand.get_double(0,1);

		if (aux < a[1]) {
			double tmp = a[1];
			a.pop_back();
			a.push_back(aux);
			a.push_back(tmp);
		} else {
			a.push_back(aux);
		}

		a.push_back(1.0);
		vector<double> b;
		for(int i = 0; i <= 2; ++i) {
			b.push_back(a[i+1]-a[i]);
		}
		return b;
	}

#else
	double* generator_2(CGAL::Random& rand) {
		double a[4];
		a[0] = 0.0;
		a[1] = rand.get_double(0,1);
		double aux = rand.get_double(0,1);

		if (aux < a[1]) {
			double tmp = a[1];
			a[1] = aux;
			a[2] = tmp;
		} else {
			a[2] = aux;
		}

		a[3] = 1.0;
		double* b = new double[3];
		for(int i = 0; i <= 2; i++) {
			b[i] = a[i+1]-a[i];
		}
		return b;
	}
#endif

	// optimized implementation for 3D space
	double* generator_3(CGAL::Random& rand) {
		double a[5];
		a[0] = 0.0;
		a[1] = rand.get_double(0,1);
		double aux = rand.get_double(0,1);

		if (aux < a[1]) {
			double tmp = a[1];
			a[1] = aux;
			a[2] = tmp;
		} else {
			a[2] = aux;
		}

		aux = rand.get_double(0,1);

		if (aux < a[1]) {
			double tmp = a[1];
			a[1] = aux;
			a[2] = tmp;
		} else {
			if (aux > a[2]) {
				a[3] = aux;
			} else {
				double tmp = a[2];
				a[2] = aux;
				a[3] = tmp;
			}
		}

		a[4] = 1.0;
		double* b = new double[4];
		for(int i = 0; i <= 3; ++i) {
			b[i] = a[i+1]-a[i];
		}
		return b;
	}
}; // namespace internal
}; // namespace CGAL
#endif // CGAL_INTERNAL_GENERATE_BARYCOORDS_H

