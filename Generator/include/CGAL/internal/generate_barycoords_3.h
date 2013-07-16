#ifndef CGAL_INTERNAL_GENERATE_BARYCOORDS_H
#define CGAL_INTERNAL_GENERATE_BARYCOORDS_H

#include <iostream>
#include <vector>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <iterator>
#include <cstdlib>
#include <time.h>

//#define WITH_STD_VECTOR

using namespace std;
using namespace CGAL;

typedef CGAL::Simple_cartesian<double>      K;
typedef K::Point_3                          Point;

namespace CGAL { namespace internal {

	#ifdef WITH_STD_VECTOR // not functional yet
	std::vector<double> generator_3(CGAL::Random& rand) {
		vector<double> a, tmp;
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

		aux = rand.get_double(0,1);

		int ok = 1;
		for (int i = 1; i < 3; i++) {
			tmp.reserve(3);
			if (aux < a[i]) {
				for (int j = 2; j >= i; j--) {
					tmp.push_back(a[j]);
					a.pop_back();
				}
				a.push_back(aux);
				for (int j = 0; j < tmp.size(); j++) {
					a.push_back(tmp[j]);
				}
				tmp.clear();
				ok = 0;
				break;
			}
		}
		if (ok) { // if no change has been made inside the for loop
			a[3] = aux;
		}

		a.push_back(1.0);
		vector<double> b;
		for(int i = 0; i <= 3; ++i) {
			b.push_back(a[i+1]-a[i]);
		}
		return b;
	}

	#else
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

		int ok = 1;
		for (int i = 1; i < 3; i++) {
			if (aux < a[i]) {
				for (int j = i; j < 3; j++) {
					a[j+1] = a[j];
				}
				a[i] = aux;
				ok = 0;
				break;
			}
		}
		if (ok) { // if no change has been made inside the for loop
			a[3] = aux;
		}

		a[4] = 1.0;
		double* b = new double[4];
		for(int i = 0; i <= 3; ++i) {
			b[i] = a[i+1]-a[i];
		}
		return b;
	}
	#endif

	Point add(const Point &p, const Point &q) {
		double coords[3];
		coords[0] = p.x() + q.x();
		coords[1] = p.y() + q.y();
		coords[2] = p.z() + q.z();
		return Point(coords[0], coords[1], coords[2]);
	}

	Point multiply(const Point &p, double c) {
		Point aux(ORIGIN);
		return aux + (c * (p - aux));
	}

	Point assign(const Point &p) {
		return Point(p.x(), p.y(), p.z());
	}

	template <typename RandomAccessIterator, typename OutputIterator>
	void barycoords_3(RandomAccessIterator in, OutputIterator out,
			CGAL::Random& rand) {
		// in contains the coords of the simplex

#ifdef WITH_STD_VECTOR
		vector<double> random;
		random.reserve(3);
#else
		double* random;
#endif
		random = generator_3(rand);

		Point p = add(multiply(in[0], random[0]), multiply(in[1], random[1]));
		p = add(p, multiply(in[2], random[2]));
		p = add(p, multiply(in[3], random[3]));

		for (int i = 0; i < 3; i++) {
			*out = p.cartesian(i);
			out++;
		}
#ifndef WITH_STD_VECTOR
		delete[] random;
#endif
	}
}; // namespace internal
}; // namespace CGAL
#endif // CGAL_INTERNAL_GENERATE_BARYCOORDS_H

