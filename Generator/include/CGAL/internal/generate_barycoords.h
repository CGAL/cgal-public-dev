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

using namespace std;
using namespace CGAL;
typedef Cartesian_d<double>   K;
typedef Point_d< K >                 Point;
typedef Vector_d< K >                 Vector;

namespace CGAL { namespace internal {

	std::vector<double> generator(int d) {
		CGAL::Random rand;
		//srand(time(NULL));
		vector<double> a;
		a.push_back(0.0);
		for(int i = 0; i < d; ++i) {
			a.push_back(rand.get_double(0,1));
			//a.push_back((double) rand() / RAND_MAX);
		}
		a.push_back(1.0);
		sort(a.begin(),a.end());
		vector<double> b;
		for(int i = 0; i <= d; ++i) {
			b.push_back(a[i+1]-a[i]);
		}
		return b;
	}

	Point add(const Point &p, const Point &q) {
		int d = p.dimension();
		vector<double> coords;
		coords.reserve(d);
		for (int i = 0; i < d; i++) {
			coords.push_back(p.cartesian(i) + q.cartesian(i));
		}
		return Point(d, coords.begin(), coords.end());
	}

	Point multiply(const Point &p, double c) {
		int d = p.dimension();
		Point aux(d,ORIGIN);
		return aux + (c * (p - aux));
	}

	Point assign(const Point &p) {
		int d = p.dimension();
		return Point(d, p.cartesian_begin(), p.cartesian_end());
	}

	template <typename RandomAccessIterator, typename OutputIterator>
	void barycoords_d(int d, RandomAccessIterator in, OutputIterator out) {
		// in contains the coords of the simplex
		vector<double> random;
		random.clear();
		random = generator(d);
		Point p(d, ORIGIN);
//		cout << "What rand returns: " << random[0] << " " <<random[1] << '\n';
		for (int i = 0; i < d; i++) {
			p = add(p, multiply(in[i], random[i]));
			in++;
		}
		for (int i = 0; i < d; i++) {
			*out = p.cartesian(i);
			out++;
		}
	}
}; // namespace internal
}; // namespace CGAL
#endif // CGAL_INTERNAL_GENERATE_BARYCOORDS_H

