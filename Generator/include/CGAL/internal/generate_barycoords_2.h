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

typedef CGAL::Simple_cartesian<double>      K;
typedef K::Point_2                          Point;

namespace CGAL { namespace internal {

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

	Point add(const Point &p, const Point &q) {
		double coords[2];
		coords[0] = p.x() + q.x();
		coords[1] = p.y() + q.y();
		return Point(coords[0], coords[1]);
	}

	Point multiply(const Point &p, double c) {
		Point aux(ORIGIN);
		return aux + (c * (p - aux));
	}

	Point assign(const Point &p) {
		return Point(p.x(), p.y());
	}

	template <typename RandomAccessIterator, typename OutputIterator>
	void barycoords_2(RandomAccessIterator in, OutputIterator out,
			CGAL::Random& rand) {
		// in contains the coords of the simplex

#ifdef WITH_STD_VECTOR
		vector<double> random;
		random.reserve(3);
#else
		double* random;
#endif
		random = generator_2(rand);

		Point aux = multiply(in[0], random[0]);
		Point tmp = multiply(in[1], random[1]);
		//cout<<aux.x()<<" "<<aux.y()<<" & "<<tmp.x()<<" "<<tmp.y()<<'\n';

		Point p = add(multiply(in[0], random[0]), multiply(in[1], random[1]));
		p = add(p, multiply(in[2], random[2]));

		for (int i = 0; i < 2; i++) {
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

