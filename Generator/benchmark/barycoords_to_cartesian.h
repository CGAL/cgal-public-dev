#ifndef CGAL_INTERNAL_BARYCOORDS_TO_CARTESIAN_H
#define CGAL_INTERNAL_BARYCOORDS_TO_CARTESIAN_H
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
typedef Cartesian_d<double>                 Kd;
typedef CGAL::Simple_cartesian<double>      K2;
typedef CGAL::Simple_cartesian<double>      K3;
typedef Point_d< Kd >                       Pointd;
typedef K2::Point_2                         Point2;
typedef K3::Point_3                         Point3;


namespace CGAL { namespace internal {

	// dD case
	Pointd add(const Pointd &p, const Pointd &q) {
		int d = p.dimension();
		vector<double> coords;
		coords.reserve(d);
		for (int i = 0; i < d; i++) {
			coords.push_back(p.cartesian(i) + q.cartesian(i));
		}
		return Pointd(d, coords.begin(), coords.end());
	}

	Pointd multiply(const Pointd &p, double c) {
		int d = p.dimension();
		Pointd aux(d,ORIGIN);
		return aux + (c * (p - aux));
	}

	Pointd assign(const Pointd &p) {
		int d = p.dimension();
		return Pointd(d, p.cartesian_begin(), p.cartesian_end());
	}

	template <typename RandomAccessIterator, typename OutputIterator>
	void barycoords_d(int d, RandomAccessIterator in, OutputIterator out,
			CGAL::Random& rand) {
		// in contains the coords of the simplex

		vector<double> random;
		random.clear();
		random = generator(d, rand);
		Pointd p(d, ORIGIN);

		for (int i = 0; i <= d; i++) {
			p = add(p, multiply(in[i], random[i]));
		}

		for (int i = 0; i < d; i++) {
			*out = p.cartesian(i);
			out++;
		}
	}

	// 2D optimized case
	Point2 add(const Point2 &p, const Point2 &q) {
		double coords[2];
		coords[0] = p.x() + q.x();
		coords[1] = p.y() + q.y();
		return Point2(coords[0], coords[1]);
	}

	Point2 multiply(const Point2 &p, double c) {
		Point2 aux(ORIGIN);
		return aux + (c * (p - aux));
	}

	Point2 assign(const Point2 &p) {
		return Point2(p.x(), p.y());
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

		Point2 aux = multiply(in[0], random[0]);
		Point2 tmp = multiply(in[1], random[1]);
		//cout<<aux.x()<<" "<<aux.y()<<" & "<<tmp.x()<<" "<<tmp.y()<<'\n';

		Point2 p = add(multiply(in[0], random[0]), multiply(in[1], random[1]));
		p = add(p, multiply(in[2], random[2]));

		for (int i = 0; i < 2; i++) {
			*out = p.cartesian(i);
			out++;
		}
#ifndef WITH_STD_VECTOR
		delete[] random;
#endif
	}

	// 3D optimized case
	Point3 add(const Point3 &p, const Point3 &q) {
		double coords[3];
		coords[0] = p.x() + q.x();
		coords[1] = p.y() + q.y();
		coords[2] = p.z() + q.z();
		return Point3(coords[0], coords[1], coords[2]);
	}

	Point3 multiply(const Point3 &p, double c) {
		Point3 aux(ORIGIN);
		return aux + (c * (p - aux));
	}

	Point3 assign(const Point3 &p) {
		return Point3(p.x(), p.y(), p.z());
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

		Point3 p = add(multiply(in[0], random[0]), multiply(in[1], random[1]));
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
#endif // CGAL_INTERNAL_BARYCOORDS_TO_CARTESIAN
