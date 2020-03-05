#pragma once
#pragma warning(disable:4005)

#ifndef Q_MOC_RUN
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#endif

namespace Skippy {
	typedef unsigned char uchar;
	typedef unsigned int uint;

	typedef enum {
		PLUS = 1,
		ZERO = 0,
		MINUS = -1
	} Sign;

	typedef enum {
		NONE,
		CONVEX_ENVELOPS,
		OPTIMAL_RECTANGLES
	} Preprocess;

	template <typename T>
	const T & jmin(const T & a, const T & b) {
		return (a < b ? a : b);
	}

	template <typename T>
	const T & jmax(const T & a, const T & b) {
		return (a > b ? a : b);
	}

	template <typename T>
	const T & jpclamp(const T & a, const T & x, const T & b) {
		return (x < a ? a : (x > b ? b : x));
	}

	template <typename T>
	bool jin(const T & a, const T & x, const T & b) {
		return (a <= x && x <= b);
	}

	template <typename T>
	bool jinr(const T & a, const T & x, const T & b) {
		return (a < x && x < b);
	}

	namespace bg = boost::geometry;
	namespace bgi = boost::geometry::index;

	typedef bg::model::point<double, 2, bg::cs::cartesian> Boost_Point;
	typedef bg::model::box<Boost_Point> Boost_Box;
	typedef std::pair<Boost_Box, unsigned int> Boost_Value;
	typedef bgi::rtree<Boost_Value, bgi::quadratic<16> > Boost_RTree;
}


#ifndef PI
#define PI 3.141592653589783238462643383279
#endif


/*#ifndef COS_PI_360
#define COS_PI_360 0.99996192306
#endif*/


#ifdef KINETIC_PARTITION_EXPORTS
#define KINETIC_PARTITION_API __declspec(dllexport)
#else
#define KINETIC_PARTITION_API __declspec(dllimport)
#endif