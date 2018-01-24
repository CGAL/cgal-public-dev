#ifndef JEAN_PHILIPPE_DEFS_H
#define JEAN_PHILIPPE_DEFS_H

/*#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>*/
#ifndef Q_MOC_RUN
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#endif

#ifndef minjp
#define minjp(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef maxjp
#define maxjp(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef jclamp
#define jclamp(a, x, b) ((x) < (a) ? (a) : ((x) > (b) ? (b) : (x)))
#endif

#ifndef PI
#define PI 3.141592653589783238462643383279
#endif
typedef unsigned char uchar;
typedef unsigned int uint;

/*typedef double FT;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 CGAL_Point2d;
typedef K::Point_3 CGAL_Point3d;
typedef K::Vector_2 CGAL_Vec2d;
typedef K::Vector_3 CGAL_Vec3d;
typedef K::Plane_3 CGAL_Plane;
typedef K::Segment_2 CGAL_Segment_2;
typedef K::Triangle_2 CGAL_Triangle2d;
typedef K::Triangle_3 CGAL_Triangle3d;

typedef CGAL::Color CGAL_Color;
typedef CGAL::Simple_cartesian<int>::Point_2 CGAL_Point2i;
typedef CGAL::Simple_cartesian<int>::Point_3 CGAL_Point3i;*/


namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> Boost_Point;
typedef bg::model::box<Boost_Point> Boost_Box;
typedef std::pair<Boost_Box, uint> Boost_Value;
typedef bgi::rtree<Boost_Value, bgi::quadratic<16> > Boost_RTree;

#define NOT_MEASURING_PERFORMANCES 1

#endif