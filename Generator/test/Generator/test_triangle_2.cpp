#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>

#include <CGAL/double.h>

typedef double RT;

#define VERBOSE

typedef CGAL::Exact_predicates_inexact_constructions_kernel		K;
typedef K::Point_2							Point_2;
typedef K::Vector_2							Vector_2;
typedef K::Triangle_2							Triangle_2;
typedef std::vector<Point_2>						Container;
typedef CGAL::Random_points_in_triangle_2<Point_2>			Point_generator;

const double EPS = 1e-30;

template<class InputIterator>
bool inside_or_close_to_triangle(const Triangle_2& tri,InputIterator begin, InputIterator end) {
	while(begin!=end) {
		K::FT dist = squared_distance(tri,*begin);
		if(dist>EPS) {
			return false;
		}
		++begin;
	}
	return true;
}

//template<class InputIterator>
//bool is_uniform(const Triangle_2& tri, InputIterator begin, InputIterator end, double r) {
//	Vector_2 vecs[3];
//	for(int i = 0; i < 3; ++i) {
//		vecs[i]=tri[i]-CGAL::ORIGIN;
//	}
//	Vector_2 c(CGAL::NULL_VECTOR);
//	for(int i = 0; i < 3; ++i) {
//		c=c+vecs[i];
//	}	
//	c=c/3.0;
//	for(int i = 0; i < 3; ++i) {
//		vecs[i]=c+r*(vecs[i]-c);
//	}
//	Triangle_2 interior_triangle(CGAL::ORIGIN+vecs[0],CGAL::ORIGIN+vecs[1],
//									CGAL::ORIGIN+vecs[2]);
//	int inside_smaller = 0;
//	int total=0;
//	while(begin!=end) {
//		switch(interior_triangle.bounded_side(*begin)) {
//			case CGAL::ON_BOUNDED_SIDE:
//				++inside_smaller;
//			break;
//			
//			case CGAL::ON_BOUNDARY:
//				++inside_smaller;
//			break;
//		}
//		++begin,++total;
//	}
//	return fabs(inside_smaller-r*r*total)<total/100.0;
//}

int main() {
	CGAL::Random rand;
	Container point_set;
	const int MIN_TRIANGLES = 1;
	const int MAX_TRIANGLES = 10;
	const int MIN_POINTS = 1000;
	const int MAX_POINTS = 1000000;
	const int number_triangles = rand.get_int(MIN_TRIANGLES,MAX_TRIANGLES);
	const int number_points = rand.get_int(MIN_POINTS, MAX_POINTS);
	const double r = rand.get_double(0.2,0.8);
	for(int i = 0; i < number_triangles; ++i) {
		Point_2 pts[3];
		for(int j = 0; j < 3; ++j) {
			pts[j]=Point_2(rand.get_double(),rand.get_double());
		}
		Triangle_2 tri(pts[0],pts[1],pts[2]);
		Point_generator g1( pts[0], pts[1], pts[2] );
		Point_generator g2( tri );
		Point_generator g3( g1 );

		point_set.clear();
		CGAL::cpp11::copy_n( g1, number_points,
		               std::back_inserter(point_set));
		assert(inside_or_close_to_triangle(tri,point_set.begin(),point_set.end()));

		point_set.clear();
		CGAL::cpp11::copy_n( g2, number_points,
		               std::back_inserter(point_set));
		assert(inside_or_close_to_triangle(tri,point_set.begin(),point_set.end()));

		point_set.clear();
		CGAL::cpp11::copy_n( g3, number_points,
		               std::back_inserter(point_set));
		assert(inside_or_close_to_triangle(tri,point_set.begin(),point_set.end()));
	}
   return 0;
}
