#include <iostream>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <CGAL/double.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cassert>

typedef double RT;

typedef CGAL::Exact_predicates_inexact_constructions_kernel 		K;
typedef K::Point_3 							Point_3;
typedef K::Vector_3 							Vector_3;
typedef K::Triangle_3 							Triangle_3;
typedef K::Plane_3 							Plane_3;
typedef std::vector<Point_3>						Container;
typedef CGAL::Random_points_in_triangle_3<Point_3> 			Point_generator;
	
const double EPS = 1e-30;

template<class InputIterator>
bool inside_or_close_to_triangle(const Triangle_3& tri,InputIterator begin, InputIterator end) {
	Plane_3 plane = tri.supporting_plane();
	while(begin!=end) {
		K::FT dist = squared_distance(plane,*begin);
		if(dist>EPS) { 
			return false;
		}
		Triangle_3 OAB = Triangle_3(*begin,tri[0],tri[1]);
		Triangle_3 OAC = Triangle_3(*begin,tri[0],tri[2]);
		Triangle_3 OBC = Triangle_3(*begin,tri[1],tri[2]);
		K::FT OAB_area = sqrt(OAB.squared_area());
		K::FT OAC_area = sqrt(OAC.squared_area());
		K::FT OBC_area = sqrt(OBC.squared_area());
		K::FT tri_area = sqrt(tri.squared_area());
		if(fabs(OAB_area+OAC_area+OBC_area-tri_area)>1e-15) {
			return false;
		}
		++begin;
	}
	return true;
}

//template<class InputIterator>
//bool is_uniform(const Triangle_3& tri, InputIterator begin, InputIterator end, double r) {
//	Vector_3 vecs[3];
//	for(int i = 0; i < 3; ++i) {
//		vecs[i]=tri[i]-CGAL::ORIGIN;
//	}
//	Vector_3 c(CGAL::NULL_VECTOR);
//	for(int i = 0; i < 3; ++i) {
//		c=c+vecs[i];
//	}	
//	c=c/3.0;
//	for(int i = 0; i < 3; ++i) {
//		vecs[i]=c+r*(vecs[i]-c);
//	}
//	Triangle_3 interior_triangle(CGAL::ORIGIN+vecs[0],CGAL::ORIGIN+vecs[1],
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
	const int MAX_TRIANGLES = 50;
	const int MIN_POINTS = 1000;
	const int MAX_POINTS = 1000000;
	const int number_triangles = rand.get_int(MIN_TRIANGLES,MAX_TRIANGLES);
	const int number_points = rand.get_int(MIN_POINTS, MAX_POINTS);
	for(int i = 0; i < number_triangles; ++i) {
		Point_3 pts[3];
		for(int j = 0; j < 3; ++j) {
			pts[j]=Point_3(rand.get_double(),rand.get_double(),rand.get_double());
		}
		Triangle_3 tri(pts[0],pts[1],pts[2]);
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
