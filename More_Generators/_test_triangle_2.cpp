#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>

#include <CGAL/double.h>
#include "Random_points_in_triangle_2.h"

typedef double RT;

#define VERBOSE

typedef CGAL::Simple_cartesian<RT>                        K;
typedef K::Point_2                                        Point_2;
typedef K::Vector_2										  Vector_2;
typedef K::Triangle_2									  Triangle_2;
typedef std::vector<Point_2>                              Container;
typedef CGAL::Random_points_in_triangle_2<Point_2>        Point_generator;

template<class InputIterator>
bool inside_triangle(const Triangle_2& tri,InputIterator begin, InputIterator end) {
	while(begin!=end) {
		if(tri.bounded_side(*begin) == CGAL::ON_UNBOUNDED_SIDE) {
				std::cout<<"Point outside: "<<*begin<<std::endl;
				std::cout<<"Triangle: "<<tri<<std::endl;
				return false;
		}
		++begin;
	}
	return true;
}


//creates a smaller triangle inside the larger with
//area tri.area()*r*r. checks if fraction of points
//inside smaller triangle equals r*r 

//has to be fixed
template<class InputIterator>
bool is_uniform(const Triangle_2& tri, InputIterator begin, InputIterator end, double r) {
	Vector_2 vecs[3];
	for(int i = 0; i < 3; ++i) {
		vecs[i]=tri[i]-CGAL::ORIGIN;
	}
	Vector_2 c(CGAL::NULL_VECTOR);
	for(int i = 0; i < 3; ++i) {
		c=c+vecs[i];
	}	
	c=c/3.0;
	for(int i = 0; i < 3; ++i) {
		vecs[i]=c+r*(vecs[i]-c);
	}
	Triangle_2 interior_triangle(CGAL::ORIGIN+vecs[0],CGAL::ORIGIN+vecs[1],
									CGAL::ORIGIN+vecs[2]);
	int inside_smaller = 0;
	int total=0;
	while(begin!=end) {
		switch(interior_triangle.bounded_side(*begin)) {
			case CGAL::ON_BOUNDED_SIDE:
				++inside_smaller;
			break;
			
			case CGAL::ON_BOUNDARY:
				++inside_smaller;
			break;
		}
		++begin,++total;
	}
	#ifdef VERBOSE
	std::cout<<"number of points inside smaller triangle: "<<inside_smaller<<std::endl;
	std::cout<<"expected number: "<<r*r*total<<std::endl;
	#endif
	//this has to be changed
	return fabs(inside_smaller-r*r*total)<total/100.0;
}

int main() {
   CGAL::Random rand;
   Container point_set;
   int number_triangles;
   #ifdef VERBOSE
   std::cout<<"Type the number of triangles: ";
   #endif
   std::cin>>number_triangles;
   int number_points;
   #ifdef VERBOSE
   std::cout<<"Type the number of points inside each triangle: ";
   #endif 
   std::cin>>number_points;
   #ifdef VERBOSE
   std::cout<<"Type the factor: ";//the parameter r in function is_uniform
   #endif 
   double r;
   std::cin>>r;
   for(int i = 0; i < number_triangles; ++i) {
   		Point_2 pts[3];
   		for(int j = 0; j < 3; ++j) {
   			pts[j]=Point_2(rand.get_double(),rand.get_double());
   		}
   		Triangle_2 tri(pts[0],pts[1],pts[2]);
   		point_set.clear();
   		CGAL::copy_n_unique(Point_generator( pts[0], pts[1], pts[2] ), 
                       number_points,
                       std::back_inserter(point_set));
   		if(!inside_triangle(tri,point_set.begin(),point_set.end())) {
   			std::cout<<"POINT OUTSIDE TRIANGLE\n";
   		}
   		if(!is_uniform(tri,point_set.begin(),point_set.end(),r)) {
   			std::cout<<"GENERATOR IS NOT UNIFORM\n";
   		}
   }
   return 0;
}
