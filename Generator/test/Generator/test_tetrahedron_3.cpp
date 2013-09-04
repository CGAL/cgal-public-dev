#include <iostream>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <CGAL/double.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Cartesian.h>
typedef double RT;

#define VERBOSE

typedef CGAL::Cartesian<RT>					K;
typedef K::Point_3 						Point_3;
typedef K::Vector_3 						Vector_3;
typedef K::Tetrahedron_3 					Tetrahedron_3;
typedef std::vector<Point_3> 					Container;
typedef CGAL::Random_points_in_tetrahedron_3<Point_3> 		Point_generator;

template<class InputIterator>
bool inside_tetrahedron(const Tetrahedron_3& tet,InputIterator begin, InputIterator end) {
    while(begin!=end) {
		if(!tet.bounded_side(*begin) == CGAL::ON_UNBOUNDED_SIDE) {
			std::cout<<"Point outside: "<<*begin<<std::endl;
			std::cout<<"Tetrahedron: "<<tet<<std::endl;
			return false;
		}
		++begin;
    }
    return true;
}

template<class InputIterator>
void is_uniform(const Tetrahedron_3& tet, InputIterator begin, InputIterator end, double r) {
	Vector_3 vecs[4];
	for(int i = 0; i < 4; ++i) {
		vecs[i]=tet[i]-CGAL::ORIGIN;
	}
	Vector_3 c(CGAL::NULL_VECTOR);
	for(int i = 0; i < 4; ++i) {
		c=c+vecs[i];
	}	
	c=c/4.0;
	Point_3 pts[4];
	for(int i = 0; i < 4; ++i) {
		vecs[i]=c+r*(vecs[i]-c);
		pts[i]=CGAL::ORIGIN+vecs[i];
	}
	Tetrahedron_3 interior_tetrahedron(pts[0],pts[1],pts[2],pts[3]);
	int inside_smaller = 0;
	int total=0;
	while(begin!=end) {
		switch(interior_tetrahedron.bounded_side(*begin)) {
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
	std::cout<<"number of points inside smaller tetrahedron: "<<inside_smaller<<std::endl;
	std::cout<<"expected number: "<<r*r*r*total<<std::endl;
	#endif
}


int main() {
    CGAL::Random rand;
    Container point_set;
    int number_tetrahedrons;
#ifdef VERBOSE
    std::cout<<"Type the number of tetrahedrons: ";
#endif
    std::cin>>number_tetrahedrons;
    int number_points;
#ifdef VERBOSE
    std::cout<<"Type the number of points inside each tetrahedron: ";
#endif 
    std::cin>>number_points;
#ifdef VERBOSE
    std::cout<<"Type the factor: ";//the parameter r in function is_uniform
#endif 
    double r;
    std::cin>>r;
    for(int i = 0; i < number_tetrahedrons; ++i) {
		Point_3 pts[4];
		for(int j = 0; j < 4; ++j) {
			pts[j]=Point_3(rand.get_double(),rand.get_double(),rand.get_double());
		}
		Tetrahedron_3 tet(pts[0],pts[1],pts[2],pts[3]);
		Point_generator g1( pts[0], pts[1], pts[2], pts[3] ); // constructor that is given points
		Point_generator g2( tet ); // constructor that is given a tetrahedron
		Point_generator g3( g1 ); // copy-constructor

		//Testing the point-constructor
		point_set.clear();
		CGAL::cpp11::copy_n( g1, number_points,
		               std::back_inserter(point_set));
		if(!inside_tetrahedron(tet,point_set.begin(),point_set.end())) {
		    std::cout<<"POINT OUTSIDE TETRAHEDRON\n"<<std::endl;
		}
		is_uniform(tet,point_set.begin(),point_set.end(),r);

		//Testing the tetrahedron-constructor
		point_set.clear();
		CGAL::cpp11::copy_n( g2, number_points,
		               std::back_inserter(point_set));
		if(!inside_tetrahedron(tet,point_set.begin(),point_set.end())) {
		    std::cout<<"POINT OUTSIDE TETRAHEDRON\n"<<std::endl;
		}
		is_uniform(tet,point_set.begin(),point_set.end(),r);

		//Testing the copy-constructor
		point_set.clear();
		CGAL::cpp11::copy_n( g3, number_points,
		               std::back_inserter(point_set));
		if(!inside_tetrahedron(tet,point_set.begin(),point_set.end())) {
		    std::cout<<"POINT OUTSIDE TETRAHEDRON\n"<<std::endl;
		}
		is_uniform(tet,point_set.begin(),point_set.end(),r);
    }
    return 0;
}
