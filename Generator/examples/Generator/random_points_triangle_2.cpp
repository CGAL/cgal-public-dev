#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Random.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel 		K;
typedef K::Point_2 							Point_2;
typedef K::Triangle_2 							Triangle_2;
typedef std::vector<Point_2>						Container;
typedef CGAL::Random_points_in_triangle_2<Point_2> 			Point_generator;

void points_in_triangle_2(int N) {
	CGAL::Random rand;
	Container point_set;
	Point_2 pts[3];
	for(int j = 0; j < 3; ++j) {
		pts[j]=Point_2(rand.get_double(),rand.get_double());
	}
	Triangle_2 tri(pts[0],pts[1],pts[2]);
	Point_generator g(tri);
	CGAL::cpp11::copy_n(g, N, std::back_inserter(point_set));
	std::cout << "First of the "<<N<<" generated points is: "<<std::endl;
	std::cout <<point_set[0].x()<<" "<<point_set[0].y()<<std::endl;
}

int main() {
	int N;
	std::cout << "Number of points to be generated inside the triangle:"<<std::endl;
	std::cin >> N;
	points_in_triangle_2(N);
	return 0;
}

