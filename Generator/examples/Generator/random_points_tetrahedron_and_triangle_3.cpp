#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel 		K;
typedef K::Point_3 							Point_3;
typedef K::Tetrahedron_3						Tetrahedron_3;
typedef K::Triangle_3 							Triangle_3;
typedef std::vector<Point_3>						Container;
typedef CGAL::Random_points_in_tetrahedron_3<Point_3>			Point_generator_tetra;
typedef CGAL::Random_points_in_triangle_3<Point_3>			Point_generator_tri;

void points_in_triangle_3(int N) {
	CGAL::Random rand;
	Container point_set;
	Point_3 pts[3];
	for(int j = 0; j < 3; ++j) {
		pts[j]=Point_3(rand.get_double(),rand.get_double(),rand.get_double());
	}
	Triangle_3 tri(pts[0],pts[1],pts[2]);
	Point_generator_tri g(tri);
	CGAL::cpp11::copy_n(g, N, std::back_inserter(point_set));
	std::cout << "First of the "<<N<<" points generated inside the triangle is:"<<std::endl;
	std::cout <<point_set[0].x()<<" "<<point_set[0].y()<<" "<<point_set[0].z()<<std::endl;
}

void points_in_tetrahedron_3(int N) {
	CGAL::Random rand;
	Container point_set;
	Point_3 pts[3];
	for(int j = 0; j < 4; ++j) {
		pts[j]=Point_3(rand.get_double(),rand.get_double(),rand.get_double());
	}
	Tetrahedron_3 tetra(pts[0],pts[1],pts[2],pts[3]);
	Point_generator_tetra g(tetra);
	CGAL::cpp11::copy_n(g, N, std::back_inserter(point_set));
	std::cout << "First of the "<<N<<" points generated inside the tetrahedron is:"<<std::endl;
	std::cout <<point_set[0].x()<<" "<<point_set[0].y()<<" "<<point_set[0].z()<<std::endl;
}

int main() {
	int N;
	std::cout << "Number of points to be generated inside the triangle:"<<std::endl;
	std::cin >> N;
	points_in_triangle_3(N);
	std::cout << "Number of points to be generated inside the tetrahedron:"<<std::endl;
	std::cin >> N;
	points_in_tetrahedron_3(N);
	return 0;
}
