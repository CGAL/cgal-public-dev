#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/point_generators_3.h>

typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Triangle_3 Triangle_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
typedef typename C2t3::Vertex_handle Vertex;

FT sphere_function2 (Point_3 p) {
  const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  return x2+y2+z2-1;
}

void points_on_surface_mesh_3(int N) {
	Tr tr;
	C2t3 c2t3 (tr);

	Surface_3 surface(sphere_function2,
	      	    Sphere_3(CGAL::ORIGIN, 2.));

	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30., 0.1, 0.1);
	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
	CGAL::Random rand;

	std::vector<Point_3> points;
	points.reserve(N);
	CGAL::Random_points_on_surface_mesh_3<Point_3, C2t3> g(c2t3);
	CGAL::cpp11::copy_n( g, N, std::back_inserter(points));
	std::cout << "First of the "<<N<<" generated points is: "<<std::endl;
	std::cout <<points[0].x()<<" "<<points[0].y()<<" "<<points[0].z()<<std::endl;
	points.clear();
}

int main() {
	int N;
	std::cout << "Number of points to be generated on the surface mesh:"<<std::endl;
	std::cin >> N;
	points_on_surface_mesh_3(N);
	return 0;
}
