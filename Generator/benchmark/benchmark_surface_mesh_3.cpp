#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Timer.h>
#include <CGAL/point_generators_3.h>
#include <cmath>
#include <CGAL/internal/Finite_support_distribution.h>
#include <CGAL/internal/Weighted_random_generator.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3       Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr>           C2t3;

typedef Tr::Geom_traits                                  GT;
typedef GT::Sphere_3                                     Sphere_3;
typedef GT::Point_3                                      Point_3;
typedef GT::Triangle_3                                   Triangle_3;
typedef GT::FT                                           FT;

typedef FT                                               (*Function)(Point_3);

typedef CGAL::Implicit_surface_3<GT, Function>           Surface_3;

FT sphere_function (Point_3 p) {
	const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
	return x2+y2+z2-1;
}

//Random generator
typedef CGAL::Random_points_in_triangle_3<Point_3> PointGen;
typedef CGAL::internal::Weighted_random_generator<PointGen>
	GeneratorWithWeight;

int main() {
	Tr tr;            // 3D-Delaunay triangulation
	C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
	
	// defining the surface
	Surface_3 surface(sphere_function,             // pointer to function
	                  Sphere_3(CGAL::ORIGIN, 2.)); // bounding sphere
	// Note that "2." above is the *squared* radius of the bounding sphere!
	
	// defining meshing criteria
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
	                                                   0.1,  // radius bound
	                                                   0.1); // distance bound
	// meshing surface
	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

	int nr = 10000;
	std::vector<Point_3> points;
	points.reserve(nr);
	CGAL::Random_points_in_surface_mesh_3<Point_3, C2t3> g(c2t3);

	CGAL::cpp11::copy_n(g, nr, std::back_inserter(points));
	std::cout << "The generated points are: " << std::endl;
	for (int i = 0; i < nr; i++) {
		std::cout << points[i].x() << "\n" << points[i].y() << "\n" <<
			points[i].z() << std::endl;
	}
	return 0;
}

