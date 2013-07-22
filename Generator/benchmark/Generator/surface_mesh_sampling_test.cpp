#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Timer.h>
#include <CGAL/internal/ENHANCED_element_sampling.h>
#include <CGAL/point_generators_3.h>
#include <cmath>

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

class PointGen {
	private:
		Triangle_3 t;
	public:
		PointGen() {}

		PointGen(Triangle_3 T) {
			t = T;
		}

//		PointGen(const PointGen &p) {
//			this->t = p.t;
//		}

//		PointGen operator= (const PointGen x) {
//			this->t = x.t;
//			return *this;
//		}

		Point operator() (int q) {
			vector<Point_3> points;
			//CGAL::cpp11::copy_n(CGAL::Random_points_in_tetrahedron_3<Point>(t[0], t[1], t[2], t[3]),
			CGAL::cpp11::copy_n(CGAL::Random_points_in_triangle_3<Point_3>(t),
					1, std::back_inserter(points));
			return points[0];
		}
};

class VolTriangle {
	private:
		Triangle_3 t;
	public:
		VolTriangle(Triangle_3 T) {
			t = T;
		}

		double operator() () {
			return sqrt(t.squared_area());
		}
};

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

	Tr tmp = c2t3.triangulation();
	
	std::cout << "Actual number of facets in tr: " << tr.number_of_finite_facets() << "\n";
	std::cout << "Actual number of facets in c2t3: " << c2t3.number_of_facets() << "\n";
	std::cout << "Actual number of facets in tmp: " <<
		tmp.number_of_finite_facets() << '\n';

	int Nr_cells = c2t3.number_of_facets();
	Triangle_3 *triang;
	triang = new Triangle_3[Nr_cells];
	int i = 0;
	C2t3::Facet_iterator it = c2t3.facets_begin();
	for (; it != c2t3.facets_end(); it++) {
		triang[i] = Triangle_3(c2t3.triangle(*it));
		i++;
	}
//	std::cout << "Number of facets counted by me " << cont << '\n';
//	std::cout << "Areas are: " << '\n';
//	for (int i = 0; i < cont; i++) {
//		std::cout << sqrt(weights[i]) << '\n';
//	}

	CGAL::internal::EnhancedElementSampling <Triangle_3 *, VolTriangle,
		PointGen> (Nr_cells, triangle, triangle+Nr_cells);

	return 0;
}

