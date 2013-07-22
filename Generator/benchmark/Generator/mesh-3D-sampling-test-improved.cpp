#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Timer.h>
#include <CGAL/internal/ENHANCED_element_sampling.h>
#include <CGAL/point_generators_3.h>

using namespace std;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT (Function)(const Point&);
typedef CGAL::Implicit_mesh_domain_3<Function,K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

//typedef CGAL::MeshComplex_3InTriangulation_3::Triangulation Trig;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

typedef Tr::Geom_traits GT;
typedef GT::Tetrahedron_3 Tetrahedron3;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Function
FT sphere_function (const Point& p)
{ return CGAL::squared_distance(p, Point(CGAL::ORIGIN))-1; }

class PointGen {
	private:
		Tetrahedron3 t;
	public:
		PointGen() {}

		PointGen(Tetrahedron3 T) {
			t = T;
		}

		PointGen(const PointGen &p) {
			this->t = p.t;
		}

		PointGen operator= (const PointGen x) {
			this->t = x.t;
			return *this;
		}

		Point operator() (int q) {
			vector<Point> points;
			CGAL::cpp11::copy_n(CGAL::Random_points_in_tetrahedron_3<Point>(t[0], t[1], t[2], t[3]),
					1, std::back_inserter(points));
			return points[0];
		}
};

class VolTetrahedron {
	private:
		Tetrahedron3 t;
	public:
		VolTetrahedron(Tetrahedron3 T) {
			t = T;
		}

		double operator() () {
			return t.volume();
		}
};

//double VolTetrahedron(Tetrahedron3 &t) {
//	return t.volume();
//}

int main()
{
	CGAL::Timer timp;
	// Domain (Warning: Sphere_3 constructor uses squared radius !)
	Mesh_domain domain(sphere_function, K::Sphere_3(CGAL::ORIGIN, 2.));
	
	// Mesh criteria
	Mesh_criteria criteria(facet_angle=30, facet_size=0.1, facet_distance=0.025,
	                       cell_radius_edge_ratio=2, cell_size=0.1);
	
	// Mesh generation
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
	
	// Output
	std::ofstream medit_file("out.mesh");
	c3t3.output_to_medit(medit_file);
	
	cout << "Actual number of cells in c3t3 in complex: " << c3t3.number_of_cells_in_complex() << "\n";
	
	
	int Nr_cells_in_cplx = c3t3.number_of_cells_in_complex();
	int i = 0;
	C3t3::Cells_in_complex_iterator iter = c3t3.cells_in_complex_begin();
	for ( ; iter != c3t3.cells_in_complex_end(); ++iter) {
		//tetra[i] = tr.tetrahedron(iter);
		i++;
	}


	cout << "Cells in complex counted by me: " << i << '\n';
	Tr tr = c3t3.triangulation();
	int Nr_cells = tr.number_of_finite_cells();
	cout << "Actual number of cells in triangulation: " <<
		Nr_cells << "\n";

	Tetrahedron3 *tetra;
	tetra = new Tetrahedron3[Nr_cells_in_cplx];
	i = 0;
	Tr::Finite_cells_iterator iter2 = tr.finite_cells_begin();
	for ( ; iter2 != tr.finite_cells_end(); ++iter2) {
		if (c3t3.is_in_complex(iter2)) {
			tetra[i] = tr.tetrahedron(iter2);
			i++;
		} else {
			Nr_cells--;
		}
	}
	cout << "Actual number of cells in triangulation afterwards: " <<
		Nr_cells << "\n";
	
	CGAL::internal::EnhancedElementSampling <Tetrahedron3 *, VolTetrahedron,
		PointGen> (Nr_cells, tetra, tetra+Nr_cells);

//	timp.start();
//	for (i = 0; i < Nr_cells; i++) {
//		cout << tetra[i].volume() << '\n';
//	}
//	timp.stop();
//	cout << "Now the total time is: " << timp.time() << '\n';

	delete[] tetra;
	return 0;
}
