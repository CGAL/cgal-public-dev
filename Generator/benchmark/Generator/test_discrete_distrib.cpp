#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Timer.h>
#include <CGAL/internal/Discrete_distribution_with_finite_support_generator.h>
#include <CGAL/internal/Random_generator_with_weight.h>
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
{return CGAL::squared_distance(p, Point(CGAL::ORIGIN))-1; }

//Random generator
typedef CGAL::Random_points_in_tetrahedron_3<Point> PointGen;
typedef CGAL::internal::Random_generator_with_weight<PointGen>
	GeneratorWithWeight;

class WeightFunctor {
	public:
		double operator() (Tetrahedron3 t) {
			return t.volume();
		}
};

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

	WeightFunctor weightElem;
	Tr tr = c3t3.triangulation();
	std::vector<GeneratorWithWeight> containing_structure;
	containing_structure.reserve(Nr_cells_in_cplx);
	Tr::Finite_cells_iterator iter = tr.finite_cells_begin();
	for (; iter != tr.finite_cells_end(); ++iter) {
		if (c3t3.is_in_complex(iter)) {
			Tetrahedron3 aux = tr.tetrahedron(iter);
			double weight = weightElem(aux);
			PointGen randGen(aux);
			GeneratorWithWeight tmp = GeneratorWithWeight (randGen, weight);
			containing_structure.push_back(tmp);
		}
	}

	CGAL::internal::Discrete_distribution_with_finite_support_generator<GeneratorWithWeight
		> randomGen(containing_structure);

	randomGen.generate();
	return 0;
}
