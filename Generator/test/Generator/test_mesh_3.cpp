#include <iostream>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <CGAL/double.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cassert>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_3/implicit_to_labeled_function_wrapper.h>
#include <CGAL/Mesh_3/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include "implicit_functions.h"
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

using namespace std;
using namespace CGAL::parameters;

typedef double RT;
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel 		K;
typedef FT_to_point_function_wrapper<K::FT, K::Point_3> 		Function;
typedef CGAL::Mesh_3::Implicit_vector_to_labeled_function_wrapper<Function, K>
                                                        		Function_wrapper;
typedef Function_wrapper::Function_vector Function_vector;
typedef CGAL::Mesh_3::Labeled_mesh_domain_3<Function_wrapper, K> 	Mesh_domain;

typedef K::Point_3 							Point;
typedef K::Plane_3 							Plane_3;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type 			Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> 			C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr>			 		Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    		 		Facet_criteria;
typedef Mesh_criteria::Cell_criteria     		 		Cell_criteria;

typedef Tr::Geom_traits 						GT;
typedef GT::Tetrahedron_3 						Tetrahedron3;

typedef CGAL::FasterMemoryExpensiveTag 					FastPolicy;
typedef CGAL::SlowerMemoryEfficientTag 					SlowPolicy;

typedef CGAL::Random_points_in_tetrahedron_3<Point>			 PointGen;
typedef CGAL::internal::Weighted_random_generator<PointGen>		 GeneratorWithWeight;

bool inside_mesh_3(Point pt, const Tetrahedron3 *tet, int size) {
	for (int  i = 0 ; i < size; i++) {
		if(!(!tet[i].bounded_side(pt) == CGAL::ON_UNBOUNDED_SIDE)) {
			return true;
		}
	}
	return false;
}

template<typename Tetrahedron_3>
class WeightFunctor_tetrahedron_3 {
	public:
		double operator() (Tetrahedron_3 &t) {
			return t.volume();
		}
};

int main() {
	// Define functions
	Function f1(&torus_function);
	Function f2(&sphere_function<5>);
	Function f3(&tanglecube_function);
	Function f4(&heart_function);
	Function f5(&klein_function);
	Function f6(&false_knot_function);
	Function f7(&knot1_function);
	Function f8(&octic_function);

	Function_vector v;
	//v.push_back(&f1);
	//v.push_back(&f2);
	//v.push_back(&f3);
	//v.push_back(&f4);
	//v.push_back(&f5);
	v.push_back(&f6);
	//v.push_back(&f7);
	//v.push_back(&f8);

	// Domain (Warning: Sphere_3 constructor uses square radius !)
	Mesh_domain domain(v, K::Sphere_3(CGAL::ORIGIN, 5.*5.), 1e-6);

	// Set mesh criteria
	Facet_criteria facet_criteria(30, 0.2, 0.02); // angle, size, approximation
	Cell_criteria cell_criteria(2., 0.4); // radius-edge ratio, size
	Mesh_criteria criteria(facet_criteria, cell_criteria);

	// Mesh generation
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());

	// Perturbation (maximum cpu time: 10s, targeted dihedral angle: default)
	CGAL::perturb_mesh_3(c3t3, domain, time_limit = 10);
	
	// Exudation
	CGAL::exude_mesh_3(c3t3,12);

	CGAL::Random rand;

	for (int number_points = 1; number_points < 10; number_points++) {
		std::vector<Point> points;
		points.reserve(number_points);
		CGAL::Random_points_in_mesh_3<Point, C3t3> g(c3t3);
		Tetrahedron3 *aux = new
			Tetrahedron3[c3t3.number_of_cells_in_complex()];

		CGAL::cpp11::copy_n(g, number_points, std::back_inserter(points));

		for (int i = 0; i < number_points; i++) {
			typename C3t3::Cells_in_complex_iterator iter =
				c3t3.cells_in_complex_begin();
			int j = 0;
			for (; iter != c3t3.cells_in_complex_end(); ++iter) {
				aux[j] = Tetrahedron3(iter->vertex(0)->point(),
						iter->vertex(1)->point(),
						iter->vertex(2)->point(),
						iter->vertex(3)->point());
				j++;
			}
			assert(inside_mesh_3(points[i], aux, j));
		}
	}

	const int MIN_POINTS = 50;
	const int MAX_POINTS = 100;
	const int number_points = rand.get_int(MIN_POINTS, MAX_POINTS);
	std::vector<Point> points;
	points.reserve(number_points);
	CGAL::Random_points_in_mesh_3<Point, C3t3> g(c3t3);
	Tetrahedron3 *aux = new Tetrahedron3[c3t3.number_of_cells_in_complex()];

	CGAL::cpp11::copy_n(g, number_points, std::back_inserter(points));

	for (int i = 0; i < number_points; i++) {
		typename C3t3::Cells_in_complex_iterator iter =
			c3t3.cells_in_complex_begin();
		int j = 0;
		for (; iter != c3t3.cells_in_complex_end(); ++iter) {
			aux[j] = Tetrahedron3(iter->vertex(0)->point(),
					iter->vertex(1)->point(),
					iter->vertex(2)->point(),
					iter->vertex(3)->point());
			j++;
		}
		assert(inside_mesh_3(points[i], aux, j));
	}

// Testing the copy-constructor
	points.clear();
	points.reserve(number_points);
	CGAL::Random_points_in_mesh_3<Point, C3t3> g1(g);
	delete[] aux;
	aux = new Tetrahedron3[c3t3.number_of_cells_in_complex()];

	CGAL::cpp11::copy_n(g1, number_points, std::back_inserter(points));

	for (int i = 0; i < number_points; i++) {
		typename C3t3::Cells_in_complex_iterator iter =
			c3t3.cells_in_complex_begin();
		int j = 0;
		for (; iter != c3t3.cells_in_complex_end(); ++iter) {
			aux[j] = Tetrahedron3(iter->vertex(0)->point(),
					iter->vertex(1)->point(),
					iter->vertex(2)->point(),
					iter->vertex(3)->point());
			j++;
		}
		assert(inside_mesh_3(points[i], aux, j));
	}

// Testing the constructor that has FSD as argument
	points.clear();
	points.reserve(number_points);

	int Nr_cells_in_cplx = c3t3.number_of_cells_in_complex();
	WeightFunctor_tetrahedron_3<Tetrahedron3> weightElem;
	GeneratorWithWeight* containing_structure;
	containing_structure = new GeneratorWithWeight[Nr_cells_in_cplx];
	typename C3t3::Cells_in_complex_iterator iter =
		c3t3.cells_in_complex_begin();
	int i = 0;
	for (; iter != c3t3.cells_in_complex_end(); ++iter) {
		Tetrahedron3 aux(iter->vertex(0)->point(),
				iter->vertex(1)->point(),
				iter->vertex(2)->point(),
				iter->vertex(3)->point());

		double weight = weightElem(aux);
		PointGen randGen(aux);
		containing_structure[i] = GeneratorWithWeight (randGen, weight);
		i++;
	}
	int N = 1;
//	int N = 1<<10;
	CGAL::internal::Finite_support_distribution<GeneratorWithWeight> fsd =
		CGAL::internal::Finite_support_distribution<GeneratorWithWeight>
		(containing_structure, i, N);

	CGAL::Random_points_in_mesh_3<Point, C3t3> g2(fsd);
	delete[] aux;
	aux = new Tetrahedron3[c3t3.number_of_cells_in_complex()];

	CGAL::cpp11::copy_n(g2, number_points, std::back_inserter(points));

	for (int i = 0; i < number_points; i++) {
		typename C3t3::Cells_in_complex_iterator iter =
			c3t3.cells_in_complex_begin();
		int j = 0;
		for (; iter != c3t3.cells_in_complex_end(); ++iter) {
			aux[j] = Tetrahedron3(iter->vertex(0)->point(),
					iter->vertex(1)->point(),
					iter->vertex(2)->point(),
					iter->vertex(3)->point());
			j++;
		}
		assert(inside_mesh_3(points[i], aux, j));
	}

	points.clear();
	v.clear();
	return 0;
}

