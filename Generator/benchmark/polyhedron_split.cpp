#include <cmath>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

// IO
#include <CGAL/IO/Polyhedron_iostream.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

typedef K::Triangle_3 Triangle;
typedef K::Point_3 Point;

//typedef Polyhedron::Facet Facet;

// Triangulation
//typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
//typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
//typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main()
{
  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input("data/elephant.off");
  input >> polyhedron;
   
  int cont = 0;
  Polyhedron::Facet_iterator it = polyhedron.facets_begin();
  for (; it != polyhedron.facets_end(); ++it) {
	cont++;
	if (it->is_triangle()) {
		std::cout << "Triangle ";
		Point a[3];
		int i = 0;
		Polyhedron::Halfedge_around_facet_circulator iter = it->facet_begin();
		do {
			a[i] = iter->vertex()->point();
			i++;
		} while(++iter != it->facet_begin());
		Triangle aux(a[0], a[1], a[2]);
		std::cout << sqrt(aux.squared_area()) << '\n';
	}
	/*
	if (it->is_quad()) {
		std::cout << "Quadrilateral ";
		Point a[4];
		int i = 0;
		double total_area;
		Polyhedron::Halfedge_around_facet_circulator iter = it->facet_begin();
		do {
			a[i] = iter->vertex()->point();
			i++;
		} while(++iter != it->facet_begin());
		Triangle aux1(a[0], a[1], a[2]);
		Triangle aux2(a[2], a[3], a[0]);
		total_area = sqrt(aux1.squared_area());
		total_area += sqrt(aux2.squared_area());
		std::cout << total_area << '\n';
	}
	*/
	else {
		std::cout << "Not triangle ";
		Point *a = new Point[2];
		int i = 0;
		Polyhedron::Halfedge_around_facet_circulator iter = it->facet_begin();
		do {
			a[i] = iter->vertex()->point();
			i++;
		} while(++iter != it->facet_begin());
		Polyhedron::Halfedge_handle h = polyhedron.create_center_vertex ( it->facet_begin());
		a[i++] = h->vertex()->point();
		double totalArea = 0;
		for (int j = 0; j < i - 1; i++) {
			Triangle aux(a[j], a[j+1], a[i-1]);
		}
		std::cout << totalArea << '\n';
	}
  }

  std::cout << "Number of facets is: " << cont << '\n';
  return 0;
}
