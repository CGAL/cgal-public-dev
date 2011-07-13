#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/collapse.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/IO/Color.h>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#include <CGAL/Random.h>

#include <iostream>
#include <cstring>
#include <fstream>
#include <list>
#include <set>

using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel		K;
typedef CGAL::Triangulation_vertex_base_with_info_3<CGAL::Color, K>	Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>			Tds;
typedef CGAL::Triangulation_3<K, Tds>					Triangulation2;
typedef CGAL::Delaunay_triangulation_3<K, Tds>				Triangulation;

typedef K::FT FT;

typedef Triangulation::Finite_vertices_iterator 	Finite_vertices_iterator;
typedef Triangulation::Finite_edges_iterator 		Finite_edges_iterator;
typedef Triangulation::Finite_facets_iterator 		Finite_facets_iterator;
typedef Triangulation::Finite_cells_iterator 		Finite_cells_iterator;
typedef Triangulation::Simplex        			Simplex;
typedef Triangulation::Locate_type    			Locate_type;
typedef Triangulation::Point          			Point;

typedef Triangulation::Vertex_handle			Vertex_handle;
typedef Triangulation::Cell_handle			Cell_handle;

CGAL::Random my_rand(42);

Triangulation my_triangulation()
{
	std::list<Point> L;

	L.push_front(Point(0,0,0));
	L.push_front(Point(.5,0,0));
	L.push_front(Point(0,.5,0));
	L.push_front(Point(0,0,.5));
	L.push_front(Point(.5,.5,.5));
	L.push_front(Point(.2,.3,.1));
	L.push_front(Point(.3,.1,.4));
	L.push_front(Point(.9,.5,.8));
	L.push_front(Point(.7,.7,.1));
	L.push_front(Point(.1,.6,.9));
	L.push_front(Point(.3,.1,.4));

	return Triangulation(L.begin(), L.end());
}

Triangulation get_rand_triangulation(int n)
{
	int i;
	std::list<Point> L;

	for (i=0; i<n; i++) {
		double x = my_rand.get_double(0.0, 1.0);
		double y = my_rand.get_double(0.0, 1.0);
		double z = my_rand.get_double(0.0, 1.0);

		L.push_front(Point(x,y,z));
	}

	return Triangulation(L.begin(), L.end());
}

// call <number of vertices>
int main(int argn, char *args[])
{
	int n = 10;

	if (argn >= 2)
		n = atoi(args[1]);

	Triangulation T = get_rand_triangulation(n);
	Triangulation2 T2 = T;

	assert( T.is_valid() );

	CGAL::Geomview_stream gs;
	gs.set_wired(true);
	gs << T;

	//std::cout << T.dimension() << ' ' << T.number_of_vertices() << ' ' << T.number_of_cells() << std::endl;
	//std::cout << "Number of finite edges: " << T.number_of_finite_edges() << std::endl;

	//for (Finite_edges_iterator it = T.finite_edges_begin(); it != T.finite_edges_end(); it++) {
	//	std::cout << it->first->vertex(it->third)->point() << std::endl;
	//}

	//std::ofstream oFileT("output",std::ios::out);
	//oFileT << T;

	Finite_edges_iterator eit = T.finite_edges_begin();
	for (; eit != T.finite_edges_end(); eit++) {
		Cell_handle cell = eit->first;
		Vertex_handle v1 = cell->vertex( eit->second );
		Vertex_handle v2 = cell->vertex( eit->third );

		bool res = T.is_top_collapsible(*eit);
		if (res) cout << "Yahoo!" << endl;
		else cout << ".";

		//if ((v1->point() == orig && v2->point() == last) || 
		//    (v2->point() == orig && v1->point() == last))
		//	cout << "T.is_collapsible(Edge( " << v1->point() << " -- " << v2->point() << ")):       " << T.is_collapsible(*eit) << endl;
		//T.collapse_edge(*eit);
	}

	char ch;
	cout << "Enter something";
	cin >> ch;

	return 0;
}
