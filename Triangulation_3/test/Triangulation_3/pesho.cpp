#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/utility.h>

#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
//#include <CGAL/Periodic_3_triangulation_traits_3.h>
//#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>

//#include <CGAL/collapse.h>

#include <CGAL/IO/Color.h>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

//#include <CORE/CORE.h>
//#include <CGAL/CORE_Expr.h> // From CORE 1.4.1
#include <CGAL/Cartesian.h>

//#include <QGLViewer/qglviewer.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

#include <CGAL/squared_distance_3.h>
#include <CGAL/determinant.h>
#include <CGAL/internal/Static_filters/Orientation_3.h>

#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <list>
#include <set>
#include <queue>

using namespace std;

//#define PESHO_DEBUG

//typedef CORE::Expr							NT;
//typedef CGAL::Cartesian<NT>   						K;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel		K;
typedef CGAL::Exact_predicates_exact_constructions_kernel		K;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K>		traits;

typedef CGAL::Triangulation_vertex_base_with_info_3<int, K>		Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>			Tds;
//typedef CGAL::Triangulation_3<K, Tds>					Triangulation;
//typedef CGAL::Delaunay_triangulation_3<K, Tds>				Triangulation;
//typedef CGAL::Regular_triangulation_3<traits>				Triangulation;
//typedef CGAL::Periodic_3_triangulation_traits_3<K>			GT;
//typedef CGAL::Periodic_3_Delaunay_triangulation_3<GT>			Triangulation;

typedef K::FT FT;
typedef K::Vector_3					Vector;

typedef Triangulation::Finite_vertices_iterator 	Finite_vertices_iterator;
typedef Triangulation::Finite_edges_iterator 		Finite_edges_iterator;
typedef Triangulation::Finite_facets_iterator 		Finite_facets_iterator;
typedef Triangulation::Finite_cells_iterator 		Finite_cells_iterator;
typedef Triangulation::Face_circulator 			Face_circulator;
typedef Triangulation::Cell_circulator 			Cell_circulator;
typedef Triangulation::Locate_type    			Locate_type;
typedef Triangulation::Point          			Point;
typedef Triangulation::Triangle          		Triangle;
typedef Triangulation::Tetrahedron			Tetrahedron;

typedef Triangulation::Vertex_handle			Vertex_handle;
typedef Triangulation::Cell_handle			Cell_handle;
typedef Triangulation::Facet				Facet;
typedef Triangulation::Cell				Cell;
typedef Triangulation::Edge				Edge;

typedef CGAL::Creator_uniform_3<double,Point>  		Creator;
typedef CGAL::Geomview_stream				Geomview_stream;

struct Simple_edge {
	Cell_handle first;
	Vertex_handle second, third;
	FT len;
	Simple_edge(Cell_handle _first, Vertex_handle _v1, Vertex_handle _v2) {
		first = _first;
		second = _v1;
		third = _v2;
		FT len = squared_distance(second->point(), third->point());
	}
};

bool operator<(const Simple_edge &a, const Simple_edge &b)
{
	if (a.len != b.len) return a.len > b.len;
	if (a.second != b.second) return a.second < b.second;
	return a.third < b.third;
}

typedef priority_queue< Simple_edge, vector<Simple_edge> >	PQ_of_edges;

CGAL::Random my_rand(40);

std::list<Point> get_my_points()
{
	std::list<Point> L;
	
	L.push_front(Point(0,2,0));
	L.push_front(Point(0,0,1));
	L.push_front(Point(0,1,1));
	L.push_front(Point(0,0,0));
	L.push_front(Point(1,1,1));

	return L;
}

std::list<Point> get_layered_2d_points(int n, int m)
{
	std::list<Point> L;

	for(int i=0; i<n; i++) {
		for(int j=0; j<m; j++) {
			double x = my_rand.get_double(0.0, 5.0);
			double y = x;
			double z = i;

			L.push_front(Point(x,y,z));
		}
	}

	return L;
}

std::list<Point> get_layered_3d_points(int n, int m)
{
	std::list<Point> L;

	for(int i=0; i<n; i++) {
		for(int j=0; j<m; j++) {
			double x = my_rand.get_double(0.0, 5.0);
			double y = my_rand.get_double(0.0, 5.0);
			double z = i;

			L.push_front(Point(x,y,z));
		}
	}

	return L;
}

std::list<Point> get_my_2d_points()
{
	std::list<Point> L;
	
	L.push_front(Point(0,0,0));
	L.push_front(Point(0,2,0));
	L.push_front(Point(0,0,3));
	L.push_front(Point(0,1,1));

	return L;
}

std::list<Point> get_1D(int n)
{
	int i;
	std::list<Point> L;

	for (i=0; i<n; i++) {
		double x = my_rand.get_double(0.0, 1.0);
		double y = x;
		double z = 0.0; //my_rand.get_double(0.0, 1.0);

		L.push_front(Point(x,y,z));
	}

	return L;
}

std::list<Point> get_2D(int n)
{
	int i;
	std::list<Point> L;

	for (i=0; i<n; i++) {
		double x = my_rand.get_double(0.0, 1.0);
		double y = my_rand.get_double(0.0, 1.0);
		double z = 0.0; //my_rand.get_double(0.0, 1.0);

		L.push_front(Point(x,y,z));
	}

	return L;
}

list<Point> get_rand_in_sphere(int n, int r)
{
	list<Point> L;

	CGAL::Random_points_in_sphere_3<Point,Creator> g(r, my_rand);
	CGAL::copy_n(g, n, std::back_inserter(L));

	return L;
}

///////////////
// VISUALIZE //
///////////////

void geomview_show_vertices(Geomview_stream& gs, Triangulation& T)
{
	Finite_vertices_iterator vit;
	for(vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); vit++)
		gs << vit->point();
}

void geomview_show_edge(Geomview_stream &gs, Vertex_handle v1, Vertex_handle v2)
{
	gs.set_edge_color( CGAL::RED );
	gs.set_line_width(3);
	CGAL::Segment_3<K> sg(v1->point(), v2->point());
	gs << sg;
	gs.set_line_width(1);
	
	//double r = gs.set_vertex_radius(1.0);
	//gs.set_vertex_color( CGAL::YELLOW );
	gs << v2->point();	
	//gs.set_vertex_color( CGAL::GREEN );
}

void geomview_show_point (Geomview_stream &gs, double x, double y, double z, CGAL::Color color)
{
	Point p(x,y,z);
	gs.set_vertex_color(color);
	gs << p;	
}

void geomview_show_cell (Geomview_stream &gs, Point &p1, Point &p2, Point &p3, Point &p4, CGAL::Color color)
{
	gs.set_face_color(color);
	Tetrahedron t(p1,p2,p3,p4);
	gs << t;	
}

void write_to_OFF(const char *filename, Triangulation &T)
{
	ofstream oFile(filename, std::ios::out);

	// Write polyhedron in Object File Format (OFF).
	CGAL::set_ascii_mode( oFile );
	oFile << "OFF" << endl << T.number_of_vertices() << " " << T.number_of_finite_facets() << " 0" << endl;
	copy( T.points_begin(), T.points_end(), ostream_iterator<Point>( oFile, "\n"));

	cout << "Sizes: " << std::distance(T.points_begin(), T.points_end()) << " " << std::distance(T.finite_vertices_begin(), T.finite_vertices_end()) << endl;

	for(Finite_facets_iterator fit = T.finite_facets_begin(); fit != T.finite_facets_end(); ++fit) {
		Cell_handle cell = fit->first;
		int findex = fit->second;
		
		Vertex_handle v[] = { cell->vertex(0), cell->vertex(1), cell->vertex(2), cell->vertex(3) };
		swap(v[3], v[findex]);		

		oFile << '4';
		oFile << ' ' << std::distance(T.vertices_begin(), v[0]) - 1;
		oFile << ' ' << std::distance(T.vertices_begin(), v[1]) - 1;
		oFile << ' ' << std::distance(T.vertices_begin(), v[2]) - 1;
		oFile << endl;
	}
}

void push_edge(PQ_of_edges &pq, Edge &e)
{
	Cell_handle cell = e.first;
	Vertex_handle v1 = cell->vertex( e.second );
	Vertex_handle v2 = cell->vertex( e.third );
	//pq.push( make_triple(cell, v1, v2) );
	//pq.push( make_triple(cell, v2, v1) );
	pq.push( Simple_edge(cell, v1, v2) );
	pq.push( Simple_edge(cell, v2, v1) );
}

void insert_all_edges_to_queue(Triangulation &T, PQ_of_edges &pq)
{
	Triangulation::All_edges_iterator eit;	
	for (eit=T.all_edges_begin(); eit!=T.all_edges_end(); eit++) {
		Edge e = *eit;
		if (T.is_infinite(e))
			push_edge(pq, e);
	}
}

void insert_finite_edges_to_queue(Triangulation &T, PQ_of_edges &pq)
{
	Finite_edges_iterator eit;
	for (eit=T.finite_edges_begin(); eit != T.finite_edges_end(); eit++) {

		Edge e = *eit;
			push_edge(pq, e);
	}
}

void insert_adjacent_edges_to_queue(Triangulation &T, PQ_of_edges &pq, Vertex_handle v)
{
	switch( T.dimension() ) {
	case 1: {
		vector<Edge> edges;
		T.incident_edges (v, back_inserter(edges));

		for(vector<Edge>::iterator eit=edges.begin(); eit!=edges.end(); eit++) {
			Edge e = *eit;
			if(!T.is_infinite(e)) {
				push_edge(pq,e);
			}
		}

		return;
	}
	case 2: {
		Face_circulator fcirc = T.tds().incident_faces(v);
		Face_circulator fend = fcirc;

		CGAL_For_all(fcirc,fend) {
			for(int i=0; i<3; i++) {
				Edge e(fcirc, i, (i+1)%3);
				if(!T.is_infinite(e)) {
					push_edge(pq,e);
				}
			}	
		}

		return;
	}
	case 3:
		vector<Cell_handle> cells;
		T.incident_cells (v, back_inserter(cells));

		for(vector<Cell_handle>::iterator cit=cells.begin(); cit!=cells.end(); cit++) {
			for(int i=0; i<3; i++)
				for(int j=i+1; j<3; j++) {
					Edge e(*cit, i, j);
					if(!T.is_infinite(e)) {
						push_edge(pq,e);
					}
				}
		}

		return;
	}

}

void insert_edges_from_adjacent_cells_to_queue(Triangulation &T, PQ_of_edges &pq, Vertex_handle v)
{

}

int number_of_incident_cells(Triangulation &T, Edge e)
{
	int cnt = 0;
	Cell_circulator ccirc = T.incident_cells(e); 
	Cell_circulator cend  = ccirc;

	CGAL_For_all(ccirc, cend)
	{
		cnt++;
	}

	return cnt;
}

bool incident_to_inf(Triangulation &T, Vertex_handle v)
{
	vector<Edge> E;
	T.incident_edges(v, back_inserter(E));

	for(vector<Edge>::iterator eit = E.begin(); eit != E.end(); eit++)
		if (T.is_infinite(*eit))
			return true;

	return false;
}

void print_info(Triangulation &T)
{
	cout << "vertices: " << T.number_of_vertices() << endl;
	cout << "cells: " << T.number_of_cells() << endl;
	cout << "facets: " << T.number_of_facets() << endl;
	cout << "edges: " << T.number_of_edges() << endl;
	cout << "dim: " << T.dimension() << endl;
}

int collapse_all_edges(Geomview_stream &gs1, Geomview_stream &gs2, Triangulation &T, int n)
{
	int stat_start = 0;
	int stat_examined = 0;
	int stat_top_collapsible = 0;
	int stat_geom_collapsible = 0;
	int stat_collapsed = 0;
	int stat_collapsed_internal = 0;
	int stat_collapsed_inf = 0;
	int cnt = 0;

	print_info(T);
	
	gs1.clear();
	gs2.clear();
	gs1.set_vertex_color( CGAL::GREEN );
	gs2.set_vertex_color( CGAL::GREEN );
	gs1.set_wired(true);
	gs2.set_wired(true);
	
	PQ_of_edges pq;
	insert_finite_edges_to_queue(T, pq);

	stat_start = pq.size();
	T.is_valid(true);

	gs1 << T;
	//geomview_show_vertices(gs1, T);

	while(T.dimension() >= 1 && !pq.empty()) {// && T.number_of_vertices()>4) {
		Simple_edge se = pq.top();
		pq.pop();

		Cell_handle c = se.first;
		Vertex_handle v1 = se.second;
		Vertex_handle v2 = se.third;;
		int v1_index, v2_index;

	cout << ".0";
		if(!T.is_simplex(c)) continue;

		if(!c->has_vertex(v1, v1_index)) continue;
		if(!c->has_vertex(v2, v2_index)) continue;

		//Point p1 = v1->point();
		//Point p2 = v2->point();
		//Vector segm(p1,p2);
		//Point mid = p2;//p1 + segm*0.9;

		//gs1 << mid;
	cout << ".1";
		Edge e( c, v1_index, v2_index );
		bool top =  T.tds().is_collapsible(e);

	cout << ".2";
		bool geom = T.is_collapsible(e);

	cout << ".3";
		stat_examined ++;
		stat_top_collapsible += top;
		stat_geom_collapsible += geom;

		//cout << "top: " << top << " geom: " << geom << endl;

		//geomview_show_edge(gs1,v1,v2);

		if (top && geom) {
			//cout << "cnt: " << cnt << " n: " << n << endl;
			//if (cnt++ < n) continue;
			
			cout << "before collapse:";
			T.is_valid(true);

			stat_collapsed++;
			if (!incident_to_inf(T,v2)) stat_collapsed_internal++;
			if (T.is_infinite(e)) stat_collapsed_inf++;	
	cout << ".4";
			geomview_show_edge(gs1,v1,v2);

			cout << "vertices: " << T.number_of_vertices() << endl;

			//if (T.number_of_vertices() == 4) T.tds().remove_from_maximal_dimension_simplex(v1);
			//else
			if (!T.collapse(e)) cout << "Hm.. Strange!" << endl;
			cout << "collapsed dim=" << T.dimension() << endl;

	cout << ".5";
			T.is_valid(true);
	cout << ".6";
			//gs2.clear();
			//gs2 << T;
			//gs2 << v2->point();
			//cin.get();

			insert_adjacent_edges_to_queue(T, pq, v2);
		}
	}

	gs2 << T;
	geomview_show_vertices(gs2, T);

	print_info(T);

	cout << "start: " << stat_start << endl;
	cout << "examined: " << stat_examined << endl;
	cout << "edges_examined: " << stat_examined << endl;
	cout << "top_collapsible: " << stat_top_collapsible << endl;
	cout << "geom_collapsible: " << stat_geom_collapsible << endl;
	cout << "collapsed: " << stat_collapsed << endl;
	cout << "collapsed_internal: " << stat_collapsed_internal << endl;
	cout << "collapsed_inf: " << stat_collapsed_inf << endl;

	return stat_collapsed;
}

// call <number of vertices>
int main(int argn, char *args[])
{
	int n = 10;
	double r = 1.0;
	Geomview_stream gs1, gs2;//, gs3, gs4;
	
	if (argn >= 2)
		n = atoi(args[1]);

	list<Point> L;

	//L = get_layered_2d_points(n, 10*n);	
	//L = get_layered_3d_points(n, 10*n);	
	L = get_rand_in_sphere(n, r);
	//L = get_my_points();
	//L = get_my_2d_points();
	//L = get_2D(n);
	//L = get_1D(n);
	Triangulation T(L.begin(), L.end());

	gs1 << T;

	T.is_valid(true);
	/*
	Triangulation::All_vertices_iterator vit;
	for(vit=T.all_vertices_begin(); vit!=T.all_vertices_end(); vit++)
		if (!n--) {
			Vertex_handle v = vit;
			T.tds().remove_from_maximal_dimension_simplex(v);
			cout << T.is_valid(true) << endl;
			T.tds().reorient();
			cout << T.is_valid(true) << endl;
			break;
		}
	*/

	//write_to_OFF("init.off", T);
	//std::ofstream oFileT("output",std::ios::out);
	//oFileT << T;

	//while (collapse_all_edges(gs1, gs2, T));
	collapse_all_edges(gs1, gs2, T, n);
	
	cin.get();
	
	return 0;
}
