#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/utility.h>

#include <CGAL/collapse.h>

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

typedef CGAL::Triangulation_vertex_base_with_info_3<CGAL::Color, K>	Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>			Tds;
typedef CGAL::Triangulation_3<K, Tds>					Triangulation;
//typedef CGAL::Delaunay_triangulation_3<K, Tds>				Triangulation;

typedef K::FT FT;
typedef K::Vector_3					Vector;

typedef Triangulation::Finite_vertices_iterator 	Finite_vertices_iterator;
typedef Triangulation::Finite_edges_iterator 		Finite_edges_iterator;
typedef Triangulation::Finite_facets_iterator 		Finite_facets_iterator;
typedef Triangulation::Finite_cells_iterator 		Finite_cells_iterator;
typedef Triangulation::Cell_circulator 			Cell_circulator;
typedef Triangulation::Simplex        			Simplex;
typedef Triangulation::Locate_type    			Locate_type;
typedef Triangulation::Point          			Point;
typedef Triangulation::Triangle          		Triangle;
typedef Triangulation::Tetrahedron			Tetrahedron;

typedef Triangulation::Vertex_handle			Vertex_handle;
typedef Triangulation::Cell_handle			Cell_handle;
typedef Triangulation::Facet				Facet;
typedef Triangulation::Edge				Edge;

typedef CGAL::Creator_uniform_3<double,Point>  		Creator;
typedef CGAL::Geomview_stream				Geomview_stream;

typedef CGAL::Triple<Cell_handle, Vertex_handle, Vertex_handle>	Simple_edge;

// to make to choice better
struct my_edge_cmp {
	bool operator()(const Simple_edge &a, const Simple_edge &b)
	{
		return a < b;

		Vertex_handle av1 = a.second;
		Vertex_handle av2 = a.third;
		Vertex_handle bv1 = b.second;
		Vertex_handle bv2 = b.third;

		if (av1==bv1 && av2==bv2)
			return 0;

		FT da = squared_distance(av1->point(), av2->point());
		FT db = squared_distance(bv1->point(), bv2->point());
		
		if (da==db)
			return a < b;
		return da > db;
	}
};

typedef priority_queue< Simple_edge, vector<Simple_edge>, my_edge_cmp >	PQ_of_edges;

CGAL::Random my_rand(42);

std::list<Point> get_my_points()
{
	std::list<Point> L;
	
	L.push_front(Point(0,2,0));
	L.push_front(Point(0,0,1));
	L.push_front(Point(0,1,1));
	L.push_front(Point(0,0,0));
	L.push_front(Point(0,0.5,.5));

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
	
	gs << v2->point();
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
	pq.push( make_triple(cell, v1, v2) );
	pq.push( make_triple(cell, v2, v1) );
}

void insert_edges_to_queue(Triangulation &T, PQ_of_edges &pq)
{
	Finite_edges_iterator eit;
	
	for (eit=T.finite_edges_begin(); eit != T.finite_edges_end(); eit++) {
		Edge e = *eit;
		push_edge(pq, e);
	}
}

//void insert_adjacent_edges_to_queue(Triangulation &T, PQ_of_edges &pq, Vertex_handle v)
//{
//	vector<Edge> E;
//	T.incident_edges(v, back_inserter(E)); 

	//cout << "incident edges of (" << v->point() << "): " << E.size() << endl;
//	for(vector<Edge>::iterator eit = E.begin(); eit != E.end(); eit++) {
//		Edge e = *eit;
//		if(!T.is_infinite(e)) {
//			push_edge(pq, e);
//		}
//	}
//}

void insert_adjacent_edges_to_queue(Triangulation &T, PQ_of_edges &pq, Vertex_handle v)
{
	vector<Edge> edges;
	T.incident_edges (v, back_inserter(edges));

	for(vector<Edge>::iterator eit=edges.begin(); eit!=edges.end(); eit++) {
		Edge e = *eit;
		if(!T.is_infinite(e)) {
			push_edge(pq,e);
		}
	}
}

void insert_edges_from_adjacent_cells_to_queue(Triangulation &T, PQ_of_edges &pq, Vertex_handle v)
{
	vector<Cell_handle> cells;
	T.incident_cells (v, back_inserter(cells));

	//cout << "dim: " << T.dimension() << endl;
	//cout << "number of cells: " << T.number_of_cells() << endl;
	//cout << "incident cells: " << cells.size() << endl;

	for(vector<Cell_handle>::iterator cit=cells.begin(); cit!=cells.end(); cit++) {
		for(int i=0; i<3; i++)
			for(int j=i+1; j<3; j++) {
				Edge e(*cit, i, j);
				if(!T.is_infinite(e)) {
					push_edge(pq,e);
				}
			}
	}
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

int collapse_all_edges(Geomview_stream &gs1, Geomview_stream &gs2, Triangulation &T)
{
	int stat_start = 0;
	int stat_examined = 0;
	int stat_top_collapsible = 0;
	int stat_geom_collapsible = 0;
	int stat_collapsed = 0;
	int stat_collapsed_target_incident_to_inf = 0;

	cout << "vertices: " << T.number_of_vertices() << endl;
	cout << "cells: " << T.number_of_cells() << endl;
	cout << "facets: " << T.number_of_facets() << endl;
	cout << "edges: " << T.number_of_edges() << endl;
	cout << "dim: " << T.dimension() << endl;
	
	gs1.clear();
	gs2.clear();
	gs1.set_vertex_color( CGAL::GREEN );
	gs2.set_vertex_color( CGAL::GREEN );
	gs1.set_wired(true);
	gs2.set_wired(true);
	
	PQ_of_edges pq;
	insert_edges_to_queue(T, pq);

	stat_start = pq.size();
	T.is_valid(true);

	gs1 << T;
	geomview_show_vertices(gs1, T);

	while(!pq.empty() && T.number_of_vertices()>4) {
		Simple_edge se = pq.top();
		pq.pop();

		Cell_handle c = se.first;
		Vertex_handle v1 = se.second;
		Vertex_handle v2 = se.third;;
		int v1_index, v2_index;

	cout << ".-1";
		//if (!T.is_cell(c)) continue;
		if(!T.is_simplex(c)) continue; // is_cell?

	cout << ".0";
		if(!c->has_vertex(v1, v1_index)) continue;
		if(!c->has_vertex(v2, v2_index)) continue;

		Point p1 = v1->point();
		Point p2 = v2->point();
		Vector segm(p1,p2);
		Point mid = p2;// + segm/2;

		//gs1 << mid;
	cout << ".1";
		Edge e( c, v1_index, v2_index );
		bool top = T.is_top_collapsible(e);

	cout << ".2";
		bool geom = T.is_geom_collapsible(e);

	cout << ".3";
		stat_examined ++;
		stat_top_collapsible += top;
		stat_geom_collapsible += geom;

		if (top && geom) {
			//geomview_show_edge(gs1,v1,v2);
			T.collapse_edge(e);

	cout << ".4";
			//if (!T.is_valid(false)) {
			//	cout << "NOT VALID TRIANGULATION!" << endl;
			//	cout << "Dim: " << T.dimension() << endl;
			//	cout << "Vertices: " << T.number_of_vertices() << endl;
			//}

			stat_collapsed++;
			if (incident_to_inf(T,v2))
				stat_collapsed_target_incident_to_inf++;

			//gs2.clear();
			//geomview_show(gs2, T);
			//gs2 << v2->point();
			//cin.get();

			switch( T.dimension() ) {
			case 1:
				insert_adjacent_edges_to_queue(T, pq, v2);
				break;
			case 3:
				insert_edges_from_adjacent_cells_to_queue(T, pq, v2);
				break;
			}
		}
	}

	gs2 << T;
	geomview_show_vertices(gs2, T);

	cout << "start: " << stat_start << endl;
	cout << "examined: " << stat_examined << endl;
	cout << "edges_examined: " << stat_examined << endl;
	cout << "top_collapsible: " << stat_top_collapsible << endl;
	cout << "geom_collapsible: " << stat_geom_collapsible << endl;
	cout << "collapsed: " << stat_collapsed << endl;
	cout << "collapsed_target_incident_to_inf: " << stat_collapsed_target_incident_to_inf << endl;

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

	L = get_rand_in_sphere(n, r);
	//L = get_my_points();
	//L = get_rand_2D(n);
	//L = get_1D(n);
	Triangulation T(L.begin(), L.end());

	//write_to_OFF("init.off", T);
	//std::ofstream oFileT("output",std::ios::out);
	//oFileT << T;

	//while (collapse_all_edges(gs1, gs2, T));
	collapse_all_edges(gs1, gs2, T);
	
	cin.get();
	
	return 0;
}
