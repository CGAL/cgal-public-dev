#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

#include <vector>
#include <list>
#include <queue>

template <typename Tri>
void print_info(Tri &T)
{
	std::cout << "vertices: " << T.number_of_vertices() << std::endl;
	std::cout << "cells: " << T.number_of_cells() << std::endl;
	std::cout << "facets: " << T.number_of_facets() << std::endl;
	std::cout << "edges: " << T.number_of_edges() << std::endl;
	std::cout << "dim: " << T.dimension() << std::endl;
}

template <typename Point> 
std::list<Point> get_rand_in_sphere(int n, int r)
{
	std::list<Point> L;

	CGAL::Random_points_in_sphere_3<Point> g(r);
	CGAL::copy_n(g, n, std::back_inserter(L));

	return L;
}

template<typename Tri, typename PQ_of_edges>
void push_edge(PQ_of_edges &pq, const typename Tri::Edge &e) 
{
	typename Tri::Cell_handle cell = e.first;
	typename Tri::Vertex_handle v1 = cell->vertex( e.second );
	typename Tri::Vertex_handle v2 = cell->vertex( e.third );

	// Simple_edge
	pq.push( typename PQ_of_edges::value_type(cell, v1, v2) );
	pq.push( typename PQ_of_edges::value_type(cell, v2, v1) );
}

// simplificated representation of Edge as <cell,vertex_handle,vertex_handle>
// additional information (len) determins the order of the collapse operations
template <typename Tri> 
struct Simple_edge {
	typedef typename Tri::Geom_traits::FT	FT;  
	typedef typename Tri::Vertex_handle	Vertex_handle;
	typedef typename Tri::Cell_handle	Cell_handle;

	Cell_handle first;
	Vertex_handle second, third;
	FT len;
	
	Simple_edge(Cell_handle _first, Vertex_handle _v1, Vertex_handle _v2) {
		first = _first;
		second = _v1;
		third = _v2;
		FT len = squared_distance(second->point(), third->point());
	}

	bool operator<(const Simple_edge &b) const
	{
		if (this->len != b.len) return this->len > b.len;
		if (this->second != b.second) return this->second < b.second;
		return this->third < b.third;
	}
};

template<typename Tri, typename PQ_of_edges>
void insert_adjacent_edges_to_queue(Tri &T, PQ_of_edges &pq, typename Tri::Vertex_handle v)
{
	typedef typename Tri::Edge			Edge;
	typedef typename Tri::Cell_handle		Cell_handle;
	typedef typename Tri::Face_circulator 		Face_circulator;

	switch( T.dimension() ) {
	case 1: {
		std::vector<Edge> edges;
		T.incident_edges (v, back_inserter(edges));

		for(typename std::vector<Edge>::iterator eit=edges.begin(); eit!=edges.end(); eit++) {
			Edge e = *eit;
			if(!T.is_infinite(e))
				push_edge<Tri,PQ_of_edges>(pq,e);
		}

		return;
	}
	case 2: {
		Face_circulator fcirc = T.tds().incident_faces(v);
		Face_circulator fend = fcirc;

		CGAL_For_all(fcirc,fend) {
			for(int i=0; i<3; i++) {
				Edge e(fcirc, i, (i+1)%3);
				if(!T.is_infinite(e))
					push_edge<Tri,PQ_of_edges>(pq,e);
			}	
		}

		return;
	}
	case 3:
		std::vector<Cell_handle> cells;
		T.incident_cells (v, back_inserter(cells));

		for(typename std::vector<Cell_handle>::iterator cit=cells.begin(); cit!=cells.end(); cit++) {
			for(int i=0; i<3; i++)
				for(int j=i+1; j<3; j++) {
					Edge e(*cit, i, j);
					if(!T.is_infinite(e))
						push_edge<Tri,PQ_of_edges>(pq,e);
				}
		}

		return;
	}
}

template <typename Tri>
void _test_collapse() {
	typedef typename Tri::Point          			Point;
	typedef typename Tri::Finite_edges_iterator 		Finite_edges_iterator;
	typedef typename Tri::Vertex_handle	Vertex_handle;
	typedef typename Tri::Cell_handle	Cell_handle;

	typedef std::priority_queue< Simple_edge<Tri>, std::vector<Simple_edge<Tri> > >	PQ_of_edges;

	int n = 50;
	std::cout << "      Collapsing a triangulation of random " << n << " vertices in a shpere." << std::endl;

	std::list<Point> L;

	//L = get_layered_2d_points(n, 10*n);	
	//L = get_layered_3d_points(n, 10*n);	
	L = get_rand_in_sphere<Point>(n, 1.0);
	//L = get_my_points();
	//L = get_my_2d_points();
	//L = get_2D(n);
	//L = get_1D(n);

	Tri T(L.begin(), L.end());
	assert(T.is_valid());

	int stat_start = 0;
	int stat_examined = 0;
	int stat_top_collapsible = 0;
	int stat_geom_collapsible = 0;
	int stat_collapsed = 0;
	int stat_collapsed_internal = 0;
	int stat_collapsed_inf = 0;
	int cnt = 0;

	PQ_of_edges pq;

	// insert finite edges to priority queue
	Finite_edges_iterator eit;
	for (eit=T.finite_edges_begin(); eit != T.finite_edges_end(); eit++) 
		push_edge<Tri,PQ_of_edges>(pq, *eit);

	stat_start = pq.size();

	while(!pq.empty()) {
		Simple_edge<Tri> se = pq.top();
		pq.pop();

		Cell_handle c = se.first;
		Vertex_handle v1 = se.second;
		Vertex_handle v2 = se.third;;
		int v1_index, v2_index;

		if(!T.is_simplex(c)) continue;

		if(!c->has_vertex(v1, v1_index)) continue;
		if(!c->has_vertex(v2, v2_index)) continue;

		typename Tri::Edge e( c, v1_index, v2_index );
		if (T.is_collapsible(e)) {
			stat_collapsed++;

			assert(T.collapse(e));
			assert(T.is_valid());

			insert_adjacent_edges_to_queue<Tri,PQ_of_edges>(T, pq, v2);
		}
	}

	assert(T.number_of_vertices() == 1);
	assert(T.dimension() == 0);

	std::cout << "start: " << stat_start << std::endl;
	std::cout << "examined: " << stat_examined << std::endl;
	std::cout << "edges_examined: " << stat_examined << std::endl;
	std::cout << "top_collapsible: " << stat_top_collapsible << std::endl;
	std::cout << "geom_collapsible: " << stat_geom_collapsible << std::endl;
	std::cout << "collapsed: " << stat_collapsed << std::endl;
	std::cout << "collapsed_internal: " << stat_collapsed_internal << std::endl;
	std::cout << "collapsed_inf: " << stat_collapsed_inf << std::endl;
}


