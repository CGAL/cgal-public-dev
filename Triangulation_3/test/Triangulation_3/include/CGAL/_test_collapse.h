#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>

#include <vector>
#include <list>
#include <queue>

// simplificated representation of Edge as <cell,vertex_handle,vertex_handle>
// additional information (len) determins the order of the collapse operations
template <typename Triangulation> 
struct Simple_edge {
  typedef typename Triangulation::Geom_traits::FT   FT;  
  typedef typename Triangulation::Vertex_handle     Vertex_handle;
  typedef typename Triangulation::Cell_handle       Cell_handle;

  Cell_handle first;
  Vertex_handle second, third;
  FT len;
  
  Simple_edge(Cell_handle _first, Vertex_handle _v1, Vertex_handle _v2) {
    first = _first;
    second = _v1;
    third = _v2;
    len = squared_distance(second->point(), third->point());
  }

  bool operator<(const Simple_edge &b) const
  {
    if (this->len != b.len) return this->len > b.len;
    if (this->second != b.second) return this->second < b.second;
    return this->third < b.third;
  }
};

// TODO: remove the function
template <typename Triangulation>
void print_info(Triangulation &T)
{
  std::cout << "vertices: " << T.number_of_vertices() << std::endl;
  std::cout << "cells: " << T.number_of_cells() << std::endl;
  std::cout << "facets: " << T.number_of_facets() << std::endl;
  std::cout << "edges: " << T.number_of_edges() << std::endl;
  std::cout << "dim: " << T.dimension() << std::endl;
}

template <typename Point> 
std::list<Point> get_layered_2d_points(int n, int m)
{
	std::list<Point> L;

        CGAL::Random my_rand(0);
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

template <typename Point> 
std::list<Point> get_layered_3d_points(int n, int m)
{
	std::list<Point> L;

        CGAL::Random my_rand(0);
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

template <typename Point> 
std::list<Point> get_my_2d_points()
{
	std::list<Point> L;
	
	L.push_front(Point(0,0,0));
	L.push_front(Point(0,2,0));
	L.push_front(Point(0,0,3));
	L.push_front(Point(0,1,1));

	return L;
}

template <typename Point> 
std::list<Point> get_rand_in_line(int n)
{
	int i;
	std::list<Point> L;

        CGAL::Random my_rand(0);
	for (i=0; i<n; i++) {
		double x = my_rand.get_double(0.0, 1.0);
		double y = x;
		double z = 0.0; //my_rand.get_double(0.0, 1.0);

		L.push_front(Point(x,y,z));
	}

	return L;
}

template <typename Point> 
std::list<Point> get_rand_in_disc(int n)
{
	int i;
	std::list<Point> L;

        CGAL::Random my_rand(0);
	for (i=0; i<n; i++) {
		double x = my_rand.get_double(0.0, 1.0);
		double y = my_rand.get_double(0.0, 1.0);
		double z = 0.0; //my_rand.get_double(0.0, 1.0);

		L.push_front(Point(x,y,z));
	}

	return L;
}

template <typename Point> 
std::list<Point> get_rand_in_sphere(int n, int r)
{
  std::list<Point> L;

  // sets fixed seed
  CGAL::Random my_rand(0);
  CGAL::Random_points_in_sphere_3<Point> g(r, my_rand);
  CGAL::copy_n(g, n, std::back_inserter(L));

  return L;
}

template<typename Triangulation, typename PQ_of_edges>
struct push_edge_t {
  static void push_edge(Triangulation &T, PQ_of_edges &pq, const typename Triangulation::Edge &e) 
  {
    typename Triangulation::Cell_handle cell = e.first;
    typename Triangulation::Vertex_handle v1 = cell->vertex( e.second );
    typename Triangulation::Vertex_handle v2 = cell->vertex( e.third );

    if (!T.is_infinite(e)) {
      // PQ_of_edges::values_type stands for Simple_edge
      pq.push( typename PQ_of_edges::value_type(cell, v1, v2) );
      pq.push( typename PQ_of_edges::value_type(cell, v2, v1) );
    }
  }
};

// inserts all the edges incident to a vertex v to the priority queue
template<typename Triangulation, typename PQ_of_edges>
void insert_adjacent_edges_to_queue(Triangulation &T, PQ_of_edges &pq, typename Triangulation::Vertex_handle v)
{
  typedef typename Triangulation::Edge                Edge;
  typedef typename Triangulation::Cell_handle         Cell_handle;
  typedef typename Triangulation::Face_circulator     Face_circulator;

  switch( T.dimension() ) {
  case 1: {
    std::vector<Edge> edges;
    T.incident_edges (v, back_inserter(edges));

    for(typename std::vector<Edge>::iterator eit=edges.begin(); eit!=edges.end(); eit++) {
      Edge e = *eit;
      push_edge_t<Triangulation,PQ_of_edges>::push_edge(T,pq,e);
    }

    return;
  }
  case 2: {
    Face_circulator fcirc = T.tds().incident_faces(v);
    Face_circulator fend = fcirc;

    CGAL_For_all(fcirc,fend) {
      for(int i=0; i<3; i++) {
        Edge e(fcirc, i, (i+1)%3);
        push_edge_t<Triangulation,PQ_of_edges>::push_edge(T,pq,e);
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
          push_edge_t<Triangulation,PQ_of_edges>::push_edge(T,pq,e);
        }
    }

    return;
  }
}

// generates a set of random points and iteratively tries to collapse the remaining edges
// until the triangulation becames a single finite point and no more collapses are possible
template <typename Triangulation>
void collapse_to_a_point(Triangulation T)
{
  typedef typename Triangulation::Finite_edges_iterator     Finite_edges_iterator;
  typedef typename Triangulation::Vertex_handle             Vertex_handle;
  typedef typename Triangulation::Cell_handle               Cell_handle;

  typedef std::priority_queue< Simple_edge<Triangulation> >  PQ_of_edges;

  assert(T.is_valid());

  int stat_collapsed = 0;
  PQ_of_edges pq;

  // insert all the finite edges (reversed too) to priority queue
  Finite_edges_iterator eit;
  for (eit=T.finite_edges_begin(); eit != T.finite_edges_end(); eit++) 
    push_edge_t<Triangulation,PQ_of_edges>::push_edge(T, pq, *eit);

  // until a single finite and an infinite point remains
  while (!pq.empty()) {
    Simple_edge<Triangulation> se = pq.top();
    pq.pop();

    Cell_handle c = se.first;
    Vertex_handle v1 = se.second;
    Vertex_handle v2 = se.third;;
    int v1_index, v2_index;

    // the corresponding simplex and vertices could have
    // disappeared due to other neighbour collapses
    if (!T.is_simplex(c)) continue;
    if (!c->has_vertex(v1, v1_index)) continue;
    if (!c->has_vertex(v2, v2_index)) continue;

    typename Triangulation::Edge e( c, v1_index, v2_index );
    if (T.is_collapsible(e)) {
      stat_collapsed++;

      assert(T.collapse(e));
      assert(T.is_valid());

      // prepare the neighbouring edges to be collapsed
      insert_adjacent_edges_to_queue<Triangulation,PQ_of_edges>(T, pq, v2);
    }
  }

  assert(T.number_of_vertices() == 1);
  assert(T.dimension() == 0);
}

template <typename Triangulation>
void _test_collapse() {
  typedef typename Triangulation::Point                     Point;

  int n = 20;

  std::list<Point> L;

  L = get_rand_in_line<Point>(n);
  collapse_to_a_point<Triangulation>( Triangulation(L.begin(), L.end()) );
  std::cout << "        1D triangulation of " << n << " random points in a line collapsed to a single point." << std::endl;
  
  L = get_rand_in_disc<Point>(n);
  collapse_to_a_point<Triangulation>( Triangulation(L.begin(), L.end()) );
  std::cout << "        2D triantulation of " << n << " random points in a disc collapsed to a single point." << std::endl;

  L = get_rand_in_sphere<Point>(n, 1.0);
  collapse_to_a_point<Triangulation>( Triangulation(L.begin(), L.end()) );
  std::cout << "        3D triangulation of " << n << " random points in a sphere collapsed to a single point." << std::endl;

  /*
  L = get_layered_2d_points<Point>(3, n); 
  collapse_to_a_point<Triangulation>(L);
  std::cout << "        Layered 2d points collapsed." << std::endl;

  L = get_layered_3d_points<Point>(3, n);  
  collapse_to_a_point<Triangulation>(L);
  std::cout << "        Layered 3d points collapsed." << std::endl;
  */
  //L = get_my_points();
  //L = get_my_2d_points();
}

