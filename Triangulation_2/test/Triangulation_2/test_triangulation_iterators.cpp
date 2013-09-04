// #define CGAL_PROFILE
// #define CGAL_NO_ASSERTIONS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <boost/foreach.hpp>

//#define LOG {std::cout << "L:" << __LINE__ << std::endl;}
#define LOG {}

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_2<K>                            T;

int main(int argc, char *argv[]) {
	T t;
	const T &t_const = t;
	
	t.insert(T::Point(0,0));
	t.insert(T::Point(1,0));
	t.insert(T::Point(0,1));
	t.insert(T::Point(1,2));

	{	// Vertices
		size_t nVertices = t.number_of_vertices();
		size_t cnt = 0;
		// Iterator over Finite_vertices
		for (T::Finite_vertices_iterator vit = t.finite_vertices_begin(); vit != t.finite_vertices_end(); ++vit) { ++cnt; }
		for (T::Finite_vertices_iterator vit = t.finite_vertices().begin(); vit != t.finite_vertices().end(); ++vit) { ++cnt; }
		T::Finite_vertices_range fvr = t.finite_vertices();
		BOOST_FOREACH(T::Vertex &v, fvr) { ++cnt; }
		BOOST_FOREACH(T::Vertex &v, t.finite_vertices()) { ++cnt; }

		// Const iterator over Finite_vertices
		for (T::Finite_vertices_iterator vit = t_const.finite_vertices_begin(); vit != t_const.finite_vertices_end(); ++vit) { ++cnt; }
		for (T::Finite_vertices_iterator vit = t_const.finite_vertices().begin(); vit != t_const.finite_vertices().end(); ++vit) { ++cnt; }
		T::Finite_vertices_range fvr_c = t_const.finite_vertices();
		BOOST_FOREACH(T::Vertex &v, fvr_c) { ++cnt; }
		BOOST_FOREACH(T::Vertex &v, t_const.finite_vertices()) { ++cnt; }

		// Iterator over Finite_vertex_handles
		T::Finite_vertex_handles_range vhs = t.finite_vertex_handles();
		BOOST_FOREACH(T::Finite_vertices_iterator &vh, vhs) { ++cnt; }
		BOOST_FOREACH(T::Finite_vertices_iterator &vh, t.finite_vertex_handles()) { ++cnt; }
		BOOST_FOREACH(T::Vertex_handle vh, t.finite_vertex_handles()) { ++cnt; }

		// Const iterator over Finite_vertex_handles
		T::Finite_vertex_handles_range vhs_c = t_const.finite_vertex_handles();
		BOOST_FOREACH(T::Finite_vertices_iterator &vh, vhs_c) { ++cnt; }
		BOOST_FOREACH(T::Finite_vertices_iterator &vh, t_const.finite_vertex_handles()) { ++cnt; }
		BOOST_FOREACH(T::Vertex_handle vh, t_const.finite_vertex_handles()) { ++cnt; }

		CGAL_assertion(cnt == nVertices * 14);
	}

	{	// Points
		size_t nVertices = t.number_of_vertices();
		size_t cnt = 0;

		// Iterator range over the points in the triangulation
		BOOST_FOREACH(T::Point &p, t.points()) { ++cnt; }
		BOOST_FOREACH(T::Point p, t.points()) { ++cnt; }

		CGAL_assertion(cnt == nVertices * 2);
	}

	return 0;
}
