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

	{	// Finite vertices
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

	{	// (Finite) vertices
		size_t nVertices = t.number_of_vertices();
		size_t cnt = 0;
		// Iterator over Finite_vertices
		for (T::Vertex_iterator vit = t.vertices_begin(); vit != t.vertices_end(); ++vit) { ++cnt; }
		for (T::Vertex_iterator vit = t.vertices().begin(); vit != t.vertices().end(); ++vit) { ++cnt; }
		T::Vertex_range fvr = t.vertices();
		BOOST_FOREACH(T::Vertex &v, fvr) { ++cnt; }
		BOOST_FOREACH(T::Vertex &v, t.vertices()) { ++cnt; }

		// Const iterator over Vertex
		for (T::Vertex_iterator vit = t_const.vertices_begin(); vit != t_const.vertices_end(); ++vit) { ++cnt; }
		for (T::Vertex_iterator vit = t_const.vertices().begin(); vit != t_const.vertices().end(); ++vit) { ++cnt; }
		T::Vertex_range fvr_c = t_const.vertices();
		BOOST_FOREACH(T::Vertex &v, fvr_c) { ++cnt; }
		BOOST_FOREACH(T::Vertex &v, t_const.vertices()) { ++cnt; }

		// Iterator over Finite_vertex_handles
		T::Finite_vertex_handles_range vhs = t.vertex_handles();
		BOOST_FOREACH(T::Vertex_iterator &vh, vhs) { ++cnt; }
		BOOST_FOREACH(T::Vertex_iterator &vh, t.vertex_handles()) { ++cnt; }
		BOOST_FOREACH(T::Vertex_handle vh, t.vertex_handles()) { ++cnt; }

		// Const iterator over Finite_vertex_handles
		T::Vertex_handle_range vhs_c = t_const.vertex_handles();
		BOOST_FOREACH(T::Vertex_iterator &vh, vhs_c) { ++cnt; }
		BOOST_FOREACH(T::Vertex_iterator &vh, t_const.vertex_handles()) { ++cnt; }
		BOOST_FOREACH(T::Vertex_handle vh, t_const.vertex_handles()) { ++cnt; }

		CGAL_assertion(cnt == nVertices * 14);
	}

	{	// All vertices
		size_t nVertices = t.number_of_vertices();
		size_t cnt = 0;
		// Iterator over All_vertices
		for (T::All_vertices_iterator vit = t.all_vertices_begin(); vit != t.all_vertices_end(); ++vit) { ++cnt; }
		for (T::All_vertices_iterator vit = t.all_vertices().begin(); vit != t.all_vertices().end(); ++vit) { ++cnt; }
		T::All_vertices_range fvr = t.all_vertices();
		BOOST_FOREACH(T::Vertex &v, fvr) { ++cnt; }
		BOOST_FOREACH(T::Vertex &v, t.all_vertices()) { ++cnt; }

		// Const iterator over All_vertices
		for (T::All_vertices_iterator vit = t_const.all_vertices_begin(); vit != t_const.all_vertices_end(); ++vit) { ++cnt; }
		for (T::All_vertices_iterator vit = t_const.all_vertices().begin(); vit != t_const.all_vertices().end(); ++vit) { ++cnt; }
		T::All_vertices_range fvr_c = t_const.all_vertices();
		BOOST_FOREACH(T::Vertex &v, fvr_c) { ++cnt; }
		BOOST_FOREACH(T::Vertex &v, t_const.all_vertices()) { ++cnt; }

		// Iterator over All_vertex_handles
		T::All_vertex_handles_range vhs = t.all_vertex_handles();
		BOOST_FOREACH(T::All_vertices_iterator &vh, vhs) { ++cnt; }
		BOOST_FOREACH(T::All_vertices_iterator &vh, t.all_vertex_handles()) { ++cnt; }
		BOOST_FOREACH(T::Vertex_handle vh, t.all_vertex_handles()) { ++cnt; }

		// Const iterator over All_vertex_handles
		T::All_vertex_handles_range vhs_c = t_const.all_vertex_handles();
		BOOST_FOREACH(T::All_vertices_iterator &vh, vhs_c) { ++cnt; }
		BOOST_FOREACH(T::All_vertices_iterator &vh, t_const.all_vertex_handles()) { ++cnt; }
		BOOST_FOREACH(T::Vertex_handle vh, t_const.all_vertex_handles()) { ++cnt; }

		CGAL_assertion(cnt == nVertices * 14);
	}

	{	// Finite faces
		size_t nFaces = t.number_of_faces();
		size_t cnt = 0;
		// Iterator over finite_faces
		for (T::Finite_faces_iterator fit = t.finite_faces_begin(); fit != t.finite_faces_end(); ++fit) { ++cnt; }
		for (T::Finite_faces_iterator fit = t.finite_faces().begin(); fit != t.finite_faces().end(); ++fit) { ++cnt; }
		T::Finite_faces_range fvr = t.finite_faces();
		BOOST_FOREACH(T::Face &f, fvr) { ++cnt; }
		BOOST_FOREACH(T::Face &f, t.finite_faces()) { ++cnt; }

		// Const iterator over finite_faces
		for (T::Finite_faces_iterator fit = t_const.finite_faces_begin(); fit != t_const.finite_faces_end(); ++fit) { ++cnt; }
		for (T::Finite_faces_iterator fit = t_const.finite_faces().begin(); fit != t_const.finite_faces().end(); ++fit) { ++cnt; }
		T::Finite_faces_range fvr_c = t_const.finite_faces();
		BOOST_FOREACH(T::Face &f, fvr_c) { ++cnt; }
		BOOST_FOREACH(T::Face &f, t_const.finite_faces()) { ++cnt; }

		// Iterator over finite_face_handles
		T::Finite_face_handles_range fhs = t.finite_face_handles();
		BOOST_FOREACH(T::Finite_faces_iterator &fh, fhs) { ++cnt; }
		BOOST_FOREACH(T::Finite_faces_iterator &fh, t.finite_face_handles()) { ++cnt; }
		BOOST_FOREACH(T::Face_handle fh, t.finite_face_handles()) { ++cnt; }

		// Const iterator over finite_face_handles
		T::Finite_face_handles_range fhs_c = t_const.finite_face_handles();
		BOOST_FOREACH(T::Finite_faces_iterator &fh, fhs_c) { ++cnt; }
		BOOST_FOREACH(T::Finite_faces_iterator &fh, t_const.finite_face_handles()) { ++cnt; }
		BOOST_FOREACH(T::Face_handle fh, t_const.finite_face_handles()) { ++cnt; }
		
		CGAL_assertion(cnt == nFaces * 14);
	}

	{	// (Finite) faces
		size_t nFaces = t.number_of_faces();
		size_t cnt = 0;
		// Iterator over Finite_faces
		for (T::Face_iterator vit = t.faces_begin(); vit != t.faces_end(); ++vit) { ++cnt; }
		for (T::Face_iterator vit = t.faces().begin(); vit != t.faces().end(); ++vit) { ++cnt; }
		T::Face_range fvr = t.faces();
		BOOST_FOREACH(T::Face &f, fvr) { ++cnt; }
		BOOST_FOREACH(T::Face &f, t.faces()) { ++cnt; }

		// Const iterator over Face
		for (T::Face_iterator vit = t_const.faces_begin(); vit != t_const.faces_end(); ++vit) { ++cnt; }
		for (T::Face_iterator vit = t_const.faces().begin(); vit != t_const.faces().end(); ++vit) { ++cnt; }
		T::Face_range fvr_c = t_const.faces();
		BOOST_FOREACH(T::Face &f, fvr_c) { ++cnt; }
		BOOST_FOREACH(T::Face &f, t_const.faces()) { ++cnt; }

		// Iterator over Finite_face_handles
		T::Finite_face_handles_range fhs = t.face_handles();
		BOOST_FOREACH(T::Face_iterator &fh, fhs) { ++cnt; }
		BOOST_FOREACH(T::Face_iterator &fh, t.face_handles()) { ++cnt; }
		BOOST_FOREACH(T::Face_handle fh, t.face_handles()) { ++cnt; }

		// Const iterator over Finite_face_handles
		T::Face_handle_range fhs_c = t_const.face_handles();
		BOOST_FOREACH(T::Face_iterator &fh, fhs_c) { ++cnt; }
		BOOST_FOREACH(T::Face_iterator &fh, t_const.face_handles()) { ++cnt; }
		BOOST_FOREACH(T::Face_handle fh, t_const.face_handles()) { ++cnt; }

		CGAL_assertion(cnt == nFaces * 14);
	}

	{	// All faces
		size_t nFaces = t.number_of_faces();
		size_t cnt = 0;
		// Iterator over All_faces
		for (T::All_faces_iterator vit = t.all_faces_begin(); vit != t.all_faces_end(); ++vit) { ++cnt; }
		for (T::All_faces_iterator vit = t.all_faces().begin(); vit != t.all_faces().end(); ++vit) { ++cnt; }
		T::All_faces_range fvr = t.all_faces();
		BOOST_FOREACH(T::Face &f, fvr) { ++cnt; }
		BOOST_FOREACH(T::Face &f, t.all_faces()) { ++cnt; }

		// Const iterator over All_faces
		for (T::All_faces_iterator vit = t_const.all_faces_begin(); vit != t_const.all_faces_end(); ++vit) { ++cnt; }
		for (T::All_faces_iterator vit = t_const.all_faces().begin(); vit != t_const.all_faces().end(); ++vit) { ++cnt; }
		T::All_faces_range fvr_c = t_const.all_faces();
		BOOST_FOREACH(T::Face &f, fvr_c) { ++cnt; }
		BOOST_FOREACH(T::Face &f, t_const.all_faces()) { ++cnt; }

		// Iterator over All_face_handles
		T::All_face_handles_range fhs = t.all_face_handles();
		BOOST_FOREACH(T::All_faces_iterator &fh, fhs) { ++cnt; }
		BOOST_FOREACH(T::All_faces_iterator &fh, t.all_face_handles()) { ++cnt; }
		BOOST_FOREACH(T::Face_handle fh, t.all_face_handles()) { ++cnt; }

		// Const iterator over All_face_handles
		T::All_face_handles_range fhs_c = t_const.all_face_handles();
		BOOST_FOREACH(T::All_faces_iterator &fh, fhs_c) { ++cnt; }
		BOOST_FOREACH(T::All_faces_iterator &fh, t_const.all_face_handles()) { ++cnt; }
		BOOST_FOREACH(T::Face_handle fh, t_const.all_face_handles()) { ++cnt; }

		CGAL_assertion(cnt == nFaces * 14);
	}


	{	// Finite edges
		size_t cnt = 0;

		// Iterator range over the finite edges in the triangulation
		BOOST_FOREACH(T::Edge &e, t.finite_edges()) { ++cnt; }
		BOOST_FOREACH(T::Edge e, t.finite_edges()) { ++cnt; }

		// Iterator range over the finite edges in the const triangulation
		BOOST_FOREACH(T::Edge &e, t_const.finite_edges()) { ++cnt; }
		BOOST_FOREACH(T::Edge e, t_const.finite_edges()) { ++cnt; }

		CGAL_assertion(cnt == nEdges * 4);
	}

	{	// (Finite) edges
		size_t cnt = 0;

		// Iterator range over the finite edges in the triangulation
		BOOST_FOREACH(T::Edge &e, t.edges()) { ++cnt; }
		BOOST_FOREACH(T::Edge e, t.edges()) { ++cnt; }

		// Iterator range over the finite edges in the const triangulation
		BOOST_FOREACH(T::Edge &e, t_const.edges()) { ++cnt; }
		BOOST_FOREACH(T::Edge e, t_const.edges()) { ++cnt; }
	}

	{	// All edges
		size_t cnt = 0;

		// Iterator range over the finite edges in the triangulation
		BOOST_FOREACH(T::Edge &e, t.all_edges()) { ++cnt; }
		BOOST_FOREACH(T::Edge e, t.all_edges()) { ++cnt; }

		// Iterator range over the finite edges in the const triangulation
		BOOST_FOREACH(T::Edge &e, t_const.all_edges()) { ++cnt; }
		BOOST_FOREACH(T::Edge e, t_const.all_edges()) { ++cnt; }
	}

	{	// Points
		size_t nVertices = t.number_of_vertices();
		size_t cnt = 0;

		// Iterator range over the points in the triangulation
		BOOST_FOREACH(T::Point &p, t.points()) { ++cnt; }
		BOOST_FOREACH(T::Point p, t.points()) { ++cnt; }

		// Iterator range over the points in the const triangulation
		BOOST_FOREACH(T::Point &p, t_const.points()) { ++cnt; }
		BOOST_FOREACH(T::Point p, t_const.points()) { ++cnt; }

		CGAL_assertion(cnt == nVertices * 4);
	}

	return 0;
}
