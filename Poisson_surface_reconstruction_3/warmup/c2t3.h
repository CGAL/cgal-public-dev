#ifndef _My_C2T3_
#define _My_C2T3_

#include <CGAL/basic.h>

#include <QtOpenGL>
#include <QtOpenGL/qgl.h>
#include"types.h"
#include <CGAL/Surface_mesh_complex_2_in_triangulation_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#undef min
#undef max


template < class Kernel, class Tr >
class C2T3 : public CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr>
{
public:
	typedef C2T3<Kernel, Tr> C2t3;

	typedef typename Kernel::FT            FT;
	typedef typename Kernel::Point_3       Point;
	typedef typename Kernel::Vector_3      Vector;
	typedef typename Kernel::Segment_3     Segment;
	typedef typename Kernel::Line_3        Line;
	typedef typename Kernel::Triangle_3    Triangle;
	typedef typename Kernel::Tetrahedron_3 Tetrahedron;

	typedef typename CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr> Mesh3;
	typedef typename CGAL::Surface_mesh_default_triangulation_3 STriangulation;
	//typedef typename C2t3::Triangulation STriangulation;
	typedef typename STriangulation::Vertex_handle Vertex_handle;
	typedef typename STriangulation::Edge Edge;
	typedef typename STriangulation::Facet Facet;
	typedef typename STriangulation::Cell_handle Cell_handle;
	Polyhedron polyhedron;
	using Mesh3:: Mesh3;

public:

	Vertex_handle get_source_vertex(const Edge& edge) const
	{
		return edge.first->vertex(edge.second);
	}

	Vertex_handle get_target_vertex(const Edge& edge) const
	{
		return edge.first->vertex(edge.third);
	}

	void set_polyhedron(Polyhedron& p){
		polyhedron = p;
	}

	// RENDERING
	void gl_vertex(const Point& p)
	{
		::glVertex3d(p.x(),p.y(),p.z());
	}

	void render_edges(const float line_width,
		const unsigned char red,
		const unsigned char green,
		const unsigned char blue)
	{
		::glLineWidth(line_width);
		::glColor3ub(red,green,blue);
		::glBegin(GL_LINES);

		STriangulation& tr = this->triangulation();

		//typename STriangulation::Finite_edges_iterator e;
		typename C2t3::Edge_iterator e;
		for (e = this->edges_begin(); e != this->edges_end(); e++)
		{
			typename STriangulation::Edge edge = *e;
			gl_vertex(get_source_vertex(edge)->point());
			gl_vertex(get_target_vertex(edge)->point());
		}
		::glEnd();
	}

	int nb_vertices()
	{
		STriangulation& tr = this->triangulation();
		return (int)tr.number_of_vertices();
	}

	void render_vertices(const float point_size,
		const unsigned char red,
		const unsigned char green,
		const unsigned char blue)
	{
		::glPointSize(point_size);
		::glColor3ub(red, green, blue);

		STriangulation& tr = this->triangulation();
		const int nbv = (int)tr.number_of_vertices();

		if (nbv == 0)
			return;
		// std::cout << nbv << " vertices" << std::endl;

		::glBegin(GL_POINTS);
		///typename STriangulation::Finite_vertices_iterator v;
		typename C2t3::Vertex_iterator v;
		for (v = this->vertices_begin(); v != this->vertices_end(); v++)
			gl_vertex(v->point()); //just v->point() would do?
		::glEnd();

	}

	void render_facets(const FT xcut){
		::glEnable(GL_CULL_FACE);
		::glCullFace(GL_BACK);
		::glColor3ub(128, 128, 200);
		::glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		::glEnable(GL_LIGHTING);

		::glBegin(GL_TRIANGLES);
		STr& tr = this->triangulation();
		typename C2t3::Facet_iterator fi;
		int i = 0;
		std::map<Vertex_handle, int> V;
		int inum = 0;

		for(fi = this->facets_begin(); fi != this->facets_end(); fi++){

			i++;

			const Point& a = fi->first->vertex(tr.vertex_triple_index(fi->second, 0))->point();
			const Point& b = fi->first->vertex(tr.vertex_triple_index(fi->second, 1))->point();
			const Point& c = fi->first->vertex(tr.vertex_triple_index(fi->second, 2))->point();


			Point cc = CGAL::centroid(a, b, c);

			gl_shaded_triangle(a, b, c);
		}

		std::cout << "number of facets : " << i << std::endl;
		::glEnd();
		::glDisable(GL_LIGHTING);

	}

	void gl_shaded_triangle(const Point& a, const Point& b, const Point& c)
	{
		// compute normal
		Vector n = CGAL::cross_product(c-a, b-a);
		n = n / std::sqrt(n*n);

		// draw one front facing
		::glNormal3d(-n.x(), -n.y(), -n.z());
		::glVertex3d(a.x(), a.y(), a.z());
		::glVertex3d(b.x(), b.y(), b.z());
		::glVertex3d(c.x(), c.y(), c.z());

		// and the other back facing
		::glNormal3d(n.x(), n.y(), n.z());
		::glVertex3d(a.x(), a.y(), a.z());
		::glVertex3d(c.x(), c.y(), c.z());
		::glVertex3d(b.x(), b.y(), b.z());
	}

};

#endif // _My_C2T3_
