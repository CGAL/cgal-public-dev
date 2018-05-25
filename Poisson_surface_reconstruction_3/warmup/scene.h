#ifndef SCENE_H
#define SCENE_H

#include <QString>
#include <list>

#include "types.h"
#include"c2t3.h"
#undef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

class Scene
{
private:

	// input point set
	Bbox m_bbox;

	// Complex 3 in Delaunay triangulation
	STr tr;
	C2T3<Kernel, STr> m_c2t3;


	// rendering options
	bool m_view_mesh;
	bool m_view_edges;
	bool m_view_vertices;
	PointList points;

public: // life cycle

	Scene();
	virtual ~Scene();

	// file menu
	int open(QString filename);
	void read_xyz(QString filename);

	// algorithms menu
	void mesh_torus(const FT angle, const FT sizing, const FT approximation);
	void mesh_sphere(const FT angle, const FT sizing, const FT approximation);
	void mesh_ellipsoid(const FT angle, const FT sizing, const FT approximation);
	void implicit_function();

	// toggle rendering options
	void toggle_view_mesh()  { m_view_mesh = !m_view_mesh; }
	void toggle_view_edges()    { m_view_edges = !m_view_edges; }
	void toggle_view_vertices()  { m_view_vertices = !m_view_vertices; }

	// rendering
	void render();
	void render_mesh();
	void render_edges();
	void render_vertices();

};

#endif // SCENE_H
