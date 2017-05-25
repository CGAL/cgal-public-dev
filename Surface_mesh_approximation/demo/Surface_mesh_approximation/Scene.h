#ifndef SCENE_H
#define SCENE_H

#include <QtOpenGL/qgl.h>
#include <iostream>
#include <cmath>

#include "types.h"

class Scene
{
public:
  Scene();
  ~Scene();
public:
  // types
  typedef CGAL::Bbox_3 Bbox;

public:
  void update_bbox();
  Bbox bbox() { return m_bbox; }

private:
  // member data
  Bbox m_bbox;
  Line m_line;
  Plane m_plane;
  Point m_centroid;
  Polyhedron *m_pPolyhedron;

  // view options
  bool m_view_polyhedron;

public:
  // file menu
  int open(QString filename);

  Vector normalize(const Vector& v);

  // algorithms
  void refine_loop();
  void fit_triangles();
  void fit_edges();
  void fit_vertices();
  
  // toggle view options
  void toggle_view_poyhedron();

  // rendering
  void draw(); 
  void render_plane();
  void render_line();
  void render_centroid();
  void render_polyhedron();

private:

}; // end class Scene

#endif // SCENE_H
