#ifndef SCENE_H
#define SCENE_H

#include <QtOpenGL/qgl.h>
#include <iostream>
#include <cmath>

#include "types.h"

class Scene
{
public:
  // types
  typedef CGAL::Bbox_3 Bbox;

public:
  Scene();
  ~Scene();

  void update_bbox();
  Bbox bbox() { return m_bbox; }

  // file menu
  int open(QString filename);

  // algorithms
  void VSA_segmentation();

  // toggle view options
  void toggle_view_poyhedron() {
    m_view_polyhedron = !m_view_polyhedron;
  }

  void draw();

private:
  Vector normalize(const Vector& v) {
    return v / std::sqrt(v * v);
  }

  // rendering
  void render_polyhedron();

private:
  // member data
  Bbox m_bbox;
  Polyhedron *m_pPolyhedron;

  // view options
  bool m_view_polyhedron;
}; // end class Scene

#endif // SCENE_H
