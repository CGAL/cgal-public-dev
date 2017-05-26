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
  void VSA_segmentation(const std::size_t num_proxies, const std::size_t num_iterations);

  // toggle view options
  void toggle_view_wireframe() {
    m_view_wireframe = !m_view_wireframe;
  }

  void draw();

private:
  Vector normalize(const Vector& v) {
    return v / std::sqrt(v * v);
  }

  // rendering
  void render_polyhedron_fill();
  void render_polyhedron_wireframe();

private:
  // member data
  Bbox m_bbox;
  Polyhedron *m_pPolyhedron;

  std::vector<std::size_t> m_px_id;
  std::size_t m_px_num;

  // view options
  bool m_view_wireframe;
}; // end class Scene

#endif // SCENE_H
