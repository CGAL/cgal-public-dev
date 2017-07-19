#ifndef SCENE_H
#define SCENE_H

#include <QtOpenGL/qgl.h>
#include <iostream>
#include <cmath>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/property_map.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet_const_handle Facet_const_handle;
typedef Polyhedron::Halfedge_const_handle Halfedge_const_handle;
typedef Polyhedron::Facet_const_iterator Facet_const_iterator;
typedef Polyhedron::Halfedge_around_facet_const_circulator Halfedge_around_facet_const_circulator;
typedef Polyhedron::Edge_const_iterator Edge_const_iterator;
typedef CGAL::Bbox_3 Bbox;

class Scene
{
public:
  Scene();
  ~Scene();

  void update_bbox();
  Bbox bbox() { return m_bbox; }

  // file menu
  int open(QString filename);
  void save_approximation(const std::string &filename);

  // algorithms
  void variational_shape_approximation(const int &init, const std::size_t num_proxies, const std::size_t num_iterations);
  void compact_approximation(const int &init, const std::size_t num_proxies, const std::size_t num_iterations);

  // toggle view options
  void toggle_view_polyhedron() {
    m_view_polyhedron = !m_view_polyhedron;
  }

  void toggle_view_wireframe() {
    m_view_wireframe = !m_view_wireframe;
  }

  void toggle_view_seg_boundary() {
    m_view_seg_boundary = !m_view_seg_boundary;
  }

  void toggle_view_anchors() {
    m_view_anchors = !m_view_anchors;
  }

  void toggle_view_approximation() {
    m_view_approximation = !m_view_approximation;
  }

  void draw();

private:
  Vector normalize(const Vector& v) {
    return v / std::sqrt(v * v);
  }

  // rendering
  void render_polyhedron();
  void render_wireframe();
  void render_segment_boundary();
  void render_anchors();
  void render_borders();
  void render_approximation();

private:
  // member data
  Bbox m_bbox;
  Polyhedron *m_pPolyhedron;

  typedef std::map<Polyhedron::Facet_const_handle, std::size_t> FacetIdMap;
  typedef boost::associative_property_map<FacetIdMap> FacetIdPmap;
  FacetIdMap m_fidx_map;
  FacetIdPmap m_fidx_pmap; // property-map for segment-idx

  std::vector<Point> m_anchor_pos;
  std::vector<Polyhedron::Vertex_handle> m_anchor_vtx;
  std::vector<std::vector<std::size_t> > m_bdrs; // anchor borders
  std::vector<int> m_tris;

  std::size_t m_px_num;

  // view options
  bool m_view_polyhedron;
  bool m_view_wireframe;
  bool m_view_seg_boundary;
  bool m_view_anchors;
  bool m_view_approximation;
}; // end class Scene

#endif // SCENE_H
