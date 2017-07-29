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
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Edge_iterator Edge_iterator;
typedef Polyhedron::Facet_handle Facet_handle;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
typedef CGAL::Bbox_3 Bbox_3;

typedef boost::associative_property_map<std::map<Facet_handle, Vector_3> > FacetNormalMap;
typedef boost::associative_property_map<std::map<Facet_handle, FT> > FacetAreaMap;
typedef boost::associative_property_map<std::map<Facet_handle, Point_3> > FacetCenterMap;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type VertexPointMap;

class Scene
{
public:
  Scene();
  ~Scene();

  void update_bbox();
  Bbox_3 bbox() { return m_bbox; }

  // file menu
  int open(QString filename);
  void save_approximation(const std::string &filename);

  // algorithms
  void l21_approximation(const int &init, const std::size_t num_proxies, const std::size_t num_iterations);
  void compact_approximation(const int &init, const std::size_t num_proxies, const std::size_t num_iterations);
  void l2_approximation(const int &init, const std::size_t num_proxies, const std::size_t num_iterations);

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
  Vector_3 normalize(const Vector_3& v) {
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
  Bbox_3 m_bbox;
  Polyhedron *m_pPolyhedron;

  // property-map for segment-idx
  std::map<Facet_handle, std::size_t> m_fidx_map;
  boost::associative_property_map<std::map<Facet_handle, std::size_t> > m_fidx_pmap;

  // facet property maps
  std::map<Facet_handle, Vector_3> m_facet_normals;
  std::map<Facet_handle, Point_3> m_facet_centers;
  std::map<Facet_handle, FT> m_facet_areas;
  FacetNormalMap m_normal_pmap;
  FacetCenterMap m_center_pmap;
  FacetAreaMap m_area_pmap;
  VertexPointMap m_point_pmap;

  std::vector<Point_3> m_anchor_pos;
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
