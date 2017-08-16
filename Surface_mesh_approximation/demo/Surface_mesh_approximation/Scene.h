#ifndef SCENE_H
#define SCENE_H

#include <QtOpenGL/qgl.h>
#include <iostream>
#include <cmath>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/property_map.h>
#include <CGAL/vsa_mesh_approximation_traits.h>
#include <CGAL/VSA_approximation.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef Polyhedron_3::Halfedge_handle Halfedge_handle;
typedef Polyhedron_3::Edge_iterator Edge_iterator;
typedef Polyhedron_3::Facet_handle Facet_handle;
typedef Polyhedron_3::Facet_iterator Facet_iterator;
typedef Polyhedron_3::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
typedef CGAL::Bbox_3 Bbox_3;

typedef boost::associative_property_map<std::map<Facet_handle, Vector_3> > FacetNormalMap;
typedef boost::associative_property_map<std::map<Facet_handle, FT> > FacetAreaMap;
typedef boost::associative_property_map<std::map<Facet_handle, Point_3> > FacetCenterMap;
typedef boost::property_map<Polyhedron_3, boost::vertex_point_t>::type VertexPointMap;


typedef CGAL::PlaneProxy<Polyhedron_3> PlaneProxy;

typedef CGAL::L21Metric<Polyhedron_3, FacetNormalMap, FacetAreaMap> L21Metric;
typedef CGAL::L21ProxyFitting<Polyhedron_3, FacetNormalMap, FacetAreaMap> L21ProxyFitting;
typedef CGAL::VSA_approximation<Polyhedron_3, PlaneProxy, L21Metric, L21ProxyFitting> VSAL21;

typedef CGAL::L2Metric<Polyhedron_3, FacetAreaMap> L2Metric;
typedef CGAL::L2ProxyFitting<Polyhedron_3> L2ProxyFitting;
typedef CGAL::VSA_approximation<Polyhedron_3, PlaneProxy, L2Metric, L2ProxyFitting> VSAL2;

// user defined compact metric
struct PointProxy {
  Facet_handle seed;
  Point_3 center;
};

struct CompactMetric {
  typedef PointProxy Proxy;

  CompactMetric(const FacetCenterMap &_center_pmap)
    : center_pmap(_center_pmap) {}

  FT operator()(const Facet_handle &f, const PointProxy &px) const {
    return FT(std::sqrt(CGAL::to_double(
      CGAL::squared_distance(center_pmap[f], px.center))));
  }

  const FacetCenterMap center_pmap;
};

struct PointProxyFitting {
  typedef PointProxy Proxy;

  PointProxyFitting(const FacetCenterMap &_center_pmap,
    const FacetAreaMap &_area_pmap)
    : center_pmap(_center_pmap),
    area_pmap(_area_pmap) {}

  template<typename FacetIterator>
  PointProxy operator()(const FacetIterator beg, const FacetIterator end) const {
    CGAL_assertion(beg != end);

    // fitting center
    Vector_3 center = CGAL::NULL_VECTOR;
    FT area(0);
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      center = center + (center_pmap[*fitr] - CGAL::ORIGIN) * area_pmap[*fitr];
      area += area_pmap[*fitr];
    }
    center = center / area;

    // construct proxy
    PointProxy px;
    px.center = CGAL::ORIGIN + center;

    return px;
  }

  const FacetCenterMap center_pmap;
  const FacetAreaMap area_pmap;
};
typedef CGAL::VSA_approximation<Polyhedron_3, PointProxy, CompactMetric, PointProxyFitting> VSACompact;

class Scene
{
  enum Metric { L21, L2, Compact };

public:
  Scene() :
    m_pmesh(NULL),
    m_fidx_pmap(m_fidx_map),
    m_normal_pmap(m_facet_normals),
    m_center_pmap(m_facet_centers),
    m_area_pmap(m_facet_areas),
    m_pl21_metric(NULL),
    m_pl21_proxy_fitting(NULL),
    m_pl2_metric(NULL),
    m_pl2_proxy_fitting(NULL),
    m_pcompact_metric(NULL),
    m_pcompact_proxy_fitting(NULL),
    m_px_num(0),
    m_view_polyhedron(false),
    m_view_wireframe(false),
    m_view_seg_boundary(false),
    m_view_anchors(false) {}

  ~Scene() {
    delete m_pmesh;
    if (m_pl21_metric)
      delete m_pl21_metric;
    if (m_pl21_proxy_fitting)
      delete m_pl21_proxy_fitting;
    if (m_pl2_metric)
      delete m_pl2_metric;
    if (m_pl2_proxy_fitting)
      delete m_pl2_proxy_fitting;
    if (m_pcompact_metric)
      delete m_pcompact_metric;
    if (m_pcompact_proxy_fitting)
      delete m_pcompact_proxy_fitting;
  }

  void update_bbox();
  Bbox_3 bbox() { return m_bbox; }

  // file menu
  int open(QString filename);
  void save_approximation(const std::string &filename);

  // algorithms
  void l21_approximation(const int &init, const std::size_t num_proxies, const std::size_t num_iterations);
  void l2_approximation(const int &init, const std::size_t num_proxies, const std::size_t num_iterations);
  void compact_approximation(const int &init, const std::size_t num_proxies, const std::size_t num_iterations);
  void meshing();

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
  Polyhedron_3 *m_pmesh;

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

  // algorithm instance
  Metric m_metric; // current metric
  L21Metric *m_pl21_metric;
  L21ProxyFitting *m_pl21_proxy_fitting;
  VSAL21 m_vsa_l21;

  L2Metric *m_pl2_metric;
  L2ProxyFitting *m_pl2_proxy_fitting;
  VSAL2 m_vsa_l2;

  CompactMetric *m_pcompact_metric;
  PointProxyFitting *m_pcompact_proxy_fitting;
  VSACompact m_vsa_compact;

  std::vector<Point_3> m_anchor_pos;
  std::vector<Polyhedron_3::Vertex_handle> m_anchor_vtx;
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
