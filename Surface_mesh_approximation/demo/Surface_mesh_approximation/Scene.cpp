#include "Scene.h"

#include <iostream>
#include <fstream>

#include <QString>
#include <QTextStream>
#include <QFileInfo>
#include <QInputDialog>

#include <CGAL/Timer.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/centroid.h>

#include <CGAL/internal/Surface_mesh_approximation/VSA_segmentation.h>

#include "ColorCheatSheet.h"

Scene::Scene() :
  m_fidx_pmap(m_fidx_map)
{
  m_pPolyhedron = NULL;

  // view options
  m_view_wireframe = false;
  m_view_seg_boundary = false;
  m_view_anchors = false;

  m_px_num = 0;
}

Scene::~Scene()
{
  delete m_pPolyhedron;
}

void Scene::update_bbox()
{
  std::cout << "Compute bbox...";
  m_bbox = Bbox();

  if(m_pPolyhedron == NULL) {
    std::cout << "failed (no polyhedron)." << std::endl;
    return;
  }

  if(m_pPolyhedron->empty()) {
    std::cout << "failed (empty polyhedron)." << std::endl;
    return;
  }

  Point_iterator it = m_pPolyhedron->points_begin();
  m_bbox = (*it).bbox();
  for(; it != m_pPolyhedron->points_end();it++)
    m_bbox = m_bbox + (*it).bbox();
  std::cout << "done (" << m_pPolyhedron->size_of_facets()
    << " facets)" << std::endl;
}

int Scene::open(QString filename)
{
  QTextStream cerr(stderr);
  cerr << QString("Opening file \"%1\"\n").arg(filename);
  QApplication::setOverrideCursor(QCursor(::Qt::WaitCursor));

  QFileInfo fileinfo(filename);
  std::ifstream in(filename.toUtf8());

  if(!in || !fileinfo.isFile() || ! fileinfo.isReadable()) {
    std::cerr << "unable to open file" << std::endl;
    QApplication::restoreOverrideCursor();
    return -1;
  }

  if(m_pPolyhedron != NULL)
    delete m_pPolyhedron;

  // allocate new polyhedron
  m_pPolyhedron = new Polyhedron;
  in >> *m_pPolyhedron;
  if(!in) {
    std::cerr << "invalid OFF file" << std::endl;
    QApplication::restoreOverrideCursor();

    delete m_pPolyhedron;
    m_pPolyhedron = NULL;

    return -1;
  }

  QApplication::restoreOverrideCursor();
  return 0;
}

void Scene::VSA_segmentation(const std::size_t num_proxies, const std::size_t num_iterations)
{
  if(!m_pPolyhedron)
    return;

  std::cout << "VSA..." << std::endl;

  m_fidx_map.clear();
  for(Facet_const_iterator fitr = m_pPolyhedron->facets_begin();
    fitr != m_pPolyhedron->facets_end();
    ++fitr) {
    m_fidx_map.insert(
      std::pair<Polyhedron::Facet_const_handle, std::size_t>(fitr, 0));
  }

  PointPropertyMap ppmap = get(boost::vertex_point, const_cast<Polyhedron &>(*m_pPolyhedron));

  VSA vsa_seg(*m_pPolyhedron, ppmap, Kernel());
  vsa_seg.partition(num_proxies, num_iterations, m_fidx_pmap);

  std::cerr << "extract mesh" << std::endl;
  vsa_seg.extract_mesh(m_fidx_pmap);
  m_anchors = vsa_seg.collect_anchors();
  m_bdrs = vsa_seg.collect_borders(m_fidx_pmap);
  std::cerr << "#anchors " << m_anchors.size() << std::endl;
  std::cerr << "#borders " << m_bdrs.size() << std::endl;

  m_px_num = num_proxies;
  m_view_seg_boundary = true;

  std::cout << "done" << std::endl;
}

void Scene::VSA_incremental(const std::size_t num_proxies, const std::size_t num_iterations)
{
  if(!m_pPolyhedron)
    return;

  std::cout << "Incremental VSA...";

  m_fidx_map.clear();
  for(Facet_const_iterator fitr = m_pPolyhedron->facets_begin();
    fitr != m_pPolyhedron->facets_end();
    ++fitr) {
    m_fidx_map.insert(
      std::pair<Polyhedron::Facet_const_handle, std::size_t>(fitr, 0));
  }

  typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type PointPropertyMap;
  PointPropertyMap ppmap = get(boost::vertex_point, const_cast<Polyhedron &>(*m_pPolyhedron));

  CGAL::internal::VSA_segmentation<Polyhedron, Kernel, PointPropertyMap> vsa_seg(*m_pPolyhedron, ppmap, Kernel());
  vsa_seg.partition_incre(num_proxies, num_iterations, m_fidx_pmap);

  m_px_num = num_proxies;
  m_view_seg_boundary = true;

  std::cout << "done" << std::endl;
}

void Scene::draw()
{
  if(m_view_wireframe || m_view_seg_boundary) {
    ::glEnable(GL_POLYGON_OFFSET_FILL);
    ::glPolygonOffset(3.0f, 1.0f);
  }
  ::glEnable(GL_LIGHTING);
  render_polyhedron();

  if(m_view_wireframe)
    render_wireframe();
  
  if(m_view_seg_boundary)
    render_segment_boundary();

  if (m_view_anchors) {
    render_anchors();
    render_borders();
  }
}

void Scene::render_polyhedron()
{
  if(!m_pPolyhedron)
    return;

  ::glColor3ub(200, 200, 200);
  ::glBegin(GL_TRIANGLES);
  std::size_t fidx = 0;
  for(Facet_const_iterator fitr = m_pPolyhedron->facets_begin();
    fitr != m_pPolyhedron->facets_end();
    ++fitr) {
    Halfedge_around_facet_const_circulator he = fitr->facet_begin();
    const Point &a = he->vertex()->point();
    const Point &b = he->next()->vertex()->point();
    const Point &c = he->prev()->vertex()->point();

    //Vector norm = CGAL::normal(a, b, c);
    Vector norm = CGAL::unit_normal(a, b, c);
    ::glNormal3d(norm.x(), norm.y(), norm.z());

    if(m_px_num) {
      std::size_t cidx = std::floor(static_cast<double>(m_fidx_pmap[fitr]) / static_cast<double>(m_px_num) * 256.0);
      ::glColor3ub(ColorCheatSheet::r(cidx), ColorCheatSheet::g(cidx), ColorCheatSheet::b(cidx));
    }

    ::glVertex3d(a.x(), a.y(), a.z());
    ::glVertex3d(b.x(), b.y(), b.z());
    ::glVertex3d(c.x(), c.y(), c.z());
  }
  ::glEnd();
}

void Scene::render_wireframe()
{
  if(!m_pPolyhedron)
    return;
  
  // draw black edges
  ::glDisable(GL_LIGHTING);
  ::glColor3ub(0, 0, 0);
  ::glLineWidth(1.0f);
  ::glBegin(GL_LINES);
  for(Edge_iterator he = m_pPolyhedron->edges_begin();
    he != m_pPolyhedron->edges_end();
    he++) {
    const Point& a = he->vertex()->point();
    const Point& b = he->opposite()->vertex()->point();
    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(b.x(),b.y(),b.z());
  }
  ::glEnd();
}

void Scene::render_segment_boundary()
{
  if(!m_pPolyhedron || !m_px_num)
    return;

  ::glDisable(GL_LIGHTING);
  ::glColor3ub(0, 0, 0);
  ::glLineWidth(1.0);
  ::glBegin(GL_LINES);
  for(Edge_const_iterator eitr = m_pPolyhedron->edges_begin();
    eitr != m_pPolyhedron->edges_end();
    ++eitr) {
    std::size_t segid0 = std::numeric_limits<std::size_t>::max();
    if(!eitr->is_border())
      segid0 = m_fidx_pmap[eitr->facet()];
    std::size_t segid1 = std::numeric_limits<std::size_t>::max();
    if(!eitr->opposite()->is_border())
      segid1 = m_fidx_pmap[eitr->opposite()->facet()];

    if(segid0 != segid1) {
      const Point &p0 = eitr->vertex()->point();
      const Point &p1 = eitr->opposite()->vertex()->point();
      ::glVertex3d(p0.x(), p0.y(), p0.z());
      ::glVertex3d(p1.x(), p1.y(), p1.z());
    }
  }
  ::glEnd();
}

void Scene::render_anchors()
{
  ::glDisable(GL_LIGHTING);
  ::glColor3ub(0, 0, 0);
  ::glPointSize(5.0f);
  ::glBegin(GL_POINTS);
  for (std::vector<Anchor>::iterator vitr = m_anchors.begin(); vitr != m_anchors.end(); ++vitr) {
    const Point &pt = vitr->pos;
    ::glVertex3d(pt.x(), pt.y(), pt.z());
  }
  ::glEnd();

  ::glColor3ub(255, 255, 255);
  ::glPointSize(5.0f);
  ::glBegin(GL_POINTS);
  for (std::vector<Anchor>::iterator vitr = m_anchors.begin(); vitr != m_anchors.end(); ++vitr) {
    const Point &pt = vitr->vtx->point();
    ::glVertex3d(pt.x(), pt.y(), pt.z());
  }
  ::glEnd();

  ::glLineWidth(1.0f);
  ::glColor3ub(0, 0, 255);
  ::glBegin(GL_LINES);
  for (std::vector<Anchor>::iterator vitr = m_anchors.begin(); vitr != m_anchors.end(); ++vitr) {
    const Point &ps = vitr->vtx->point();
    ::glVertex3d(ps.x(), ps.y(), ps.z());
    const Point &pt = vitr->pos;
    ::glVertex3d(pt.x(), pt.y(), pt.z());
  }
  ::glEnd();
}

void Scene::render_borders()
{
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1.0f);
  ::glColor3ub(255, 0, 0);
  for (std::vector<std::vector<std::size_t> >::iterator bitr = m_bdrs.begin(); bitr != m_bdrs.end(); ++bitr) {
    ::glBegin(GL_LINE_LOOP);
    for (std::vector<std::size_t>::iterator aitr = bitr->begin(); aitr != bitr->end(); ++aitr) {
      const Point &pt = m_anchors[*aitr].pos;
      ::glVertex3d(pt.x(), pt.y(), pt.z());
    }
    ::glEnd();
  }
}
