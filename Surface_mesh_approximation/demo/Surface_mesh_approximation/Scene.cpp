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

Scene::Scene()
{
  m_pPolyhedron = NULL;

  // view options
  m_view_polyhedron = true;
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

void Scene::VSA_segmentation()
{
  std::cout << "VSA...";

  std::cout << "done" << std::endl;
}

void Scene::draw()
{
  if(m_view_polyhedron) {
    render_polyhedron_fill();
    render_polyhedron_wireframe();
  }
}

void Scene::render_polyhedron_fill()
{
  if(!m_pPolyhedron)
    return;

  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(3.0f, 1.0f);
  ::glColor3ub(200, 200, 200);
  ::glBegin(GL_TRIANGLES);
  for(Facet_iterator fitr = m_pPolyhedron->facets_begin();
    fitr != m_pPolyhedron->facets_end();
    ++fitr) {
    Halfedge_around_facet_circulator he = fitr->facet_begin();
    const Point &a = he->vertex()->point();
    const Point &b = he->next()->vertex()->point();
    const Point &c = he->prev()->vertex()->point();

    //Vector norm = CGAL::normal(a, b, c);
    Vector norm = CGAL::unit_normal(a, b, c);
    ::glNormal3d(norm.x(), norm.y(), norm.z());

    ::glVertex3d(a.x(), a.y(), a.z());
    ::glVertex3d(b.x(), b.y(), b.z());
    ::glVertex3d(c.x(), c.y(), c.z());
  }
  ::glEnd();
}

void Scene::render_polyhedron_wireframe()
{
  // draw black edges
  if(m_pPolyhedron != NULL)
  {
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
}
