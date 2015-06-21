// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Adrien Pilleboue

#ifndef CGAL_LCC_3_VIEWER_QT_H
#define CGAL_LCC_3_VIEWER_QT_H

#include <QApplication>
#include <QKeyEvent>

#include <QGLViewer/qglviewer.h>
#include <GL/gl.h>
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_converter.h>

typedef CGAL::Cartesian<double> Local_kernel;
typedef typename Local_kernel::Point_3  Local_point;
typedef typename Local_kernel::Vector_3 Local_vector;

template<class LCC, int dim=LCC::ambient_dimension>
struct Geom_utils;

template<class LCC>
struct Geom_utils<LCC,3>
{
  Local_point get_point(LCC& lcc, typename LCC::Vertex_attribute_const_handle vh)
  { return converter(lcc.point_of_vertex_attribute(vh)); }

  Local_point get_point(LCC& lcc, typename LCC::Dart_const_handle dh)
  { return converter(lcc.point(dh)); }

  Local_vector get_facet_normal(LCC& lcc, typename LCC::Dart_const_handle dh)
  {
    Local_vector n = converter(CGAL::compute_normal_of_cell_2<LCC>(lcc,dh));
    n = n/(CGAL::sqrt(n*n));
    return n;
  }

  Local_vector get_vertex_normal(LCC& lcc, typename LCC::Dart_const_handle dh)
  {
    Local_vector n = converter(CGAL::compute_normal_of_cell_0<LCC>(lcc,dh));
    n = n/(CGAL::sqrt(n*n));
    return n;
  }
protected:
  CGAL::Cartesian_converter<typename LCC::Traits, Local_kernel> converter;
};

template<class LCC>
struct Geom_utils<LCC,2>
{
  Local_point get_point(LCC& lcc, typename LCC::Vertex_attribute_const_handle vh)
  {
    Local_point p(converter(lcc.point_of_vertex_attribute(vh).x()),0,
                  converter(lcc.point_of_vertex_attribute(vh).y()));
    return p;
  }

  Local_point get_point(LCC& lcc, typename LCC::Dart_const_handle dh)
  { return get_point(lcc, lcc.vertex_attribute(dh)); }

  Local_vector get_facet_normal(LCC&, typename LCC::Dart_const_handle)
  {
    Local_vector n(0,1,0);
    return n;
  }

  Local_vector get_vertex_normal(LCC&, typename LCC::Dart_const_handle)
  {
    Local_vector n(0,1,0);
    return n;
  }
protected:
  CGAL::Cartesian_converter<typename LCC::Traits, Local_kernel> converter;
};

template<class LCC>
CGAL::Bbox_3 bbox(LCC& lcc)
{
  CGAL::Bbox_3 bb;
  Geom_utils<LCC> geomutils;

  typename LCC::Vertex_attribute_range::const_iterator
    it=lcc.vertex_attributes().begin(), itend=lcc.vertex_attributes().end();
  if ( it!=itend )
  {
    bb = geomutils.get_point(lcc, it).bbox();
    for( ++it; it!=itend; ++it)
    {
      bb = bb + geomutils.get_point(lcc, it).bbox();
    }
  }

  return bb;
}

template<class LCC>
class SimpleLCCViewerQt : public QGLViewer
{
  typedef typename LCC::Dart_handle Dart_handle;

public:

  // Constructor/Destructor
  SimpleLCCViewerQt(LCC& alcc) :
    lcc(alcc),
    wireframe(false),
    flatShading(true),
    edges(true),
    vertices(true),
    m_displayListCreated(false)
  {
    setWindowTitle("3D lcc viewer");
    resize(500, 450);

    QGLFormat newFormat = this->format();
    newFormat.setSampleBuffers(true);
    newFormat.setSamples(16);
    this->setFormat(newFormat);
  }

protected :
  void drawAllFaces(bool flat)
  {
    for(typename LCC::template One_dart_per_cell_range<2>::iterator
          dartIter=lcc.template one_dart_per_cell<2>().begin();
        dartIter.cont(); ++dartIter)
    {
      Dart_handle& dart = dartIter;
      ::glBegin(GL_POLYGON);
      ::glColor3f(1.0f, .7f, .7f);

      if(flat)
      {
        Local_vector normal = geomutils.get_facet_normal(lcc,dart);
        ::glNormal3d(normal.x(), normal.y(), normal.z());
      }

      for (typename LCC::template Dart_of_orbit_range<1>::const_iterator
             orbitIter = lcc.template darts_of_orbit<1>(dart).begin();
           orbitIter.cont(); ++orbitIter)
      {
        if(!flat)
        {
          // If Gouraud shading: 1 normal per vertex
          Local_vector n = geomutils.get_vertex_normal(lcc,orbitIter);
          ::glNormal3d(n.x(),n.y(),n.z());
        }

        Local_point p = geomutils.get_point(lcc, orbitIter);
        ::glVertex3d(p.x(),p.y(),p.z());
      }

      ::glEnd();
    }
  }

  void drawAllEdges()
  {
    //    ::glDepthRange(0.0, 1.0-0.005);
    ::glBegin(GL_LINES);
    ::glColor3f(0.0f, 0.0f, 0.0f);

    for(typename LCC::template One_dart_per_cell_range<1>::iterator
          dartIter=lcc.template one_dart_per_cell<1>().begin();
        dartIter.cont(); ++dartIter)
    {
      Dart_handle& dart = dartIter;

      Local_point p =  geomutils.get_point(lcc, dartIter);
      Dart_handle d2 = lcc.other_extremity(dartIter);
      if ( d2!=NULL )
      {
        Local_point p2 = geomutils.get_point(lcc, d2);
        glVertex3f( p.x(),p.y(),p.z());
        glVertex3f( p2.x(),p2.y(),p2.z());
      }
    }

    ::glEnd();
    //    ::glDepthRange(0.005, 1.0);
  }

  void drawAllVertices()
  {
    //    ::glDepthRange(0.0, 1.0-0.005);
    ::glPointSize(7.0);
    ::glBegin(GL_POINTS);
    ::glColor3f(0.2f, 0.2f, 0.7f);

    for (typename LCC::Vertex_attribute_const_range::iterator
           v=lcc.vertex_attributes().begin(),
           vend=lcc.vertex_attributes().end();
         v!=vend; ++v)
    {
      Local_point p = geomutils.get_point(lcc, v);
      glVertex3f(p.x(), p.y(), p.z());
    }

    ::glEnd();
    //    ::glDepthRange(0.005, 1.0);
  }

  void initDraw()
  {
    //Compile drawFacet
    std::cout << "Compile Display Lists : Faces," << std::flush;
    m_dlFaces = ::glGenLists(1);
    ::glNewList(m_dlFaces, GL_COMPILE);
    drawAllFaces(false);
    ::glEndList();

    //Compile drawFacet with flat shading
    std::cout << "Faces (flat shading), " << std::flush;
    m_dlFacesFlat = ::glGenLists(1);
    ::glNewList(m_dlFacesFlat, GL_COMPILE);
    drawAllFaces(true);
    ::glEndList();

    //Compile drawEdge
    std::cout << "edges, " << std::flush;
    m_dlEdges = ::glGenLists(1);
    ::glNewList(m_dlEdges, GL_COMPILE);
    drawAllEdges();
    ::glEndList();

    //Compile drawvertices
    std::cout << "vertices" << std::flush;
    m_dlVertices = ::glGenLists(1);
    ::glNewList(m_dlVertices, GL_COMPILE);
    drawAllVertices();
    ::glEndList();

     std::cout << ". DONE." << std::endl;
     m_displayListCreated = true;
  }

  virtual void draw()
  {
    if(!m_displayListCreated) initDraw();

    if ( !wireframe )
    {
      if(flatShading) ::glCallList(m_dlFacesFlat);
      else ::glCallList(m_dlFaces);
    }

    if(edges) ::glCallList(m_dlEdges);

    if(vertices) ::glCallList(m_dlVertices);
  }

  virtual void init()
  {
    // Restore previous viewer state.
    restoreStateFromFile();

    // Define 'Control+Q' as the new exit shortcut (default was 'Escape')
    setShortcut(EXIT_VIEWER, Qt::CTRL+Qt::Key_Q);

    // Add custom key description (see keyPressEvent).
    setKeyDescription(Qt::Key_W, "Toggles wire frame display");
    setKeyDescription(Qt::Key_F, "Toggles flat shading display");
    setKeyDescription(Qt::Key_E, "Toggles edges display");
    setKeyDescription(Qt::Key_V, "Toggles vertices display");

    // Light default parameters
    ::glLineWidth(2.4f);
    ::glPointSize(7.f);
    //    ::glEnable(GL_POLYGON_OFFSET_FILL);
    //    ::glPolygonOffset(1.f,1.f);
    ::glClearColor(1.0f,1.0f,1.0f,0.0f);
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

    ::glEnable(GL_LIGHTING);

    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    // ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

 //   ::glDepthRange(0.005, 1.0);
    if (flatShading)
    {
      ::glShadeModel(GL_FLAT);
      ::glDisable(GL_BLEND);
      ::glDisable(GL_LINE_SMOOTH);
      ::glDisable(GL_POLYGON_SMOOTH_HINT);
      ::glBlendFunc(GL_ONE, GL_ZERO);
      ::glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
    }
    else
    {
      ::glShadeModel(GL_SMOOTH);
      ::glEnable(GL_BLEND);
      ::glEnable(GL_LINE_SMOOTH);
      ::glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
      ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }

    CGAL::Bbox_3 bb = bbox(lcc);

    this->camera()->setSceneBoundingBox(qglviewer::Vec(bb.xmin(),
                                                       bb.ymin(),
                                                       bb.zmin()),
                                        qglviewer::Vec(bb.xmax(),
                                                       bb.ymax(),
                                                       bb.zmax()));

    this->showEntireScene();
    initDraw();
  }

  void keyPressEvent(QKeyEvent *e)
  {
    const Qt::KeyboardModifiers modifiers = e->modifiers();

    bool handled = false;
    if ((e->key()==Qt::Key_W) && (modifiers==Qt::NoButton))
    {
      wireframe = !wireframe;
      if (wireframe)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      handled = true;
      updateGL();
    }
    else if ((e->key()==Qt::Key_F) && (modifiers==Qt::NoButton))
    {
      flatShading = !flatShading;
      if (flatShading)
      {
        ::glShadeModel(GL_FLAT);
        ::glDisable(GL_BLEND);
        ::glDisable(GL_LINE_SMOOTH);
        ::glDisable(GL_POLYGON_SMOOTH_HINT);
        ::glBlendFunc(GL_ONE, GL_ZERO);
        ::glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
      }
      else
      {
        ::glShadeModel(GL_SMOOTH);
        ::glEnable(GL_BLEND);
        ::glEnable(GL_LINE_SMOOTH);
        ::glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      }
      handled = true;
      updateGL();
    }
    else if ((e->key()==Qt::Key_E) && (modifiers==Qt::NoButton))
    {
      edges = !edges;
      handled = true;
      updateGL();
    }
    else if ((e->key()==Qt::Key_V) && (modifiers==Qt::NoButton))
    {
      vertices = !vertices;
      handled = true;
      updateGL();
    }

    if (!handled)
      QGLViewer::keyPressEvent(e);
  }


  virtual QString helpString() const
  {
    QString text("<h2>L C C   V i e w e r</h2>");
    text += "Use the mouse to move the camera around the object. ";
    text += "You can respectively revolve around, zoom and translate with "
      "the three mouse buttons. ";
    text += "Left and middle buttons pressed together rotate around the "
      "camera view direction axis<br><br>";
    text += "Pressing <b>Alt</b> and one of the function keys "
      "(<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
    text += "Simply press the function key again to restore it. "
      "Several keyFrames define a ";
    text += "camera path. Paths are saved when you quit the application "
      "and restored at next start.<br><br>";
    text += "Press <b>F</b> to display the frame rate, <b>A</b> for the "
      "world axis, ";
    text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> "
      "to save a snapshot. ";
    text += "See the <b>Keyboard</b> tab in this window for a complete "
      "shortcut list.<br><br>";
    text += "Double clicks automates single click actions: A left button "
      "double click aligns the closer axis with the camera (if close enough). ";
    text += "A middle button double click fits the zoom of the camera and "
      "the right button re-centers the scene.<br><br>";
    text += "A left button double click while holding right button pressed "
      "defines the camera <i>Revolve Around Point</i>. ";
    text += "See the <b>Mouse</b> tab and the documentation web pages for "
      "details.<br><br>";
    text += "Press <b>Escape</b> to exit the viewer.";
    return text;
  }
private:
  LCC& lcc;
  bool wireframe;
  bool flatShading;
  bool edges;
  bool vertices;
  Geom_utils<LCC> geomutils;

  GLuint m_dlFaces;
  GLuint m_dlFacesFlat;
  GLuint m_dlEdges;
  GLuint m_dlVertices;
  bool m_displayListCreated;
};

template<class LCC>
void display_lcc(LCC& alcc)
{
  int argc=1;
  typedef char* s;

  const char* argv[2]={"lccviewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));

  SimpleLCCViewerQt<LCC> mainwindow(alcc);
  mainwindow.show();

  app.exec();
};

#endif // CGAL_LCC_3_VIEWER_QT_H
