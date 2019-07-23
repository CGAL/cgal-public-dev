// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_ARRANGEMENT_2_H
#define CGAL_DRAW_ARRANGEMENT_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

//#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Random.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

namespace CGAL
{
  
// Viewer class for Arrangement_2
template <class Arrangement_2>
class SimpleArrangementViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt Base;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle Face_const_handle;
  typedef typename Arrangement_2::Edge_const_iterator Edge_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator Ccb_halfedge_const_circulator;

  typedef typename Arrangement_2::Geometry_traits_2 Geometry_traits_2;
  typedef typename Arrangement_2::Topology_traits Topology_traits;

  typedef typename Geometry_traits_2::Point_2 Point;
  typedef typename Geometry_traits_2::Vector_2 Vector;

  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2> Traits_adapter_2;
  
public:
  /// Construct the viewer.
  /// @param a_arr the arrangement to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  SimpleArrangementViewerQt(QWidget* parent,
                            const Arrangement_2& a_arr,
                            const char* title="Basic Arrangement Viewer",
                            bool anofaces=false,
                            bool aisolatedvertices=true) :
    // First draw: vertices; edges, faces; multi-color; no inverse normal
    Base(parent, title, true, true, true, false, false), 
    arr(a_arr),
    m_nofaces(anofaces),
    m_isolated_vertices(aisolatedvertices)
  {
    geom_traits = static_cast<const Traits_adapter_2*>(arr.geometry_traits());
    top_traits = arr.topology_traits();

    setKeyDescription(::Qt::Key_I, "Toggle view isolated vertices");

    compute_elements();
  }

protected:
  void print_ccb (typename Arrangement_2::Ccb_halfedge_const_circulator circ)
  {
    typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
    //std::cout << "(" << curr->source()->point() << ")";
    do
    {
      Halfedge_const_handle he = curr;
      /*    std::cout << " [" << he->curve() << "] "
            << "(" << he->target()->point() << ")";*/
    }
    while (++curr != circ);
  }
  void compute_face(Face_const_handle fh)
  {
    if (fh->is_unbounded())
    { return; }
    
    CGAL::Random random((unsigned long)(&*fh));
    CGAL::Color c=get_random_color(random);
    
    face_begin(c);

    Ccb_halfedge_const_circulator circ = fh->outer_ccb();
    Ccb_halfedge_const_circulator curr = circ;
    do {
      add_point_in_face(curr->source()->point());
    } while(++curr != circ);

    face_end();

    print_ccb (fh->outer_ccb());
    typename Arrangement_2::Hole_const_iterator hi;
    for (hi=fh->holes_begin(); hi!=fh->holes_end(); ++hi)
    { print_ccb (*hi); }
  }

  void compute_edge(Edge_const_iterator ei)
  {
    add_segment(ei->source()->point(), ei->target()->point());
    return;
  }

  void compute_elements()
  {
    clear();

    // Draw the arrangement vertices.
    typename Arrangement_2::Vertex_const_iterator vit;
    for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
    {
      if (vit->is_isolated())
      { // Isolated vertices are shown in black color
        if (m_isolated_vertices) {
          add_point(vit->point(), CGAL::Color(CGAL::black()));
        }
        continue;
      }
      add_point(vit->point());
    }

    // Draw the arrangement edges.
    typename Arrangement_2::Edge_const_iterator eit;
    for (eit=arr.edges_begin(); eit!=arr.edges_end(); ++eit)
    {
      compute_edge(eit);
      //std::cout << "[" << eit->curve() << "]" << std::endl;
    }

    // Draw the arrangement faces.
    typename Arrangement_2::Face_const_iterator fit;
    for (fit=arr.faces_begin(); fit!=arr.faces_end(); ++fit)
    {
      compute_face(fit);
    }
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    // Test key pressed:
    //    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    //    if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::NoButton)) { ... }
    
    // Call: * compute_elements() if the model changed, followed by
    //       * redraw() if some viewing parameters changed that implies some
    //                  modifications of the buffers
    //                  (eg. type of normal, color/mono)
    //       * update() just to update the drawing

    // Call the base method to process others/classicals key
    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    if ( (e->key() == ::Qt::Key_I)  && (modifiers == ::Qt::NoButton))
    {
      m_isolated_vertices = !m_isolated_vertices;
      displayMessage(QString("Isolated Vertices=%1.").arg(m_isolated_vertices?"true":"false"));
      compute_elements();
      redraw();
    } else {
      Base::keyPressEvent(e);
    }
  }

protected:
  const Arrangement_2& arr;
  const Traits_adapter_2* geom_traits;
  const Topology_traits* top_traits;
  bool m_nofaces;
  bool m_isolated_vertices;
};
  
template<class GeomTraits_, class TopTraits_>
void draw(const Arrangement_2<GeomTraits_, TopTraits_>& a_arr,
          const char* title="Basic Arrangement Viewer",
          bool nofill=false)
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite)
  {
    int argc=1;
    const char* argv[2]={"arr_viewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    SimpleArrangementViewerQt<Arrangement_2<GeomTraits_, TopTraits_> >
      mainwindow(app.activeWindow(), a_arr, title, nofill);
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_ARRANGEMENT_2_H
