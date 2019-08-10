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

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>

#include <CGAL/Arr_circular_arc_traits_2.h>

namespace CGAL
{

// Traits Adaptor (demo/Utils.h)
// --------------
template < class ArrTraits >
class ArrTraitsAdaptor
{ };

template < class Kernel_ >
class ArrTraitsAdaptor< CGAL::Arr_segment_traits_2< Kernel_ > >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_segment_traits_2< Kernel > ArrTraits;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template < class Kernel_ >
class ArrTraitsAdaptor< CGAL::Arr_linear_traits_2< Kernel_ > >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_linear_traits_2< Kernel > ArrTraits;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template < class SegmentTraits >
class ArrTraitsAdaptor< CGAL::Arr_polyline_traits_2< SegmentTraits > >
{
public:
  typedef CGAL::Arr_polyline_traits_2< SegmentTraits > ArrTraits;
  typedef typename SegmentTraits::Kernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template < class CircularKernel >
class ArrTraitsAdaptor< CGAL::Arr_circular_arc_traits_2< CircularKernel > >
{
public:
  typedef CGAL::Arr_circular_arc_traits_2< CircularKernel > ArrTraits;
  typedef CircularKernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::Root_of_2 CoordinateType;
};

template <class CircularKernel >
class ArrTraitsAdaptor<CGAL::Arr_circle_segment_traits_2<CircularKernel>>
{
  public:
    typedef CGAL::Arr_circle_segment_traits_2<CircularKernel> ArrTraits;
    typedef CircularKernel Kernel;
    typedef typename ArrTraits::Point_2 Point_2;
    typedef typename Kernel::FT CoordinateType;
};

template < class RatKernel, class AlgKernel, class NtTraits >
class ArrTraitsAdaptor< CGAL::Arr_conic_traits_2< RatKernel, AlgKernel,
                                                  NtTraits > >
{
public:
  typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > ArrTraits;
  typedef AlgKernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template < class Coefficient_ >
class ArrTraitsAdaptor< CGAL::Arr_algebraic_segment_traits_2< Coefficient_ > >
{
public:
  typedef Coefficient_ Coefficient;
  typedef typename CGAL::Arr_algebraic_segment_traits_2<Coefficient>
                                                        ArrTraits;
  typedef typename ArrTraits::Point_2                   Point_2; // CKvA_2
  typedef typename ArrTraits::Algebraic_real_1          CoordinateType;
  typedef CGAL::Cartesian< typename ArrTraits::Bound >  Kernel;
  //typedef typename ArrTraits::CKvA_2                  Kernel;
};
// ---------------
  
// Viewer class for Arrangement_2
template <class Arrangement_2>
class SimpleArrangementViewerQtBase : public Basic_viewer_qt
{
public:
  typedef typename Arrangement_2::Geometry_traits_2 Traits;
  typedef typename ArrTraitsAdaptor<Traits>::Kernel Kernel;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Triangle_2 Triangle_2;
  typedef typename Kernel::Iso_rectangle_2 Iso_rectangle_2;
  typedef typename Kernel::Circle_2 Circle_2;

  typedef Basic_viewer_qt Base;
  //  typedef typename Arrangement_2::Halfedge_const_handle
  //  Halfedge_const_handle; typedef typename Arrangement_2::Face_const_handle
  //  Face_const_handle; typedef typename Arrangement_2::Edge_const_iterator
  //  Edge_const_iterator; typedef typename
  //  Arrangement_2::Ccb_halfedge_const_circulator
  //  Ccb_halfedge_const_circulator;

  //  typedef typename Arrangement_2::Geometry_traits_2 Geometry_traits_2;
  //  typedef typename Arrangement_2::Topology_traits Topology_traits;

  //  typedef typename Geometry_traits_2::Point_2 Point;

  //  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2> Traits_adapter_2;
  
public:
  /// Construct the viewer.
  /// @param a_arr the arrangement to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  SimpleArrangementViewerQtBase(QWidget *parent,
                                const char *title = "Basic Arrangement Viewer",
                                bool anofaces = false,
                                bool aisolatedvertices = true)
      : // First draw: vertices; edges, faces; multi-color; no inverse normal
        Base(parent, title, true, true, true, false, false),
        m_nofaces(anofaces), m_isolated_vertices(aisolatedvertices) {

    setKeyDescription(::Qt::Key_I, "Toggle view isolated vertices");

    compute_elements(); // to be overrided by child class
  }

protected:


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

  virtual void compute_elements() {}

protected:
  bool m_nofaces;
  bool m_isolated_vertices;
}; // class SimpleArrangementViewerQtBase

template <typename Arrangement_2>
class SimpleArrangementViewerQt
    : public SimpleArrangementViewerQtBase<Arrangement_2> {
public:
  SimpleArrangementViewerQt(QWidget *parent)
      : SimpleArrangementViewerQtBase<Arrangement_2>(parent) {}
};

template <typename Kernel_, typename Dcel>
    class SimpleArrangementViewerQt <
    CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel_>, Dcel>>
    : public SimpleArrangementViewerQtBase<
          CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel_>, Dcel>> {
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_segment_traits_2<Kernel> Traits;
  typedef SimpleArrangementViewerQtBase<CGAL::Arrangement_2<Traits, Dcel>> Superclass;

  typedef CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel_>, Dcel> Arrangement_2;

  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle Face_const_handle;
  typedef typename Arrangement_2::Edge_const_iterator Edge_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
      Ccb_halfedge_const_circulator;

  typedef typename Arrangement_2::Geometry_traits_2 Geometry_traits_2;
  typedef typename Arrangement_2::Topology_traits Topology_traits;

  typedef typename Geometry_traits_2::Point_2 Point;

  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2> Traits_adapter_2;

  typedef typename Superclass::Point_2 Point_2;
  typedef typename Superclass::Segment_2 Segment_2;
  typedef typename Superclass::Ray_2 Ray_2;
  typedef typename Superclass::Line_2 Line_2;
  typedef typename Superclass::Triangle_2 Triangle_2;
  typedef typename Superclass::Iso_rectangle_2 Iso_rectangle_2;
  typedef typename Superclass::Circle_2 Circle_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

public:
  SimpleArrangementViewerQt(QWidget *parent,
                            const Arrangement_2 &a_arr,
                            const char *title)
      : Superclass(parent, title), arr(a_arr)
  {
    compute_elements();
  }
public:

  void compute_face(Face_const_handle fh)
  {
    if (fh->is_unbounded())
    { return; }

    CGAL::Random random((unsigned long)(&*fh));
    CGAL::Color c=get_random_color(random);

    this->face_begin(c);

    Ccb_halfedge_const_circulator circ = fh->outer_ccb();
    Ccb_halfedge_const_circulator curr = circ;
    do {
      this->add_point_in_face(curr->source()->point());
    } while(++curr != circ);

    this->face_end();

    print_ccb (fh->outer_ccb());
    typename Arrangement_2::Hole_const_iterator hi;
    for (hi=fh->holes_begin(); hi!=fh->holes_end(); ++hi)
    { print_ccb (*hi); }
  }

  void compute_edge(Edge_const_iterator ei)
  {
    this->add_segment(ei->source()->point(), ei->target()->point());
    return;
  }

  void compute_elements()
  {
    this->clear();

    // Draw the arrangement vertices.
    typename Arrangement_2::Vertex_const_iterator vit;
    for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
    {
      if (vit->is_isolated())
      { // Isolated vertices are shown in black color
        //if (m_isolated_vertices) {
        //  this->add_point(vit->point(), CGAL::Color(CGAL::black()));
        //}
        //continue;
      }
      this->add_point(vit->point());
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

  void print_ccb(typename Arrangement_2::Ccb_halfedge_const_circulator circ) {
    typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
    //std::cout << "(" << curr->source()->point() << ")";
    do
    {
      Halfedge_const_handle he = curr;
      /*    std::cout << " [" << he->curve() << "] "
            << "(" << he->target()->point() << ")";*/
    } while (++curr != circ);
  }

private:
  const Arrangement_2 &arr;
};

template <typename CircularKernel, typename Dcel>
class SimpleArrangementViewerQt<
    CGAL::Arrangement_2<CGAL::Arr_circle_segment_traits_2<CircularKernel>, Dcel>>
    : public SimpleArrangementViewerQtBase<CGAL::Arrangement_2<
          CGAL::Arr_circle_segment_traits_2<CircularKernel>, Dcel>> {
public:
  typedef CircularKernel Kernel;
  typedef CGAL::Arr_circle_segment_traits_2<Kernel> Traits;
  typedef SimpleArrangementViewerQtBase<CGAL::Arrangement_2<Traits, Dcel>>
      Superclass;

  typedef CGAL::Arrangement_2<CGAL::Arr_circle_segment_traits_2<CircularKernel>, Dcel>
      Arrangement_2;

  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle Face_const_handle;
  typedef typename Arrangement_2::Edge_const_iterator Edge_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
      Ccb_halfedge_const_circulator;

  typedef typename Arrangement_2::Geometry_traits_2 Geometry_traits_2;
  typedef typename Arrangement_2::Topology_traits Topology_traits;

  typedef typename Geometry_traits_2::Point_2 Point;

  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2> Traits_adapter_2;

  typedef typename Superclass::Point_2 Point_2;
  typedef typename Superclass::Segment_2 Segment_2;
  typedef typename Superclass::Ray_2 Ray_2;
  typedef typename Superclass::Line_2 Line_2;
  typedef typename Superclass::Triangle_2 Triangle_2;
  typedef typename Superclass::Iso_rectangle_2 Iso_rectangle_2;
  typedef typename Superclass::Circle_2 Circle_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  typedef CGAL::Exact_predicates_exact_constructions_kernel Viewer_kernel;
  typedef typename Kernel::FT NT;

public:
  SimpleArrangementViewerQt(QWidget *parent, const Arrangement_2 &a_arr,
                            const char *title)
      : Superclass(parent, title), arr(a_arr) {
    compute_elements();
  }

public:
  void compute_face(Face_const_handle fh) {
    if (fh->is_unbounded()) {
      return;
    }

    CGAL::Random random((unsigned long)(&*fh));
    CGAL::Color c = get_random_color(random);

    //this->face_begin(c);

    Ccb_halfedge_const_circulator circ = fh->outer_ccb();
    Ccb_halfedge_const_circulator curr = circ;
    do {
      //this->add_point_in_face(curr->source()->point());
    } while (++curr != circ);

    //this->face_end();

    print_ccb(fh->outer_ccb());
    typename Arrangement_2::Hole_const_iterator hi;
    for (hi = fh->holes_begin(); hi != fh->holes_end(); ++hi) {
      print_ccb(*hi);
    }
  }

  void compute_edge(Edge_const_iterator ei) {

    Viewer_kernel::Point_3 source(to_double(ei->source()->point().x()), 0,
                               to_double(ei->source()->point().y())),
        target(to_double(ei->target()->point().x()), 0,
            to_double(ei->target()->point().y()));

    if (ei->curve().is_linear()) {
      this->add_segment(source, target);
    } else {
      CGAL_precondition(ei->curve().is_circular());
      const Circle_2& circ = ei->curve().supporting_circle();
      Viewer_kernel::Point_3 center(
          to_double(circ.center().x()), 0., to_double(circ.center().y()));
      double radius = to_double(std::sqrt(to_double(circ.squared_radius())));

      this->add_ellipse(center, 0., CGAL_PI, radius, radius);
    }
    return;
  }

  void compute_elements() {
    this->clear();

    // Draw the arrangement vertices.
    typename Arrangement_2::Vertex_const_iterator vit;
    for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
      if (vit->is_isolated()) { // Isolated vertices are shown in black color
                                // if (m_isolated_vertices) {
                                //  this->add_point(vit->point(),
                                //  CGAL::Color(CGAL::black()));
                                //}
                                // continue;
      }
      // maybe convert one root point 2 to a local point before calling this method
      Viewer_kernel::Point_3 p(to_double(vit->point().x()), 0,
                               to_double(vit->point().y()));

      this->add_point(p);
    }

    // Draw the arrangement edges.
    typename Arrangement_2::Edge_const_iterator eit;
    for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
      compute_edge(eit);
      std::cout << "[" << eit->curve() << "]" << std::endl;
    }

    // Draw the arrangement faces.
    typename Arrangement_2::Face_const_iterator fit;
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
      compute_face(fit);
    }
  }

  void print_ccb(typename Arrangement_2::Ccb_halfedge_const_circulator circ) {
    typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
    // std::cout << "(" << curr->source()->point() << ")";
    do {
      Halfedge_const_handle he = curr;
      /*    std::cout << " [" << he->curve() << "] "
            << "(" << he->target()->point() << ")";*/
    } while (++curr != circ);
  }

private:
  const Arrangement_2 &arr;
};

template<class GeomTraits_, class TopTraits_>
void draw(const Arrangement_2<GeomTraits_, TopTraits_>& a_arr,
          const char* title="Basic Arrangement Viewer")
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
      mainwindow(app.activeWindow(), a_arr, title);
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_ARRANGEMENT_2_H
