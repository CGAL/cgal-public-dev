// Copyright (c) 2012  Tel-Aviv University (Israel).
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
// $URL: $
// $Id: $
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#ifndef CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
#define CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H

#include <QObject>
#include <QRectF>
#include <QPainterPath>
#include <vector>

// TODO: should be included in PainterOstream.h
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>

#include "Utils.h"

#include <CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h>

class QPainter;

namespace CGAL {
namespace Qt {

template < typename ArrTraits >
class ArrangementPainterOstreamBase : public QGraphicsSceneMixin
{
public:
  // typedefs
  typedef ArrTraits Traits;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Kernel::Ray_2                        Ray_2;
  typedef typename Kernel::Line_2                       Line_2;
  typedef typename Kernel::Triangle_2                   Triangle_2;
  typedef typename Kernel::Iso_rectangle_2              Iso_rectangle_2;
  typedef typename Kernel::Circle_2                     Circle_2;

public:
  /*! Constructor */
  ArrangementPainterOstreamBase( QPainter* p,
                                 QRectF clippingRectangle = QRectF( ) ) :
    painterOstream( p, clippingRectangle ),
    qp( p ),
    convert( clippingRectangle ),
    // scene( NULL ),
    clippingRect( QRectF( ) ), // null rectangle
    scale( 1.0 )
  {
    if ( p != 0 )
    {
      this->scale = p->worldTransform( ).m11( );
    }
  }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstreamBase() {}

  // methods
  template < typename T >
  ArrangementPainterOstreamBase& operator<<( const T& t )
  {
    this->painterOstream << t;
    return *this;
  }

  void setScene( QGraphicsScene* scene_ )
  {
    this->scene = scene_;

    // set the clipping rectangle
    if ( scene_ == NULL )
    {
      return;
    }
    this->clippingRect = this->viewportRect( );
    this->convert = Converter< Kernel >( this->clippingRect );
  }

  void drawPoint( const QPointF& qpt )
  {
    QPen savePen = this->qp->pen( );
    this->qp->setBrush( QBrush( savePen.color( ) ) );
    double radius = savePen.width( ) / 2.0;
    radius /= this->scale;

    this->qp->drawEllipse( qpt, radius, radius );

    this->qp->setBrush( QBrush( ) );
    this->qp->setPen( savePen );
  }

protected: // methods

protected: // fields
  PainterOstream< Kernel > painterOstream;
  QPainter* qp;
  Converter< Kernel > convert;
  QRectF clippingRect;
  double scale;
}; // class ArrangementPainterOstreamBase

template < typename ArrTraits >
class ArrangementPainterOstream :
    public ArrangementPainterOstreamBase< ArrTraits >
{
public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()):
    ArrangementPainterOstreamBase< ArrTraits >( p, clippingRectangle )
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}
};

template < typename Kernel_ >
class ArrangementPainterOstream<CGAL::Arr_segment_traits_2< Kernel_> >:
  public ArrangementPainterOstreamBase<CGAL::Arr_segment_traits_2<Kernel_> >
{
public: // typedefs
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_segment_traits_2< Kernel > Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()) :
    Superclass(p, clippingRectangle)
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
  {
    const Point_2& p1 = curve.source( );
    const Point_2& p2 = curve.target( );
    Segment_2 seg( p1, p2 );

    // skip segments outside our view
    QRectF seg_bb = this->convert( seg.bbox( ) );
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.intersects( seg_bb ) )
    {
      return *this;
    }

    this->painterOstream << seg;
    return *this;
  }

  ArrangementPainterOstream& operator<<( const Point_2& p )
  {
    QPointF qpt = this->convert( p );
    // clip the point if possible
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.contains( qpt ) )
    {
      return *this;
    }

    this->drawPoint( qpt );
    return *this;
  }

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
};

template < typename SegmentTraits >
class ArrangementPainterOstream<CGAL::Arr_polyline_traits_2<SegmentTraits> > :
  public ArrangementPainterOstreamBase<CGAL::Arr_polyline_traits_2<
                                         SegmentTraits> >
{
public: // typedefs
  typedef ArrangementPainterOstreamBase<CGAL::Arr_polyline_traits_2<
                                          SegmentTraits> > Superclass;
  typedef typename Superclass::Traits                   Traits;
  typedef typename Superclass::Kernel                   Kernel;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()):
    Superclass( p, clippingRectangle )
  { }

    /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
  {
    for ( unsigned int i = 0; i < curve.size( ); ++i )
    {
      Segment_2 segment = curve[ i ];
      this->painterOstream << segment;
    }
    // TODO: implement polyline painting
#if 0
    const Point_2& p1 = curve.source( );
    const Point_2& p2 = curve.target( );
    Segment_2 seg( p1, p2 );
    this->painterOstream << seg;
#endif
    return *this;
  }

  // cloned from segtraits painter
  ArrangementPainterOstream& operator<<( const Point_2& p )
  {
    QPointF qpt = this->convert( p );
    this->drawPoint( qpt );
    return *this;
  }

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
};

template < typename RatKernel, class AlgKernel, class NtTraits >
class ArrangementPainterOstream<CGAL::Arr_conic_traits_2<RatKernel, AlgKernel,
                                                         NtTraits > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_conic_traits_2<RatKernel,
                                                                 AlgKernel,
                                                                 NtTraits> >
{
public: // typedefs
  typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Construct_x_monotone_curve_2
    Construct_x_monotone_curve_2;
  typedef typename Traits::Point_2                      Intersection_point_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::FT                           FT;

public: // inner classes
  // utility class to use with std::sort on an Intersect_2 result set.
  class Compare_intersection_point_result
  {
  public:
    typedef std::pair< Intersection_point_2, Multiplicity > Result;
    // returns whether the point1 < point2, using x-coord to compare
    bool operator()( const Result& o1, const Result& o2 )
    {
      Point_2 p1 = o1.first;
      Point_2 p2 = o2.first;
      return ( p1.x( ) < p2.x( ) );
    }
  };

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()) :
    Superclass( p, clippingRectangle ),
    //intersect_2( this->traits.intersect_2_object( ) ),
    // Why doesn't this work?
    construct_x_monotone_curve_2(this->
                                 traits.construct_x_monotone_curve_2_object())
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
  {
    CGAL::Bbox_2 bb = curve.bbox( );
    QRectF qbb = this->convert( bb );
    // quick cull
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.intersects( qbb ) )
    {
      //std::cout << "quick culled curve" << std::endl;
      return *this;
    }

#if 0
    std::cout << "bottom: ("
              << this->clippingRect.bottomLeft( ).x( )
              << " "
              << this->clippingRect.bottomLeft( ).y( )
              << " "
              << this->clippingRect.bottomRight( ).x( )
              << " "
              << this->clippingRect.bottomRight( ).y( )
              << ")"
              << std::endl;
#endif

    if ( this->clippingRect.isValid( ) )
    {
      std::vector< X_monotone_curve_2 > visibleParts;
      if ( this->clippingRect.contains( qbb ) )
      {
        visibleParts.push_back( curve );
      }
      else
      {
        visibleParts = this->visibleParts( curve );
      }
      for ( unsigned int i = 0; i < visibleParts.size( ); ++i )
      {
        X_monotone_curve_2 subcurve = visibleParts[ i ];
        int n;
        if ( this->scene == NULL )
          n = 100; // TODO: get an adaptive approximation
        else
        {
          QGraphicsView* view = this->scene->views( ).first( );
          int xmin, xmax;
          xmin = view->mapFromScene( bb.xmin( ), bb.ymin( ) ).x( );
          xmax = view->mapFromScene( bb.xmax( ), bb.ymin( ) ).x( );
          n = xmax - xmin;
        }
        if ( n == 0 )
        {
          return *this;
        }

        std::pair<double, double>* app_pts =
          new std::pair<double, double>[n + 1];
        std::pair<double, double>* end_pts =
          subcurve.polyline_approximation(n, app_pts);
        std::pair<double, double>* p_curr = app_pts;
        std::pair<double, double>* p_next = p_curr + 1;
        int count = 0;
        do
        {
          QPointF p1( p_curr->first, p_curr->second );
          QPointF p2( p_next->first, p_next->second );
#if 0
          Segment_2 seg( p1, p2 );
          this->painterOstream << seg;
#endif
          this->qp->drawLine( p1, p2 );
          p_curr++;
          p_next++;
          ++count;
        }
        while ( p_next != end_pts );
      }
    }
    else
    { // draw the whole curve
      int n;
      if ( this->scene == NULL )
        n = 100; // TODO: get an adaptive approximation
      else
      {
        QGraphicsView* view = this->scene->views( ).first( );
        int xmin, xmax;
        xmin = view->mapFromScene( bb.xmin( ), bb.ymin( ) ).x( );
        xmax = view->mapFromScene( bb.xmax( ), bb.ymin( ) ).x( );
        n = xmax - xmin;
      }
      if ( n == 0 )
      {
        return *this;
      }

      std::pair<double, double>* app_pts = new std::pair<double, double>[n + 1];
      std::pair<double, double>* end_pts =
        curve.polyline_approximation(n, app_pts);
      std::pair<double, double>* p_curr = app_pts;
      std::pair<double, double>* p_next = p_curr + 1;
      int count = 0;
      do
      {
        QPointF p1( p_curr->first, p_curr->second );
        QPointF p2( p_next->first, p_next->second );
#if 0
        Segment_2 seg( p1, p2 );
        this->painterOstream << seg;
#endif
        this->qp->drawLine( p1, p2 );
        p_curr++;
        p_next++;
        ++count;
      }
      while ( p_next != end_pts );
      //std::cout << count << " approximation points" << std::endl;
    }

    return *this;
  }

  // cloned from segtraits painter
  ArrangementPainterOstream& operator<<( const Point_2& p )
  {
    QPointF qpt = this->convert( p );
    this->drawPoint( qpt );
    return *this;
  }

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }

protected: // methods
  // Returns subcurves of curve that are actually visible in the view.
  // Assumes that clippingRect is valid.
  std::vector< X_monotone_curve_2 > visibleParts( X_monotone_curve_2 curve )
  {
    // see if we intersect the bottom edge of the viewport
    Intersect_2 intersect_2 = this->traits.intersect_2_object( );
    Point_2 bottomLeft = this->convert( this->clippingRect.bottomLeft( ) );
    Point_2 bottomRight = this->convert( this->clippingRect.bottomRight( ) );
    Point_2 topLeft = this->convert( this->clippingRect.topLeft( ) );
    Point_2 topRight = this->convert( this->clippingRect.topRight( ) );
    X_monotone_curve_2 bottom =
      this->construct_x_monotone_curve_2( bottomLeft, bottomRight );
    X_monotone_curve_2 left =
      this->construct_x_monotone_curve_2( bottomLeft, topLeft );
    X_monotone_curve_2 top =
      this->construct_x_monotone_curve_2( topLeft, topRight );
    X_monotone_curve_2 right =
      this->construct_x_monotone_curve_2( topRight, bottomRight );

    std::vector< CGAL::Object > bottomIntersections;
    std::vector< CGAL::Object > leftIntersections;
    std::vector< CGAL::Object > topIntersections;
    std::vector< CGAL::Object > rightIntersections;
    std::vector< CGAL::Object > intersections;

    intersect_2( bottom, curve, std::back_inserter( bottomIntersections ) );
    intersect_2( left, curve, std::back_inserter( leftIntersections ) );
    intersect_2( top, curve, std::back_inserter( topIntersections ) );
    intersect_2( right, curve, std::back_inserter( rightIntersections ) );
    // int total = bottomIntersections.size( )
    //   + leftIntersections.size( )
    //   + topIntersections.size( )
    //   + rightIntersections.size( );

    intersect_2( bottom, curve, std::back_inserter( intersections ) );
    intersect_2( left, curve, std::back_inserter( intersections ) );
    intersect_2( top, curve, std::back_inserter( intersections ) );
    intersect_2( right, curve, std::back_inserter( intersections ) );

    this->filterIntersectionPoints( intersections );
    //std::cout << "total intersections: " << intersections.size( )
    //          << std::endl;
    //this->printIntersectResult( intersections );

    Point_2 leftEndpt = curve.source( );
    Point_2 rightEndpt = curve.target( );
    if ( leftEndpt.x( ) > rightEndpt.x( ) )
    {
      std::swap( leftEndpt, rightEndpt );
    }
    QPointF qendpt1 = this->convert( leftEndpt );
    QPointF qendpt2 = this->convert( rightEndpt );
    std::list< Point_2 > pointList;
    for ( unsigned int i = 0; i < intersections.size( ); ++i )
    {
      CGAL::Object o = intersections[ i ];
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, o ) )
      {
        Point_2 pt = pair.first;
        pointList.push_back( pt );
      }
    }
    bool includeLeftEndpoint = this->clippingRect.contains( qendpt1 );
    bool includeRightEndpoint = this->clippingRect.contains( qendpt2 );
    if ( includeLeftEndpoint )
      pointList.push_front( leftEndpt );
    if ( includeRightEndpoint )
      pointList.push_back( rightEndpt );

    Construct_x_monotone_subcurve_2< Traits > construct_x_monotone_subcurve_2;
    std::vector< X_monotone_curve_2 > clippings;
    typename std::list< Point_2 >::iterator pointListItr = pointList.begin( );
    for ( unsigned int i = 0; i < pointList.size( ); i += 2 )
    {
      Point_2 p1 = *pointListItr++;
      Point_2 p2 = *pointListItr++;
      X_monotone_curve_2 subcurve =
        construct_x_monotone_subcurve_2( curve, p1, p2 );
      clippings.push_back( subcurve );
    }

#if 0
    // std::cout << "pointList size: " << pointList.size( ) << std::endl;
    // if ( intersections.size( ) % 2 == 0 )
    // {
    //   // either both curve endpoints are in view or both are out
    //     if ( this->clippingRect.contains( qendpt1 ) )
    //     {
    //         if ( this->clippingRect.contains( qendpt2 ) )
    //         {
    //             std::cout << "both endpoints are in view" << std::endl;
    //         }
    //     }
    //     else if ( !this->clippingRect.contains( qendpt2 ) )
    //     {
    //         std::cout << "both endpoints are out of view" << std::endl;
    //     }
    // }
    // else
    // { // one curve endpoint is in view
    //     if ( this->clippingRect.contains( qendpt1 ) )
    //     {
    //         std::cout << "left endpoint is in view" << std::endl;
    //     }
    //     else if ( this->clippingRect.contains( qendpt2 ) )
    //     {
    //         std::cout << "right endpoint is in view" << std::endl;
    //     }
    // }

    std::vector< X_monotone_curve_2 > res;
    res.push_back( curve );
    return res;
#endif
    return clippings;
  }

  // keep only the intersection points ie. throw out overlapping curve segments
  void filterIntersectionPoints( std::vector< CGAL::Object >& res )
  {
    std::vector< std::pair< Intersection_point_2, Multiplicity > > tmp;

    // filter out the non-intersection point results
    for ( unsigned int i = 0; i < res.size( ); ++i )
    {
      CGAL::Object obj = res[ i ];
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, obj ) )
      {
        tmp.push_back( pair );
      }
    }
    res.clear( );

    // sort the intersection points by x-coord
    Compare_intersection_point_result compare_intersection_point_result;
    std::sort( tmp.begin( ), tmp.end( ), compare_intersection_point_result );

    // box up the sorted elements
    for ( unsigned int i = 0; i < tmp.size( ); ++i )
    {
      std::pair< Intersection_point_2, Multiplicity > pair = tmp[ i ];
      CGAL::Object o = CGAL::make_object( pair );
      res.push_back( o );
    }
  }

  void printIntersectResult( const std::vector< CGAL::Object >& res )
  {
    for ( std::vector< CGAL::Object >::const_iterator it = res.begin( );
          it != res.end( ); ++it )
    {
      CGAL::Object obj = *it;
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, obj ) )
      {
        Point_2 pt = pair.first;
        QPointF qpt = this->convert( pt );
        // std::cout << "(" << pt.x( ) << " " << pt.y( ) < ")" << std::endl;
      }
    }
  }

protected: // members
  Traits traits;
  //Intersect_2 intersect_2;
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2;
};

template < typename Kernel_ >
class ArrangementPainterOstream< CGAL::Arr_linear_traits_2< Kernel_ > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_linear_traits_2< Kernel_ > >
{
public: // typedefs
  typedef Kernel_                                       Kernel;
  typedef CGAL::Arr_linear_traits_2< Kernel >           Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()) :
    Superclass( p, clippingRectangle )
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
  {
    if ( curve.is_segment( ) )
    {
      Segment_2 seg = curve.segment( );

      // skip segments outside our view
      QRectF seg_bb = this->convert( seg.bbox( ) );
      if ( this->clippingRect.isValid( ) &&
           ! this->clippingRect.intersects( seg_bb ) )
      {
        return *this;
      }

      this->painterOstream << seg;
    }
    else if ( curve.is_ray( ) )
    {
      Ray_2 ray = curve.ray( );
      QLineF qseg = this->convert( ray );
      if ( qseg.isNull( ) )
      { // it's out of view
        return *this;
      }
      Segment_2 seg = this->convert( qseg );
      this-> painterOstream << seg;
    }
    else // curve.is_line( )
    {
      Line_2 line = curve.line( );
      QLineF qseg = this->convert( line );
      if ( qseg.isNull( ) )
      { // it's out of view
        return *this;
      }
      Segment_2 seg = this->convert( qseg );
      this-> painterOstream << seg;
    }
    return *this;
  }

  ArrangementPainterOstream& operator<<( const Point_2& p )
  {
    QPointF qpt = this->convert( p );
    // clip the point if possible
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.contains( qpt ) )
    {
      return *this;
    }

    this->drawPoint( qpt );
    return *this;
  }

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
};

template < typename CircularKernel >
class ArrangementPainterOstream< CGAL::Arr_circular_arc_traits_2<
                                   CircularKernel > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_circular_arc_traits_2<
                                          CircularKernel > >
{
public:
  typedef CircularKernel                                Kernel;
  typedef CGAL::Arr_circular_arc_traits_2< Kernel >     Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()):
    Superclass( p, clippingRectangle )
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
  {
    this->painterOstream << curve;
    return *this;
  }

  ArrangementPainterOstream& operator<<( const Point_2& p )
  {
    QPointF qpt = this->convert( p );
    // clip the point if possible
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.contains( qpt ) )
    {
      return *this;
    }

    this->drawPoint( qpt );
    return *this;
  }

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
};

/**
Specialization of ArrangementPainterOstream for Bezier curve traits.
*/
template < typename RatKernel, class AlgKernel, class NtTraits >
class ArrangementPainterOstream<CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel,
                                                         NtTraits > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_Bezier_curve_traits_2<
                                          RatKernel, AlgKernel, NtTraits > >
{
public:
  typedef CGAL::Arr_Bezier_curve_traits_2< RatKernel, AlgKernel, NtTraits >
    Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Point_2 ArrPointType;
  typedef typename RatKernel::Point_2 Rat_point_2;
  typedef std::map< Curve_2, double,
    Compare_Bezier_curve_2< RatKernel, AlgKernel, NtTraits > >
    Curve_2_double_map;

protected:
  Converter< RatKernel > m_ratConverter;
  Curve_2_double_map m_sampleFactors;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()):
    Superclass( p, clippingRectangle ),
    m_ratConverter( clippingRectangle )
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const Curve_2& curve )
  {
#if 0 // only works for cubic splines
    std::vector< QPointF > controlPts;
    for ( int i = 0; i < curve.number_of_control_points(); ++i )
    {
      Rat_point_2 controlPt = curve.control_point( i );
      QPointF qpt = m_ratConverter( controlPt );
      controlPts.push_back( qpt );
    }
    if (controlPts.size() < 4)
      return *this;
    std::cout << "control points: " <<
      curve.number_of_control_points() << "\n";

    QPainterPath path;
    path.moveTo(controlPts.at(0));
    path.cubicTo(controlPts.at(1),
      controlPts.at(2),
      controlPts.at(3));
    this->qp->strokePath(path, this->qp->pen( ) );
#endif
    double P = 2 * (curve.number_of_control_points() + 1);
    // TODO: Calculate this number appropriately
    int N = P;

    std::vector< std::pair<double, double> > samples;
    curve.sample(0, 1, N, std::back_inserter( samples ) );
    for (int i = 0; i < samples.size() - 1; ++i)
    {
      QPointF p1( samples[i].first, samples[i].second );
      QPointF p2( samples[i+1].first, samples[i+1].second );
      this->qp->drawLine( p1, p2 );
    }

    return *this;
  }

  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
  {
    Curve_2 supportingCurve = curve.supporting_curve();
    std::pair< double, double > range = curve.parameter_range( );
    double P = 2 * (supportingCurve.number_of_control_points() + 1);
    int N = P;

    // TODO: Scale the number of samples by the norm of the view matrix
    double alpha = estimate_alpha( curve );
    const double Epsilon = 0.01;
    const int Factor = P / 2 * sqrt( alpha / Epsilon );
    N = ( Factor > N )? Factor : N;

    std::vector< std::pair<double, double> > samples;
    supportingCurve.sample(range.first, range.second, N, std::back_inserter( samples ) );
    for (int i = 0; i < samples.size() - 1; ++i)
    {
      QPointF p1( samples[i].first, samples[i].second );
      QPointF p2( samples[i+1].first, samples[i+1].second );
      this->qp->drawLine( p1, p2 );
    }

    return *this;
  }

  double estimate_alpha( const X_monotone_curve_2& curve )
  {
    typedef CGAL::Simple_cartesian< double > App_kernel;
    typedef App_kernel::Point_2 App_point_2;
    typedef App_kernel::Segment_2 App_segment_2;

    std::vector< std::pair<double, double> > samples;
    std::pair<double, double> range = curve.parameter_range( );
    Curve_2 supportingCurve = curve.supporting_curve();
    int P = 2 * (supportingCurve.number_of_control_points() + 1);
    supportingCurve.sample(range.first, range.second,
      2*P + 1, // number of samples
      std::back_inserter( samples ) );

    std::vector< App_point_2 > query_pts;
    std::vector< App_segment_2 > segments;
    double max_error = 0.0;
    for (int i = 0; i < P; ++i)
    {
      App_point_2 p1( samples[2*i].first, samples[2*i].second );
      App_point_2 p2( samples[2*i+2].first, samples[2*i+2].second );
      App_segment_2 seg( p1, p2 );
      segments.push_back( seg );
    }

    for (int i = 0; i < P; ++i)
    {
      App_point_2 query_pt( samples[2*i+1].first, samples[2*i+1].second );
      double error = CGAL::squared_distance( query_pt, segments[i] );
      if ( error > max_error )
        max_error = error;
    }

    return sqrt(max_error);
  }

  ArrangementPainterOstream& operator<<( const
    typename Construct_x_monotone_subcurve_2< Traits >::Subcurve& subcurve )
  {
    X_monotone_curve_2 curve = subcurve.m_cv;
    Curve_2 supportingCurve = curve.supporting_curve();
    std::pair< double, double > range;
    range.first = subcurve.m_t1;
    range.second = subcurve.m_t2;
    double P = 2 * (supportingCurve.number_of_control_points() + 1);
    int N = P;

    // TODO: Scale the number of samples by the norm of the view matrix
    double alpha = estimate_alpha( curve );
    const double Epsilon = 0.01;
    const int Factor = P / 2 * sqrt( alpha / Epsilon );
    N = ( Factor > N )? Factor : N;

    std::vector< std::pair<double, double> > samples;
    supportingCurve.sample(range.first, range.second, N, std::back_inserter( samples ) );
    for (int i = 0; i < samples.size() - 1; ++i)
    {
      QPointF p1( samples[i].first, samples[i].second );
      QPointF p2( samples[i+1].first, samples[i+1].second );
      this->qp->drawLine( p1, p2 );
    }

    return *this;
  }

  ArrangementPainterOstream& operator<<( const ArrPointType& pt )
  {
    CORE::BigRat x_min, x_max, y_min, y_max;
    pt.get_bbox( x_min, y_min, x_max, y_max );
    Point_2 approx( CGAL::to_double((x_min + x_max)/2), CGAL::to_double((y_min + y_max)/2) );
    //std::cout << "pt("
    //  << approx.x() << ", "
    //  << approx.y() << ")\n";
    QPointF qpt = this->convert( approx );
    // clip the point if possible
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.contains( qpt ) )
    {
      return *this;
    }

    this->drawPoint( qpt );
    return *this;
  }

  ArrangementPainterOstream& operator<<( const Point_2& p )
  {
    QPointF qpt = this->convert( p );
    // clip the point if possible
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.contains( qpt ) )
    {
      return *this;
    }

    this->drawPoint( qpt );
    return *this;
  }

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
};

#if 1
template < typename Coefficient_ >
class ArrangementPainterOstream< CGAL::Arr_algebraic_segment_traits_2<
                                   Coefficient_ > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_algebraic_segment_traits_2<
                                          Coefficient_ > >
{
public:
  typedef Coefficient_                                  Coefficient;
  typedef typename CGAL::Arr_algebraic_segment_traits_2< Coefficient >
                                                        Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Traits::CKvA_2                       CKvA_2;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

  typedef Curve_renderer_facade<CKvA_2> Facade;
  typedef std::pair< int, int > Coord_2;
  typedef std::vector< Coord_2 > Coord_vec_2;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()):
    Superclass( p, clippingRectangle ),
    m_width( 0 ),
    m_height( 0 )
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods

  /**
  Rasterize an algebraic curve.
  */
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
  {
    this->setupFacade( );

    std::list< Coord_vec_2 > points;
    this->rasterize( curve, &points );

    // stroke a path from the raster points
    QGraphicsView* view = this->scene->views( ).first( );
    for ( typename std::list<Coord_vec_2>::const_iterator lit = points.begin();
      lit != points.end(); ++lit )
    {
      const Coord_vec_2& vec = *lit;
      typename Coord_vec_2::const_iterator vit = vec.begin();

      QPainterPath path;
      QPoint coord( vit->first, vit->second );
      QPointF qpt = view->mapToScene( coord );
      if ( vit != vec.end() )
      {
        path.moveTo( qpt );
      }
      while ( vit != vec.end() )
      {
        path.lineTo( qpt );
        vit++;
        coord = QPoint( vit->first, vit->second );
        qpt = view->mapToScene( coord );
      }
      this->qp->drawPath( path );
    }

    return *this;
  }

  /**
  Draw an algebraic point.
  */
  ArrangementPainterOstream& operator<<( const Point_2& p )
  {
    this->setupFacade( );

    // draw an ellipse at the point
    QGraphicsView* view = this->scene->views( ).first( );
    QPoint coords;
    if ( this->rasterize( p, &coords ) )
    {
      QPointF qpt = view->mapToScene( coords );
      this->drawPoint( qpt );
    }

    return *this;
  }

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }

protected:
  void setupFacade( )
  {
    typedef Curve_renderer_facade<CKvA_2> Facade;
    QGraphicsView* view = this->scene->views( ).first( );
    QRectF viewport = this->viewportRect( );
    CGAL::Bbox_2 bbox = this->convert( viewport ).bbox( );
    Facade::setup(bbox, view->width(), view->height());
    m_width = view->width();
    m_height = view->height();
  }

  // reflect across the equator of the viewport
  void flipY( std::list< Coord_vec_2 >* points )
  {
    QGraphicsView* view = this->scene->views( ).first( );
    QMatrix flipMat( 1, 0, 0, -1, 0, 0 );

    for ( std::list< Coord_vec_2 >::iterator it = points->begin( );
      it != points->end( ); ++it )
    {
      for ( Coord_vec_2::iterator pit = it->begin( );
        pit != it->end( ); ++pit )
      {
        Coord_vec_2& point = *it;
        Coord_2 before = *pit;

        Coord_2 after( before.first, m_height - before.second );
        *pit = after;
      }
    }
  }

  // reflect across the equator of the viewport
  QPoint flipY( const QPoint& pt )
  {
    return QPoint( pt.x(), m_height - pt.y() );
  }

  void rasterize( const X_monotone_curve_2& curve,
    std::list< Coord_vec_2 >* points )
  {
    // rasterize the curve
    boost::optional < Coord_2 > p1, p2;
    Facade::instance().draw( curve, *points, &p1, &p2 );

    // raster points are flipped wrt y-axis, so flip them.
    flipY( points );
  }

  bool rasterize( const Point_2& point,
    QPoint* coords )
  {
    Coord_2 coord;
    bool ok = Facade::instance().draw( point, coord );
    if ( ok )
    {
      QPoint tmp( coord.first, coord.second );
      // raster point is flipped wrt y-axis, so flip it.
      *coords = flipY( tmp );
    }
    return ok;
  }

protected:
  int m_width;
  int m_height;

};
#endif

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
