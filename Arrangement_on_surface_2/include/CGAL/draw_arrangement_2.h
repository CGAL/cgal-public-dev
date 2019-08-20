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
#include <CGAL/Qt/Converter.h>

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>

#include <CGAL/Arr_circular_arc_traits_2.h>
#include <CGAL/Arr_circular_arc_traits_2.h>.h>

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

template <class CircularKernel>
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

template < class ArrTraits >
class Arr_compute_y_at_x_2 {
public:
  typedef ArrTraits Traits;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename ArrTraitsAdaptor< Traits >::CoordinateType CoordinateType;
  // typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Line_2                       Line_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef std::pair< typename Traits::Point_2, Multiplicity >
                                                        IntersectionResult;

  /*! Constructor */
  Arr_compute_y_at_x_2( ) :
    intersectCurves( this->traits.intersect_2_object( ) )
  { }

  /*! Destructor (virtual) */
  virtual ~Arr_compute_y_at_x_2() {}

  CoordinateType operator() ( const X_monotone_curve_2& curve,
                              const CoordinateType& x )
  {
    typename Traits::Left_side_category category;
    return this->operator()( curve, x, this->traits, category );
  }

  double approx( const X_monotone_curve_2& curve, const CoordinateType& x )
  {
    return CGAL::to_double( (*this)( curve, x ) );
  }

protected:
  template < class TTraits >
  CoordinateType operator() ( const X_monotone_curve_2& curve,
                              const CoordinateType& x, TTraits traits_,
                              CGAL::Arr_oblivious_side_tag )
  {
    typedef typename TTraits::Construct_x_monotone_curve_2
      Construct_x_monotone_curve_2;
    Construct_x_monotone_curve_2 construct_x_monotone_curve_2 =
      traits_.construct_x_monotone_curve_2_object( );
    CoordinateType res( 0 );
    CGAL::Bbox_2 clipRect = curve.bbox( );
    Point_2 p1c1( x, CoordinateType( clipRect.ymin( ) - 1 ) ); // clicked point
    // upper bounding box
    Point_2 p2c1( x, CoordinateType( clipRect.ymax( ) + 1 ) );

    const X_monotone_curve_2 verticalLine =
      construct_x_monotone_curve_2( p1c1, p2c1 );
    CGAL::Object o;
    CGAL::Oneset_iterator< CGAL::Object > oi( o );

    this->intersectCurves( curve, verticalLine, oi );

    IntersectionResult pair;
    if ( CGAL::assign( pair, o ) )
    {
      Point_2 pt = pair.first;
      res = pt.y( );
    }
    return res;
  }

  template < class TTraits >
  CoordinateType operator() ( const X_monotone_curve_2& curve,
                              const CoordinateType& x, TTraits traits_,
                              CGAL::Arr_open_side_tag )
  {
    typename TTraits::Construct_x_monotone_curve_2
      construct_x_monotone_curve_2 =
      traits_.construct_x_monotone_curve_2_object( );
    CoordinateType res( 0 );
    // QRectF clipRect = this->viewportRect( );
    Line_2 line = curve.supporting_line( );
    // FIXME: get a better bounding box for an unbounded segment
    Point_2 p1c1( x, CoordinateType( -10000000 ) ); // clicked point
    Point_2 p2c1( x, CoordinateType(  10000000 ) ); // upper bounding box

    const X_monotone_curve_2 verticalLine =
      construct_x_monotone_curve_2( p1c1, p2c1 );
    CGAL::Object o;
    CGAL::Oneset_iterator< CGAL::Object > oi( o );

    this->intersectCurves( curve, verticalLine, oi );

    IntersectionResult pair;
    if ( CGAL::assign( pair, o ) )
    {
      Point_2 pt = pair.first;
      res = pt.y( );
    }
    return res;
  }

protected:
  Traits traits;
  Intersect_2 intersectCurves;
};

template < class CircularKernel >
class Arr_compute_y_at_x_2< CGAL::Arr_circular_arc_traits_2<CircularKernel> > {
public:
  typedef CGAL::Arr_circular_arc_traits_2< CircularKernel > Traits;
  typedef CircularKernel                                Kernel;
  typedef typename Kernel::FT                           FT;
  typedef typename Kernel::Root_of_2                    Root_of_2;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Traits::Point_2                      Arc_point_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Kernel::Line_arc_2                   Line_arc_2;
  // Circular_arc_2
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef std::pair< typename Traits::Point_2, Multiplicity >
                                                        IntersectionResult;

  /*! Constructor */
  Arr_compute_y_at_x_2( ) :
    intersectCurves( this->traits.intersect_2_object( ) )
  { }

  /*! Destructor (virtual) */
  virtual ~Arr_compute_y_at_x_2() {}

  Root_of_2 operator() ( const X_monotone_curve_2& curve, const FT& x )
  {
    Root_of_2 res( 0 );
    CGAL::Bbox_2 clipRect = curve.bbox( );
    Point_2 p1c1( x, FT( clipRect.ymin( ) - 1 ) ); // clicked point
    Point_2 p2c1( x, FT( clipRect.ymax( ) + 1 ) ); // upper bounding box
    Line_arc_2 verticalLine( Segment_2( p1c1, p2c1 ) );

    CGAL::Object o;
    CGAL::Oneset_iterator< CGAL::Object > oi( o );

    this->intersectCurves( curve, verticalLine, oi );

    IntersectionResult pair;
    if ( CGAL::assign( pair, o ) )
    {
      Arc_point_2 pt = pair.first;
      res = pt.y( );
    }
    return res;
  }

  double approx( const X_monotone_curve_2& curve, const FT& x )
  {
    return CGAL::to_double( (*this)( curve, x ) );
  }

  // FIXME: inexact projection
  Root_of_2 operator() ( const X_monotone_curve_2& curve, const Root_of_2& x )
  {
    FT approx( CGAL::to_double( x ) );
    return this->operator()( curve, approx );
  }

  double approx( const X_monotone_curve_2& curve, const Root_of_2& x )
  {
    return CGAL::to_double( (*this)( curve, x ) );
  }

protected:
  Traits traits;
  Intersect_2 intersectCurves;
};

template <class Coefficient_>
class Arr_compute_y_at_x_2<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>> {
public:
  typedef Coefficient_ Coefficient;
  typedef CGAL::Arr_algebraic_segment_traits_2< Coefficient > Traits;
  typedef typename Traits::Algebraic_real_1             CoordinateType;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

  CoordinateType operator() ( const X_monotone_curve_2& curve,
                              const CoordinateType& x )
  {
    CGAL::Object o;
    CGAL::Oneset_iterator< CGAL::Object > oi( o );
    Intersect_2 intersect = traits.intersect_2_object( );
    X_monotone_curve_2 c2 = this->makeVerticalLine( x );
    intersect( curve, c2, oi );
    std::pair< Point_2, Multiplicity > res;
    if ( CGAL::assign( res, o ) ) // TODO: handle failure case
    {
      Point_2 p = res.first;
      // std::cout << "approx y: " << p.to_double( ).second << std::endl;
      CoordinateType coord = p.y( );
      return coord;
    }
    else
    {
      std::cout << "Warning: vertical projection failed" << std::endl;
      return CoordinateType( 0 );
    }
  }

  double approx( const X_monotone_curve_2& curve, const CoordinateType& x )
  {
    CGAL::Object o;
    CGAL::Oneset_iterator< CGAL::Object > oi( o );
    Intersect_2 intersect = traits.intersect_2_object( );
    X_monotone_curve_2 c2 = this->makeVerticalLine( x );
    intersect( curve, c2, oi );
    std::pair< Point_2, Multiplicity > res;
    if ( CGAL::assign( res, o ) ) // TODO: handle failure case
    {
      Point_2 p = res.first;
      std::pair< double, double > tmp = p.to_double();
      return tmp.second;
    }
    else
    {
      std::cout << "Warning: vertical projection failed" << std::endl;
      return 0;
    }
  }

protected:
  X_monotone_curve_2 makeVerticalLine( const CoordinateType& x )
  {
    typename Traits::Construct_point_2 constructPoint =
      traits.construct_point_2_object( );
    typename Traits::Construct_x_monotone_segment_2 constructSegment =
      traits.construct_x_monotone_segment_2_object( );

    std::vector< X_monotone_curve_2 > curves;
    Point_2 p1 = constructPoint( x, CoordinateType( -1000000 ) );
    Point_2 p2 = constructPoint( x, CoordinateType( +1000000 ) );
    constructSegment( p1, p2, std::back_inserter( curves ) );
    return curves[ 0 ]; // by construction, there is one curve in curves
  }
  Traits traits;
};

// TODO: Make Construct_x_monotone_subcurve_2 more generic
template < class ArrTraits >
class Construct_x_monotone_subcurve_2
{
public:
  typedef typename ArrTraitsAdaptor<ArrTraits>::Kernel  Kernel;
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename ArrTraits::Split_2                   Split_2;
  typedef typename ArrTraits::Intersect_2               Intersect_2;
  typedef typename ArrTraits::Multiplicity              Multiplicity;
  typedef typename ArrTraits::Construct_min_vertex_2    Construct_min_vertex_2;
  typedef typename ArrTraits::Construct_max_vertex_2    Construct_max_vertex_2;
  typedef typename ArrTraits::Compare_x_2               Compare_x_2;
  typedef typename Kernel::FT                           FT;
  typedef typename ArrTraitsAdaptor< ArrTraits >::CoordinateType
                                                        CoordinateType;
  typedef typename ArrTraits::Point_2                   Point_2;
  typedef typename Kernel::Point_2                      Kernel_point_2;
  //typedef typename Kernel::Line_2 Line_2;
  //typedef typename Kernel::Compute_y_at_x_2 Compute_y_at_x_2;

  Construct_x_monotone_subcurve_2( ):
    intersect_2( this->traits.intersect_2_object( ) ),
    split_2( this->traits.split_2_object( ) ),
    compare_x_2( this->traits.compare_x_2_object( ) ),
    construct_min_vertex_2( this->traits.construct_min_vertex_2_object( ) ),
    construct_max_vertex_2( this->traits.construct_max_vertex_2_object( ) )
  { }

  /*
    Return the subcurve of curve bracketed by pLeft and pRight.

    We assume pLeft and pRight don't lie on the curve and always do a vertical
    projection.
  */
  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const Point_2& pLeft, const Point_2& pRight )
  {
    Point_2 pMin = this->construct_min_vertex_2( curve );
    Point_2 pMax = this->construct_max_vertex_2( curve );
    X_monotone_curve_2 subcurve;
    X_monotone_curve_2 unusedTrimmings;
    X_monotone_curve_2 finalSubcurve;
    if ( this->compare_x_2( pLeft, pMin ) == CGAL::LARGER )
    {
      CoordinateType y1 = this->compute_y_at_x( curve, pLeft.x( ) );
      Point_2 splitPoint( pLeft.x( ), y1 );
      this->split_2( curve, splitPoint, unusedTrimmings, subcurve );
    }
    else
    {
      subcurve = curve;
    }

    if ( this->compare_x_2( pRight, pMax ) == CGAL::SMALLER )
    {
      CoordinateType y2 = this->compute_y_at_x( subcurve, pRight.x( ) );
      Point_2 splitPoint( pRight.x( ), y2 );
      this->split_2( subcurve, splitPoint, finalSubcurve, unusedTrimmings );
    }
    else
    {
      finalSubcurve = subcurve;
    }

    return finalSubcurve;
  }

protected:
  ArrTraits traits;
  Intersect_2 intersect_2;
  Split_2 split_2;
  Compare_x_2 compare_x_2;
  Arr_compute_y_at_x_2< ArrTraits > compute_y_at_x;
  Construct_min_vertex_2 construct_min_vertex_2;
  Construct_max_vertex_2 construct_max_vertex_2;
}; // class Construct_x_monotone_subcurve_2

template < class CircularKernel >
class Construct_x_monotone_subcurve_2< CGAL::Arr_circular_arc_traits_2<
                                         CircularKernel > >
{
public:
  typedef CGAL::Arr_circular_arc_traits_2<CircularKernel> ArrTraits;
  typedef typename ArrTraits::Intersect_2               Intersect_2;
  typedef typename ArrTraits::Split_2                   Split_2;
  typedef typename ArrTraits::Compare_x_2               Compare_x_2;
  typedef typename ArrTraits::Construct_min_vertex_2    Construct_min_vertex_2;
  typedef typename ArrTraits::Construct_max_vertex_2    Construct_max_vertex_2;
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename CircularKernel::Point_2              Non_arc_point_2;
  typedef typename ArrTraits::Point_2                   Arc_point_2;
  typedef typename CircularKernel::FT                   FT;
  typedef typename CircularKernel::Root_of_2            Root_of_2;
  typedef typename CircularKernel::Root_for_circles_2_2 Root_for_circles_2_2;

public:
  Construct_x_monotone_subcurve_2( ):
    intersect_2( this->traits.intersect_2_object( ) ),
    split_2( this->traits.split_2_object( ) ),
    compare_x_2( this->traits.compare_x_2_object( ) ),
    construct_min_vertex_2( this->traits.construct_min_vertex_2_object( ) ),
    construct_max_vertex_2( this->traits.construct_max_vertex_2_object( ) )
  { }

  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const Arc_point_2& pLeft,
                                  const Arc_point_2& pRight )
  {
    Arc_point_2 pMin = this->construct_min_vertex_2( curve );
    Arc_point_2 pMax = this->construct_max_vertex_2( curve );
    X_monotone_curve_2 subcurve;
    X_monotone_curve_2 unusedTrimmings;
    X_monotone_curve_2 finalSubcurve;
    if ( this->compare_x_2( pLeft, pMin ) == CGAL::LARGER )
    {
      Arr_compute_y_at_x_2< ArrTraits > compute_y_at_x;
      FT x_approx( CGAL::to_double( pLeft.x( ) ) );
      Root_of_2 y1 = compute_y_at_x( curve, x_approx );
      Root_for_circles_2_2 intersectionPoint( x_approx, y1 );
      Arc_point_2 splitPoint( intersectionPoint );
      this->split_2( curve, splitPoint, unusedTrimmings, subcurve );
    }
    else
    {
      subcurve = curve;
    }

    if ( this->compare_x_2( pRight, pMax ) == CGAL::SMALLER )
    {
      Arr_compute_y_at_x_2< ArrTraits > compute_y_at_x;
      FT x_approx( CGAL::to_double( pRight.x( ) ) );
      Root_of_2 y2 = compute_y_at_x( subcurve, x_approx );
      Root_for_circles_2_2 intersectionPoint( x_approx, y2 );
      Arc_point_2 splitPoint( intersectionPoint );
      this->split_2( subcurve, splitPoint, finalSubcurve, unusedTrimmings );
    }
    else
    {
      finalSubcurve = subcurve;
    }

    return finalSubcurve;
  }

protected:
  ArrTraits traits;
  Intersect_2 intersect_2;
  Split_2 split_2;
  Compare_x_2 compare_x_2;
  Construct_min_vertex_2 construct_min_vertex_2;
  Construct_max_vertex_2 construct_max_vertex_2;
};

/*
 * This specialization for conic traits makes use of X_monotone_curve_2::trim,
 * which is not necessarily available.
 */
template < class RatKernel, class AlgKernel, class NtTraits >
class Construct_x_monotone_subcurve_2< CGAL::Arr_conic_traits_2< RatKernel,
                                                                 AlgKernel,
                                                                 NtTraits > >
{
public:
  typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > ArrTraits;
  typedef typename ArrTraits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename AlgKernel::Point_2 Point_2;

  /*
    Return the subcurve of curve bracketed by pLeft and pRight.
  */
  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const Point_2& pLeft, const Point_2& pRight )
  {
    // find the points on the curve
    Point_2 left = curve.point_at_x( pLeft );
    Point_2 right = curve.point_at_x( pRight );

    // make sure the points are oriented in the direction that the curve is
    // going
    AlgKernel ker;
    if (! (((curve.is_directed_right( )) &&
            ker.compare_xy_2_object() ( left, right ) == CGAL::SMALLER) ||
           ((! curve.is_directed_right( )) &&
            ker.compare_xy_2_object() ( left, right ) == CGAL::LARGER)))
    {
      Point_2 tmp = left;
      left = right;
      right = tmp;
    }

    X_monotone_curve_2 res = curve.trim( left, right );
    return res;
  }
}; // class Construct_x_monotone_subcurve_2 for Arr_conic_traits_2

template < class Kernel_ >
class Construct_x_monotone_subcurve_2< CGAL::Arr_linear_traits_2< Kernel_ > >
{
public: // typedefs
  typedef CGAL::Arr_linear_traits_2< Kernel_ > ArrTraits;
  typedef typename ArrTraits::X_monotone_curve_2 X_monotone_curve_2;
  typedef Kernel_ Kernel;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;

public: // methods
  // curve can be unbounded. if curve is unbounded to the left,
  // pLeft is a point on the left edge of viewport.
  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const Point_2& pLeft, const Point_2& pRight )
  {
    if ( curve.is_segment( ) )
    {
      Segment_2 subsegment =
        this->constructSubsegment( curve.segment( ), pLeft, pRight );
      return X_monotone_curve_2( subsegment );
    }
    else if ( curve.is_ray( ) )
    {

    }
    return curve;
  }

protected:
  Construct_x_monotone_subcurve_2< CGAL::Arr_segment_traits_2< Kernel_ > >
    constructSubsegment;
};

template < class Coefficient_ >
class Construct_x_monotone_subcurve_2< CGAL::Arr_algebraic_segment_traits_2<
                                         Coefficient_ > >
{
public: // typedefs
  typedef Coefficient_ Coefficient;
  typedef CGAL::Arr_algebraic_segment_traits_2< Coefficient > ArrTraits;
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename ArrTraitsAdaptor< ArrTraits >::Kernel Kernel;
  typedef typename ArrTraits::Point_2                   Point_2;
  //typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2                    Segment_2;

public: // methods
  // curve can be unbounded. if curve is unbounded to the left, pLeft is a
  // point on the left edge of viewport.
  X_monotone_curve_2 operator()(const X_monotone_curve_2& curve,
                                const Point_2& /* pLeft */,
                                const Point_2& /* pRight */)
  {
    // TODO: trim the algebraic curve
    return curve;
  }

protected:
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
  Qt::Converter< Kernel > convert;
}; // class SimpleArrangementViewerQtBase

template <typename Arrangement_2>
class SimpleArrangementViewerQt
    : public SimpleArrangementViewerQtBase<Arrangement_2> {
public:
  SimpleArrangementViewerQt(QWidget *parent)
      : SimpleArrangementViewerQtBase<Arrangement_2>(parent) {}
};

template <typename Kernel_>
class SimpleArrangementViewerQt<
    CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel_>>>
    : public SimpleArrangementViewerQtBase<
          CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel_>>> {
public:
  typedef Kernel_                                       Kernel;
  typedef CGAL::Arr_segment_traits_2<Kernel>            Traits;
  typedef CGAL::Arrangement_2<Traits>                   Arrangement_2;

  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Edge_const_iterator   Edge_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                        Ccb_halfedge_const_circulator;

  typedef Arr_traits_basic_adaptor_2<Traits>            Traits_adapter_2;

  typedef SimpleArrangementViewerQtBase<Arrangement_2>  Superclass;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

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

template <typename CircularKernel>
class SimpleArrangementViewerQt<
    CGAL::Arrangement_2<CGAL::Arr_circle_segment_traits_2<CircularKernel>>>
    : public SimpleArrangementViewerQtBase<CGAL::Arrangement_2<
          CGAL::Arr_circle_segment_traits_2<CircularKernel>>> {
public:
  typedef CircularKernel                                        Kernel;
  typedef CGAL::Arr_circle_segment_traits_2<Kernel>             Traits;

  typedef CGAL::Arrangement_2<Traits>                           Arrangement_2;
  typedef typename Arrangement_2::Halfedge_const_handle         Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle             Face_const_handle;
  typedef typename Arrangement_2::Edge_const_iterator           Edge_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                                Ccb_halfedge_const_circulator;
  typedef Arr_traits_basic_adaptor_2<Traits>                    Traits_adapter_2;
  typedef SimpleArrangementViewerQtBase<Arrangement_2>          Superclass;

  typedef typename Traits::Point_2                              Point_2;
  typedef typename Traits::Curve_2                              Curve_2;
  typedef typename Traits::X_monotone_curve_2                   X_monotone_curve_2;

  typedef typename Kernel::Circle_2 Circle_2;

  typedef CGAL::Exact_predicates_exact_constructions_kernel     Viewer_kernel;
  typedef typename Kernel::FT                                   NT;

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

    Point_2 pSource = ei->source()->point();
    Point_2 pTarget = ei->target()->point();

    std::cout << pSource << std::endl;
    std::cout << pTarget << std::endl;

    Viewer_kernel::Point_3 source(to_double(pSource.x()), 0,
                                  to_double(pSource.y()));
    Viewer_kernel::Point_3 target(to_double(pTarget.x()), 0,
                                  to_double(pTarget.y()));

    if (ei->curve().is_linear()) {
      this->add_segment(source, target);
    } else
    {
      CGAL_precondition(ei->curve().is_circular());
      const Circle_2& circ = ei->curve().supporting_circle();
      Viewer_kernel::Point_3 center(
          to_double(circ.center().x()), 0., to_double(circ.center().y()));
      double radius = to_double(std::sqrt(to_double(circ.squared_radius())));

      // Calculate angles with +x axis
      double start_angle, end_angle;
      double dot = to_double((source.x() - center.x()) * 1. +
                             (source.z() - center.z()) * 0.);
      double det = to_double((source.x() - center.x()) * 0. +
                             (source.z() - center.z()) * 1.);
      start_angle = std::atan2(det, dot);

      dot = to_double((target.x() - center.x()) * 1. +
                      (target.z() - center.z()) * 0.);
      det = to_double((target.x() - center.x()) * 0. +
                      (target.z() - center.z()) * 1.);
      end_angle = std::atan2(det, dot);

      std::cout << "source: \t"<< source << std::endl
           << "target: \t"<< target << std::endl
           << "center: \t"<< center << std::endl
           << "start_angle: \t"<< start_angle << std::endl
           << "end_angle: \t"<< end_angle << std::endl;

      this->add_circle_segment(center, start_angle, end_angle, radius);
    }
    return;
  }

  void compute_elements() {
    this->clear();

    // Draw the arrangement vertices.
    typename Arrangement_2::Vertex_const_iterator vit;
    for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
      Viewer_kernel::Point_3 p(to_double(vit->point().x()), 0,
                               to_double(vit->point().y()));

      this->add_point(p);
    }

    // Draw the arrangement edges.
    typename Arrangement_2::Edge_const_iterator eit;
    for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
      //std::cout << "[" << eit->curve() << "]" << std::endl;
      compute_edge(eit);
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

template < class RatKernel, class AlgKernel, class NtTraits >
class SimpleArrangementViewerQt<
    CGAL::Arrangement_2<CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>>
    : public SimpleArrangementViewerQtBase<CGAL::Arrangement_2<
          CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>> {
public:
  typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > Traits;
  typedef AlgKernel                                             Kernel;

  typedef CGAL::Arrangement_2<Traits>                           Arrangement_2;
  typedef SimpleArrangementViewerQtBase<Arrangement_2>          Superclass;

  typedef typename Arrangement_2::Halfedge_const_handle         Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle             Face_const_handle;
  typedef typename Arrangement_2::Edge_const_iterator           Edge_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                                Ccb_halfedge_const_circulator;
  typedef Arr_traits_basic_adaptor_2<Traits>                    Traits_adapter_2;

  typedef typename Traits::Point_2                              Point_2;
  typedef typename Traits::Curve_2                              Curve_2;
  typedef typename Traits::X_monotone_curve_2                   X_monotone_curve_2;

  typedef typename Traits::Intersect_2 Intersect_2;
  typedef typename Traits::Construct_x_monotone_curve_2         Construct_x_monotone_curve_2;
  typedef typename Traits::Point_2                              Intersection_point_2;
  typedef typename Traits::Multiplicity                         Multiplicity;

  typedef typename Kernel::Circle_2 Circle_2;

  typedef CGAL::Exact_predicates_exact_constructions_kernel     Viewer_kernel;
  typedef typename Kernel::FT                                   NT;

public:
  // utility class to use with std::sort on an Intersect_2 result set.
  class Compare_intersection_point_result {
  public:
    typedef std::pair<Intersection_point_2, Multiplicity> Result;
    // returns whether the point1 < point2, using x-coord to compare
    bool operator()(const Result &o1, const Result &o2) {
      Point_2 p1 = o1.first;
      Point_2 p2 = o2.first;
      return (p1.x() < p2.x());
    }
  };

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
    std::cout << this->bounding_box() << std::endl;
    CGAL::Bbox_2 bb = ei->curve().bbox();
    std::cout << bb << std::endl;

    std::vector<X_monotone_curve_2> visible_parts = visibleParts(ei->curve());

    for(unsigned int i = 0; i < visible_parts.size(); ++i) {
      X_monotone_curve_2 subcurve = visible_parts[i];
      int n = 100;

      std::pair<double, double> *app_pts = new std::pair<double, double>[n + 1];
      std::pair<double, double> *end_pts =
          ei->curve().polyline_approximation(n, app_pts);
      std::pair<double, double> *p_curr = app_pts;
      std::pair<double, double> *p_next = p_curr + 1;
      int count = 0;
      do {
        Viewer_kernel::Point_3 p1(p_curr->first, 0.0, p_curr->second);
        Viewer_kernel::Point_3 p2(p_next->first, 0.0, p_next->second);

        this->add_segment(p1, p2);
        p_curr++;
        p_next++;
        ++count;
      } while (p_next != end_pts);
    }
  }

  std::vector<X_monotone_curve_2> visibleParts(X_monotone_curve_2 curve) {
    const Traits *traits = arr.geometry_traits();
    CGAL::Bbox_3 bbox = this->bounding_box();
    Intersect_2 intersect_2 = traits->intersect_2_object();
    Point_2 bottomLeft(bbox.xmin(), bbox.zmin());
    Point_2 bottomRight(bbox.xmax(), bbox.zmin());
    Point_2 topLeft(bbox.xmin(), bbox.zmax());
    Point_2 topRight(bbox.xmax(), bbox.zmax());
    X_monotone_curve_2 bottom =
        this->construct_x_monotone_curve_2(bottomLeft, bottomRight);
    X_monotone_curve_2 left =
        this->construct_x_monotone_curve_2(bottomLeft, topLeft);
    X_monotone_curve_2 top =
        this->construct_x_monotone_curve_2(topLeft, topRight);
    X_monotone_curve_2 right =
        this->construct_x_monotone_curve_2(topRight, bottomRight);

    std::vector<CGAL::Object> bottomIntersections;
    std::vector<CGAL::Object> leftIntersections;
    std::vector<CGAL::Object> topIntersections;
    std::vector<CGAL::Object> rightIntersections;
    std::vector<CGAL::Object> intersections;

    intersect_2(bottom, curve, std::back_inserter(bottomIntersections));
    intersect_2(left, curve, std::back_inserter(leftIntersections));
    intersect_2(top, curve, std::back_inserter(topIntersections));
    intersect_2(right, curve, std::back_inserter(rightIntersections));

    intersect_2(bottom, curve, std::back_inserter(intersections));
    intersect_2(left, curve, std::back_inserter(intersections));
    intersect_2(top, curve, std::back_inserter(intersections));
    intersect_2(right, curve, std::back_inserter(intersections));

    this->filterIntersectionPoints(intersections);

    Point_2 leftEndpt = curve.source();
    Point_2 rightEndpt = curve.target();
    if (leftEndpt.x() > rightEndpt.x()) {
      std::swap(leftEndpt, rightEndpt);
    }

    std::list<Point_2> pointList;
//    for (unsigned int i = 0; i < intersections.size(); ++i) {
//      CGAL::Object o = intersections[i];
//      std::pair<Intersection_point_2, Multiplicity> pair;
//      if (CGAL::assign(pair, o)) {
//        Point_2 pt = pair.first;
//        pointList.push_back(pt);
//      }
//    }

    CGAL::Bbox_3 leftBox(to_double(leftEndpt.x()), 0.0,
                         to_double(leftEndpt.y()), to_double(leftEndpt.x()),
                         0.0, to_double(leftEndpt.y()));
    CGAL::Bbox_3 rightBox(to_double(rightEndpt.x()), 0.0,
                          to_double(rightEndpt.y()), to_double(leftEndpt.x()),
                          0.0, to_double(leftEndpt.y()));

    bool includeLeftEndpoint = CGAL::do_overlap(leftBox, bbox);
    bool includeRightEndpoint = CGAL::do_overlap(rightBox, bbox);

    if (includeLeftEndpoint)
      pointList.push_front(leftEndpt);
    if (includeRightEndpoint)
      pointList.push_back(rightEndpt);

    Construct_x_monotone_subcurve_2<Traits> construct_x_monotone_subcurve_2;
    std::vector<X_monotone_curve_2> clippings;
    typename std::list<Point_2>::iterator pointListItr = pointList.begin();
    for (unsigned int i = 0; i < pointList.size(); i += 2) {
      Point_2 p1 = *pointListItr++;
      Point_2 p2 = *pointListItr++;
      X_monotone_curve_2 subcurve =
          construct_x_monotone_subcurve_2(curve, p1, p2);
      clippings.push_back(subcurve);
    }
    return clippings;
  }

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

  void compute_elements() {
    this->clear();

    // Draw the arrangement vertices.
    typename Arrangement_2::Vertex_const_iterator vit;
    for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
      Viewer_kernel::Point_3 p(to_double(vit->point().x()), 0,
                               to_double(vit->point().y()));
      this->add_point(p);
    }

    // Draw the arrangement edges.
    typename Arrangement_2::Edge_const_iterator eit;
    for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
      //std::cout << "[" << eit->curve() << "]" << std::endl;
      compute_edge(eit);
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
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2;
};

template <typename Kernel_>
class SimpleArrangementViewerQt<
    CGAL::Arrangement_2<CGAL::Arr_linear_traits_2<Kernel_>>>
    : public SimpleArrangementViewerQtBase<
          CGAL::Arrangement_2<CGAL::Arr_linear_traits_2<Kernel_>>> {
public:
    typedef Kernel_ Kernel;
    typedef CGAL::Arr_linear_traits_2<Kernel> Traits;
    typedef CGAL::Arrangement_2<Traits> Arrangement_2;

    typedef typename Arrangement_2::Halfedge_const_handle         Halfedge_const_handle;
    typedef typename Arrangement_2::Face_const_handle             Face_const_handle;
    typedef typename Arrangement_2::Edge_const_iterator           Edge_const_iterator;
    typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                                  Ccb_halfedge_const_circulator;

    typedef SimpleArrangementViewerQtBase<Arrangement_2>  Superclass;
    typedef typename Superclass::Point_2                  Point_2;
    typedef typename Superclass::Segment_2                Segment_2;
    typedef typename Superclass::Ray_2                    Ray_2;
    typedef typename Superclass::Line_2                   Line_2;
    typedef typename Superclass::Triangle_2               Triangle_2;
    typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
    typedef typename Superclass::Circle_2                 Circle_2;
    typedef typename Traits::Curve_2                      Curve_2;
    typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

    typedef CGAL::Exact_predicates_exact_constructions_kernel Viewer_kernel;
public:
    SimpleArrangementViewerQt(QWidget *parent,
                              const Arrangement_2 &a_arr,
                              const char *title)
        : Superclass(parent, title), arr(a_arr)
    {
      compute_elements();
    }

    void compute_edge(Edge_const_iterator ei) {
      Curve_2 curve = ei->curve();
      if (curve.is_segment()) {
        this->add_segment(curve.source(), curve.target());
      } else if (curve.is_ray()) {
        this->add_ray(curve.source(), curve.supporting_line().to_vector());
      } else if (curve.is_line()) {
        Line_2 line = curve.supporting_line();
        this->add_line(line.point(), line.to_vector());
      }
      return;
    }

    void update_bounding_box(Edge_const_iterator ei) {
      Curve_2 curve = ei->curve();
      if (curve.is_ray()) {
        this->update_bounding_box_for_ray(curve.source(),
                                          curve.supporting_line().to_vector());
      } else if (curve.is_line()) {
        Orientation o = CGAL::POSITIVE;
        Line_2 line = curve.supporting_line();

        this->update_bounding_box_for_line(line.point(), line.to_vector(),
                                           line.to_vector().perpendicular(o));
      }
      return;
    }

    void compute_elements() {
      this->clear();

      // Draw the arrangement vertices.
      typename Arrangement_2::Vertex_const_iterator vit;
      for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
        Viewer_kernel::Point_3 p(to_double(vit->point().x()), 0,
                                 to_double(vit->point().y()));

        this->add_point(p);
      }

      // Update bounding box of viewer to include rays and lines
      typename Arrangement_2::Edge_const_iterator eit;
      for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {

            update_bounding_box(eit);;
      }

      // Draw the arrangement edges
      for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
        compute_edge(eit);
      }

//      // Draw the arrangement faces.
//      typename Arrangement_2::Face_const_iterator fit;
//      for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
//        compute_face(fit);
//      }
      return;
    }
private:
  const Arrangement_2 &arr;
};

template <class GeomTraits_, class TopTraits_>
void draw(const Arrangement_2<GeomTraits_, TopTraits_> &a_arr,
          const char *title = "Basic Arrangement Viewer") {
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
